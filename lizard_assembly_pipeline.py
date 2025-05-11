#!/usr/bin/env python3
"""
Genome Assembly and Evaluation Pipeline for Scincus mitranus (Sandfish Lizard)

Task 2.1: Genome assembly using Hifiasm
Task 2.2: Assembly evaluation through various metrics

Usage:
    python lizard_assembly_pipeline.py --outdir <output_directory> [--threads <num_threads>] [--subsample <fraction>]
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
from pathlib import Path
import multiprocessing
import time

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler("assembly_pipeline.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Path to input data files
INPUT_FILES = {
    "hifi": "lizard_liver_seq.fastq.gz",     # PacBio HiFi reads
    "ont": "lizard_ont.fastq.gz",            # ONT long reads
    "hic_r1": "lizard_hic_R1.fastq.gz",      # Hi-C reads forward
    "hic_r2": "lizard_hic_R2.fastq.gz",      # Hi-C reads reverse
    "rna_eye": "lizard_rna_eye.fastq.gz",    # RNA-Seq eye tissue
    "rna_liver": "lizard_rna_liver.fastq.gz" # RNA-Seq liver tissue
}

# Default parameters
DEFAULT_THREADS = min(32, multiprocessing.cpu_count())
DEFAULT_MEMORY = "64G"
DEFAULT_OUTDIR = "sandfish_assembly"

def run_command(command, description, check=True):
    """Run a shell command and log the output"""
    logger.info(f"Running {description}...")
    logger.info(f"Command: {' '.join(command)}")
    
    try:
        result = subprocess.run(
            command, 
            check=check, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        logger.info(f"{description} completed with exit code {result.returncode}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"{description} failed with exit code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        if not check:
            return e
        raise

def check_dependencies():
    """Check if all required tools are installed"""
    tools = ["hifiasm", "quast.py", "busco", "meryl", "merqury.sh", "Inspector/inspector.py"]
    missing = []
    
    for tool in tools:
        try:
            subprocess.run(["which", tool], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            missing.append(tool)
    
    if missing:
        logger.error(f"Missing dependencies: {', '.join(missing)}")
        logger.error("Please install missing tools before proceeding")
        sys.exit(1)
    else:
        logger.info("All dependencies are installed")

def subsample_reads(input_file, output_file, fraction):
    """Subsample FASTQ reads using seqtk"""
    logger.info(f"Subsampling {input_file} to {fraction} of original reads")
    try:
        command = ["seqtk", "sample", input_file, str(fraction)]
        with open(output_file, 'w') as outfile:
            subprocess.run(command, check=True, stdout=outfile)
        
        # Compress the output file
        subprocess.run(["gzip", output_file], check=True)
        return f"{output_file}.gz"
    except subprocess.CalledProcessError:
        logger.error(f"Failed to subsample {input_file}")
        return input_file

def run_hifiasm(
    hifi_reads, 
    ont_reads, 
    hic_r1, 
    hic_r2, 
    output_prefix, 
    threads
):
    """Run Hifiasm for genome assembly with Hi-C data for phasing"""
    # Estimate haploid genome size based on input files (used for memory estimation)
    # This is a rough estimate for a reptile genome
    hifiasm_cmd = [
        "hifiasm",
        "-o", output_prefix,
        "-t", str(threads),
        "--h1", hic_r1,
        "--h2", hic_r2,
        "--ul", ont_reads,
        "--hg-size", "2g",  # Expected haploid genome size for the lizard
        "--primary",
        hifi_reads
    ]
    
    run_command(hifiasm_cmd, "Hifiasm assembly")
    
    # Wait a moment to ensure files are fully written
    time.sleep(5)
    
    # Convert GFA to FASTA for both haplotypes
    # Primary/maternal assembly
    primary_gfa = f"{output_prefix}.hic.p_ctg.gfa"
    primary_fasta = f"{output_prefix}.hic.p_ctg.fasta"
    if os.path.exists(primary_gfa):
        logger.info("Converting primary assembly GFA to FASTA...")
        gfa_to_fasta(primary_gfa, primary_fasta)
    else:
        logger.error(f"Primary assembly GFA file not found: {primary_gfa}")
        primary_fasta = None
    
    # Alternate/paternal assembly
    alt_gfa = f"{output_prefix}.hic.a_ctg.gfa"
    alt_fasta = f"{output_prefix}.hic.a_ctg.fasta"
    if os.path.exists(alt_gfa):
        logger.info("Converting alternate assembly GFA to FASTA...")
        gfa_to_fasta(alt_gfa, alt_fasta)
    else:
        logger.error(f"Alternate assembly GFA file not found: {alt_gfa}")
        alt_fasta = None
    
    return {
        "primary": primary_fasta,
        "alternate": alt_fasta
    }

def gfa_to_fasta(gfa_file, fasta_file):
    """Convert GFA format to FASTA format"""
    try:
        with open(gfa_file, 'r') as infile, open(fasta_file, 'w') as outfile:
            for line in infile:
                if line.startswith('S'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        seq_id = parts[1]
                        sequence = parts[2]
                        outfile.write(f">{seq_id}\n{sequence}\n")
        logger.info(f"Successfully converted {gfa_file} to {fasta_file}")
    except Exception as e:
        logger.error(f"Error converting GFA to FASTA: {str(e)}")
        raise

def run_quast(assembly_file, output_dir, threads):
    """Run QUAST for basic assembly metrics"""
    os.makedirs(output_dir, exist_ok=True)
    
    quast_cmd = [
        "quast.py",
        "-o", output_dir,
        "-t", str(threads),
        assembly_file
    ]
    
    run_command(quast_cmd, "QUAST analysis")
    
    # Parse and return key metrics
    report_file = os.path.join(output_dir, "report.txt")
    metrics = {}
    
    if os.path.exists(report_file):
        with open(report_file, 'r') as f:
            for line in f:
                line = line.strip()
                if "Total length" in line:
                    metrics["total_length"] = line.split()[-1]
                elif "# contigs" in line and "# contigs (>= 0 bp)" not in line:
                    metrics["num_contigs"] = line.split()[-1]
                elif "GC (%)" in line:
                    metrics["gc_content"] = line.split()[-1]
                elif "N50" in line:
                    metrics["n50"] = line.split()[-1]
                elif "N90" in line:
                    metrics["n90"] = line.split()[-1]
                elif "L50" in line:
                    metrics["l50"] = line.split()[-1]
                elif "Largest contig" in line:
                    metrics["largest_contig"] = line.split()[-1]
    
    return metrics

def run_busco(assembly_file, output_dir, lineage, threads, force=True):
    """Run BUSCO to assess gene completeness"""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    
    # Use vertebrata lineage for Scincus mitranus
    busco_cmd = [
        "busco",
        "-i", assembly_file,
        "-o", os.path.basename(output_dir),
        "-m", "genome",
        "-l", lineage,
        "-c", str(threads),
        "--out_path", os.path.dirname(output_dir)
    ]
    
    # Add force flag to overwrite existing files if needed
    if force:
        busco_cmd.append("-f")
    
    result = run_command(busco_cmd, "BUSCO analysis", check=False)
    
    # Parse and return BUSCO results
    busco_results = {}
    
    # Different versions of BUSCO might put the summary file in different locations
    possible_summary_files = [
        os.path.join(output_dir, "short_summary.txt"),
        os.path.join(output_dir, "short_summary.specific.vertebrata_odb10.txt"),
        os.path.join(output_dir, "run_vertebrata_odb10", "short_summary.txt")
    ]
    
    summary_found = False
    for summary_file in possible_summary_files:
        if os.path.exists(summary_file):
            summary_found = True
            with open(summary_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if "Complete BUSCOs" in line:
                        busco_results["complete"] = line.split("(")[0].strip().split()[-1]
                    elif "Complete and single-copy BUSCOs" in line:
                        busco_results["single_copy"] = line.split("(")[0].strip().split()[-1]
                    elif "Complete and duplicated BUSCOs" in line:
                        busco_results["duplicated"] = line.split("(")[0].strip().split()[-1]
                    elif "Fragmented BUSCOs" in line:
                        busco_results["fragmented"] = line.split("(")[0].strip().split()[-1]
                    elif "Missing BUSCOs" in line:
                        busco_results["missing"] = line.split("(")[0].strip().split()[-1]
            break
    
    if not summary_found:
        logger.warning("BUSCO summary file not found. Check the BUSCO output directory.")
    
    return busco_results

def run_merqury(assembly_file, hifi_reads, output_dir, threads):
    """Run Merqury to evaluate k-mer distribution and QV score"""
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create meryl database from HiFi reads
    meryl_db = os.path.join(output_dir, "hifi_meryl_db")
    meryl_cmd = [
        "meryl",
        "k=31",
        f"threads={threads}",
        "count",
        "output", meryl_db,
        hifi_reads
    ]
    
    run_command(meryl_cmd, "Building Meryl database")
    
    # Run Merqury for QV statistics
    merqury_out_prefix = os.path.join(output_dir, "merqury")
    assembly_name = os.path.basename(assembly_file).split('.')[0]
    
    merqury_cmd = [
        "merqury.sh",
        meryl_db,
        assembly_file,
        merqury_out_prefix
    ]
    
    run_command(merqury_cmd, "Merqury QV analysis", check=False)
    
    # Parse QV results - try different possible filenames
    qv_results = {}
    possible_qv_files = [
        # f"{merqury_out_prefix}.qv",
        f"{merqury_out_prefix}.{assembly_name}.qv"
    ]
    
    qv_found = False
    for qv_file in possible_qv_files:
        if os.path.exists(qv_file):
            qv_found = True
            with open(qv_file, 'r') as f:
                for line in f:
                    if not line.startswith("#"):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            qv_results["qv_score"] = parts[3]
            break
    
    if not qv_found:
        logger.warning("Merqury QV file not found. Check the Merqury output directory.")
        qv_results["qv_score"] = "N/A"
    
    return qv_results

def run_inspector(assembly_file, hifi_reads, ont_reads, output_dir, threads):
    """Run Inspector to identify potential misassemblies"""
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    inspector_cmd = [
        "Inspector/inspector.py",
        "-c", assembly_file,
        "-r", hifi_reads,
        "--datatype hifi",
        "--output", output_dir,
        # "--threads", str(threads)
    ]
    
    run_command(inspector_cmd, "Inspector misassembly analysis", check=False)
    
    # Parse results
    inspector_results = {}
    summary_file = os.path.join(output_dir, "summary_statistics")
    
    if os.path.exists(summary_file):
        with open(summary_file, 'r') as f:
            for line in f:
                line = line.strip()
                if "Total number of misassemblies" in line:
                    inspector_results["misassemblies"] = line.split(":")[-1].strip()
    else:
        logger.warning(f"Inspector summary file not found: {summary_file}")
        inspector_results["misassemblies"] = "N/A"
    
    return inspector_results

def generate_report(results, output_file):
    """Generate a comprehensive report with all evaluation results"""
    with open(output_file, 'w') as f:
        f.write("# Scincus mitranus Genome Assembly Evaluation Report\n\n")
        
        # Basic assembly metrics
        f.write("## QUAST Basic Assembly Metrics\n\n")
        f.write(f"- Total Length: {results.get('quast', {}).get('total_length', 'N/A')}\n")
        f.write(f"- Number of Contigs: {results.get('quast', {}).get('num_contigs', 'N/A')}\n")
        f.write(f"- Largest Contig: {results.get('quast', {}).get('largest_contig', 'N/A')}\n")
        f.write(f"- GC Content: {results.get('quast', {}).get('gc_content', 'N/A')}\n")
        f.write(f"- N50: {results.get('quast', {}).get('n50', 'N/A')}\n")
        f.write(f"- N90: {results.get('quast', {}).get('n90', 'N/A')}\n")
        f.write(f"- L50: {results.get('quast', {}).get('l50', 'N/A')}\n\n")
        
        # BUSCO results
        f.write("## BUSCO Gene Completeness\n\n")
        f.write(f"- Complete BUSCOs: {results.get('busco', {}).get('complete', 'N/A')}\n")
        f.write(f"- Complete Single-Copy BUSCOs: {results.get('busco', {}).get('single_copy', 'N/A')}\n")
        f.write(f"- Complete Duplicated BUSCOs: {results.get('busco', {}).get('duplicated', 'N/A')}\n")
        f.write(f"- Fragmented BUSCOs: {results.get('busco', {}).get('fragmented', 'N/A')}\n")
        f.write(f"- Missing BUSCOs: {results.get('busco', {}).get('missing', 'N/A')}\n\n")
        
        # Merqury results
        f.write("## Merqury QV Score\n\n")
        f.write(f"- QV Score: {results.get('merqury', {}).get('qv_score', 'N/A')}\n\n")
        
        # Inspector results
        f.write("## Inspector Misassembly Analysis\n\n")
        f.write(f"- Total Misassemblies: {results.get('inspector', {}).get('misassemblies', 'N/A')}\n\n")
        
        # Overall assessment
        f.write("## Overall Assessment\n\n")
        
        # Assess if this is a high-quality assembly
        qv_score_str = results.get('merqury', {}).get('qv_score', '0')
        qv_score = float(qv_score_str) if qv_score_str.replace('.', '', 1).isdigit() else 0
        
        complete_buscos_str = results.get('busco', {}).get('complete', '0')
        complete_buscos = int(complete_buscos_str) if complete_buscos_str.isdigit() else 0
        
        misassemblies_str = results.get('inspector', {}).get('misassemblies', '0')
        misassemblies = int(misassemblies_str) if misassemblies_str.isdigit() else 0
        
        if qv_score >= 40 and complete_buscos >= 90 and misassemblies < 100:
            f.write("This assembly is of HIGH QUALITY based on the evaluation metrics.\n")
        elif qv_score >= 30 and complete_buscos >= 80 and misassemblies < 200:
            f.write("This assembly is of MODERATE QUALITY based on the evaluation metrics.\n")
        else:
            f.write("This assembly needs improvement based on the evaluation metrics.\n")
        
        f.write("\nPotential improvements:\n")
        
        if qv_score < 40:
            f.write("- Improve base accuracy through additional sequencing or better error correction\n")
        if complete_buscos < 90:
            f.write("- Improve gene completeness by addressing missing regions or fragmented genes\n")
        if misassemblies >= 100:
            f.write("- Address structural issues to reduce the number of misassemblies\n")

def main():
    parser = argparse.ArgumentParser(description='Genome Assembly and Evaluation Pipeline')
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR, help=f'Output directory (default: {DEFAULT_OUTDIR})')
    parser.add_argument('--threads', type=int, default=DEFAULT_THREADS, help=f'Number of threads (default: {DEFAULT_THREADS})')
    parser.add_argument('--subsample', type=float, help='Subsample reads to this fraction (e.g., 0.1 for 10%)')
    parser.add_argument('--skip-assembly', action='store_true', help='Skip assembly step (use existing assembly)')
    parser.add_argument('--input-dir', default=".", help='Directory containing input files')
    parser.add_argument('--assembly-path', help='Path to existing assembly (if skipping assembly step)')
    parser.add_argument('--force', action='store_true', help='Force overwrite of existing files')
    
    args = parser.parse_args()
    
    # Start timestamp for the run
    start_time = time.time()
    logger.info(f"Starting assembly pipeline at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Make sure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "data"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "assembly"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "evaluation"), exist_ok=True)
    
    # Check dependencies
    try:
        check_dependencies()
    except Exception as e:
        logger.error(f"Failed to check dependencies: {str(e)}")
        sys.exit(1)
    
    # Prepare input files (with optional subsampling)
    input_paths = {}
    for key, filename in INPUT_FILES.items():
        input_path = os.path.join(args.input_dir, filename)
        
        if not os.path.exists(input_path):
            logger.error(f"Input file not found: {input_path}")
            sys.exit(1)
        
        if args.subsample and key in ["hifi", "ont"]:
            output_file = os.path.join(args.outdir, "data", f"subsampled_{filename.replace('.gz', '')}")
            input_paths[key] = subsample_reads(input_path, output_file, args.subsample)
        else:
            input_paths[key] = input_path
    
    assembly_files = {}
    
    # Run assembly pipeline unless --skip-assembly is specified
    if not args.skip_assembly:
        logger.info("Starting genome assembly with Hifiasm...")
        
        try:
            assembly_files = run_hifiasm(
                input_paths["hifi"],
                input_paths["ont"],
                input_paths["hic_r1"],
                input_paths["hic_r2"],
                os.path.join(args.outdir, "assembly", "sandfish"),
                args.threads
            )
        except Exception as e:
            logger.error(f"Assembly with Hifiasm failed: {str(e)}")
            sys.exit(1)
    else:
        if not args.assembly_path:
            logger.error("--assembly-path must be specified when using --skip-assembly")
            sys.exit(1)
        
        assembly_files["primary"] = args.assembly_path
        logger.info(f"Skipping assembly step, using existing assembly: {assembly_files['primary']}")
    
    # For evaluation, we'll use the primary assembly
    assembly_to_evaluate = assembly_files.get("primary")
    
    if not assembly_to_evaluate or not os.path.exists(assembly_to_evaluate):
        logger.error(f"Assembly file not found: {assembly_to_evaluate}")
        sys.exit(1)
    
    # Create evaluation directories
    quast_dir = os.path.join(args.outdir, "evaluation", "quast")
    busco_dir = os.path.join(args.outdir, "evaluation", "busco")
    merqury_dir = os.path.join(args.outdir, "evaluation", "merqury")
    inspector_dir = os.path.join(args.outdir, "evaluation", "inspector")
    
    # Run evaluation tools
    logger.info("Starting assembly evaluation...")
    
    results = {}
    
    # QUAST
    try:
        logger.info("Running QUAST for basic metrics...")
        results["quast"] = run_quast(assembly_to_evaluate, quast_dir, args.threads)
    except Exception as e:
        logger.error(f"QUAST analysis failed: {str(e)}")
        results["quast"] = {}
    
    # BUSCO
    try:
        logger.info("Running BUSCO for gene completeness...")
        results["busco"] = run_busco(
            assembly_to_evaluate, 
            busco_dir, 
            "vertebrata_odb10", 
            args.threads, 
            force=args.force
        )
    except Exception as e:
        logger.error(f"BUSCO analysis failed: {str(e)}")
        results["busco"] = {}
    
    # Merqury
    try:
        logger.info("Running Merqury for k-mer distribution and QV score...")
        results["merqury"] = run_merqury(assembly_to_evaluate, input_paths["hifi"], merqury_dir, args.threads)
    except Exception as e:
        logger.error(f"Merqury analysis failed: {str(e)}")
        results["merqury"] = {}
    
    # Inspector
    try:
        logger.info("Running Inspector to identify misassemblies...")
        results["inspector"] = run_inspector(
            assembly_to_evaluate, 
            input_paths["hifi"],
            input_paths["ont"],
            inspector_dir,
            args.threads
        )
    except Exception as e:
        logger.error(f"Inspector analysis failed: {str(e)}")
        results["inspector"] = {}
    
    # Generate comprehensive report
    report_file = os.path.join(args.outdir, "assembly_evaluation_report.md")
    generate_report(results, report_file)
    logger.info(f"Assembly evaluation complete. Report written to {report_file}")
    
    # Calculate and log total runtime
    end_time = time.time()
    runtime_hours = (end_time - start_time) / 3600
    logger.info(f"Total runtime: {runtime_hours:.2f} hours")
    
    # Output summary
    logger.info("\n=== Assembly Evaluation Summary ===")
    logger.info(f"Assembly size: {results.get('quast', {}).get('total_length', 'N/A')}")
    logger.info(f"Number of contigs: {results.get('quast', {}).get('num_contigs', 'N/A')}")
    logger.info(f"N50: {results.get('quast', {}).get('n50', 'N/A')}")
    logger.info(f"BUSCO completeness: {results.get('busco', {}).get('complete', 'N/A')}")
    logger.info(f"QV score: {results.get('merqury', {}).get('qv_score', 'N/A')}")
    logger.info(f"Detected misassemblies: {results.get('inspector', {}).get('misassemblies', 'N/A')}")
    logger.info(f"Full report available at: {report_file}")

if __name__ == "__main__":
    main()

# #!/usr/bin/env python3
# """
# Genome Assembly and Evaluation Pipeline for Scincus mitranus (Sandfish Lizard)

# Task 2.1: Genome assembly using Hifiasm
# Task 2.2: Assembly evaluation through various metrics

# Usage:
#     python lizard_assembly_pipeline.py --outdir <output_directory> [--threads <num_threads>] [--subsample <fraction>]
# """

# import os
# import sys
# import argparse
# import subprocess
# import logging
# import shutil
# from pathlib import Path
# import multiprocessing

# # Setup logging
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s [%(levelname)s] %(message)s',
#     handlers=[
#         logging.FileHandler("assembly_pipeline.log"),
#         logging.StreamHandler(sys.stdout)
#     ]
# )
# logger = logging.getLogger(__name__)

# # Path to input data files
# INPUT_FILES = {
#     "hifi": "lizard_liver_seq.fastq.gz",     # PacBio HiFi reads
#     "ont": "lizard_ont.fastq.gz",            # ONT long reads
#     "hic_r1": "lizard_hic_R1.fastq.gz",      # Hi-C reads forward
#     "hic_r2": "lizard_hic_R2.fastq.gz",      # Hi-C reads reverse
#     "rna_eye": "lizard_rna_eye.fastq.gz",    # RNA-Seq eye tissue
#     "rna_liver": "lizard_rna_liver.fastq.gz" # RNA-Seq liver tissue
# }

# # Default parameters
# DEFAULT_THREADS = min(32, multiprocessing.cpu_count())
# DEFAULT_MEMORY = "64G"
# DEFAULT_OUTDIR = "sandfish_assembly"

# def run_command(command, description, check=True):
#     """Run a shell command and log the output"""
#     logger.info(f"Running {description}...")
#     logger.info(f"Command: {' '.join(command)}")
    
#     try:
#         result = subprocess.run(
#             command, 
#             check=check, 
#             stdout=subprocess.PIPE, 
#             stderr=subprocess.PIPE,
#             universal_newlines=True
#         )
#         logger.info(f"{description} completed with exit code {result.returncode}")
#         return result
#     except subprocess.CalledProcessError as e:
#         logger.error(f"{description} failed with exit code {e.returncode}")
#         logger.error(f"STDOUT: {e.stdout}")
#         logger.error(f"STDERR: {e.stderr}")
#         raise

# def check_dependencies():
#     """Check if all required tools are installed"""
#     tools = ["hifiasm", "quast.py", "busco", "meryl", "merqury/merqury.sh", "Inspector/inspector.py"]
#     missing = []
    
#     for tool in tools:
#         try:
#             subprocess.run(["which", tool], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         except subprocess.CalledProcessError:
#             missing.append(tool)
    
#     if missing:
#         logger.error(f"Missing dependencies: {', '.join(missing)}")
#         logger.error("Please install missing tools before proceeding")
#         sys.exit(1)
#     else:
#         logger.info("All dependencies are installed")

# def subsample_reads(input_file, output_file, fraction):
#     """Subsample FASTQ reads using seqtk"""
#     logger.info(f"Subsampling {input_file} to {fraction} of original reads")
#     try:
#         command = ["seqtk", "sample", input_file, str(fraction)]
#         with open(output_file, 'w') as outfile:
#             subprocess.run(command, check=True, stdout=outfile)
        
#         # Compress the output file
#         subprocess.run(["gzip", output_file], check=True)
#         return f"{output_file}.gz"
#     except subprocess.CalledProcessError:
#         logger.error(f"Failed to subsample {input_file}")
#         return input_file

# def run_hifiasm(
#     hifi_reads, 
#     ont_reads, 
#     hic_r1, 
#     hic_r2, 
#     output_prefix, 
#     threads
# ):
#     """Run Hifiasm for genome assembly with Hi-C data for phasing"""
#     hifiasm_cmd = [
#         "hifiasm",
#         "-o", output_prefix,
#         "-t", str(threads),
#         "--h1", hic_r1,
#         "--h2", hic_r2,
#         "--ul", ont_reads,
#         "--hg-size", "400g",  # Expected haploid genome size for the lizard
#         "--primary",
#         hifi_reads
#     ]
    
#     run_command(hifiasm_cmd, "Hifiasm assembly")
    
#     # Convert GFA to FASTA for both haplotypes
#     # Primary/maternal assembly
#     logger.info("Converting primary assembly GFA to FASTA...")
#     gfa_to_fasta(f"{output_prefix}.hic.p_ctg.gfa", f"{output_prefix}.hic.p_ctg.fasta")
    
#     # Alternate/paternal assembly
#     logger.info("Converting alternate assembly GFA to FASTA...")
#     gfa_to_fasta(f"{output_prefix}.hic.a_ctg.gfa", f"{output_prefix}.hic.a_ctg.fasta")
    
#     return {
#         "primary": f"{output_prefix}.hic.p_ctg.fasta",
#         "alternate": f"{output_prefix}.hic.a_ctg.fasta"
#     }

# def gfa_to_fasta(gfa_file, fasta_file):
#     """Convert GFA format to FASTA format"""
#     with open(gfa_file, 'r') as infile, open(fasta_file, 'w') as outfile:
#         for line in infile:
#             if line.startswith('S'):
#                 parts = line.strip().split('\t')
#                 if len(parts) >= 3:
#                     seq_id = parts[1]
#                     sequence = parts[2]
#                     outfile.write(f">{seq_id}\n{sequence}\n")

# def run_quast(assembly_file, output_dir, threads):
#     """Run QUAST for basic assembly metrics"""
#     quast_cmd = [
#         "quast.py",
#         "-o", output_dir,
#         "-t", str(threads),
#         assembly_file
#     ]
    
#     run_command(quast_cmd, "QUAST analysis")
    
#     # Parse and return key metrics
#     report_file = os.path.join(output_dir, "report.txt")
#     metrics = {}
    
#     if os.path.exists(report_file):
#         with open(report_file, 'r') as f:
#             for line in f:
#                 line = line.strip()
#                 if "Total length" in line:
#                     metrics["total_length"] = line.split()[-1]
#                 elif "# contigs" in line and "# contigs (>= 0 bp)" not in line:
#                     metrics["num_contigs"] = line.split()[-1]
#                 elif "GC (%)" in line:
#                     metrics["gc_content"] = line.split()[-1]
#                 elif "N50" in line:
#                     metrics["n50"] = line.split()[-1]
#                 elif "N90" in line:
#                     metrics["n90"] = line.split()[-1]
#                 elif "L50" in line:
#                     metrics["l50"] = line.split()[-1]
#                 elif "Largest contig" in line:
#                     metrics["largest_contig"] = line.split()[-1]
    
#     return metrics

# def run_busco(assembly_file, output_dir, lineage, threads):
#     """Run BUSCO to assess gene completeness"""
#     # Use vertebrata lineage for Scincus mitranus
#     busco_cmd = [
#         "busco",
#         "-i", assembly_file,
#         "-o", os.path.basename(output_dir),
#         "-m", "genome",
#         "-l", lineage,
#         "-c", str(threads),
#         "--out_path", os.path.dirname(output_dir)
#     ]
    
#     run_command(busco_cmd, "BUSCO analysis")
    
#     # Parse and return BUSCO results
#     busco_results = {}
#     summary_file = os.path.join(output_dir, "short_summary.txt")
    
#     if os.path.exists(summary_file):
#         with open(summary_file, 'r') as f:
#             for line in f:
#                 line = line.strip()
#                 if "Complete BUSCOs" in line:
#                     busco_results["complete"] = line.split("(")[0].strip().split()[-1]
#                 elif "Complete and single-copy BUSCOs" in line:
#                     busco_results["single_copy"] = line.split("(")[0].strip().split()[-1]
#                 elif "Complete and duplicated BUSCOs" in line:
#                     busco_results["duplicated"] = line.split("(")[0].strip().split()[-1]
#                 elif "Fragmented BUSCOs" in line:
#                     busco_results["fragmented"] = line.split("(")[0].strip().split()[-1]
#                 elif "Missing BUSCOs" in line:
#                     busco_results["missing"] = line.split("(")[0].strip().split()[-1]
    
#     return busco_results

# def run_merqury(assembly_file, hifi_reads, output_dir, threads):
#     """Run Merqury to evaluate k-mer distribution and QV score"""
#     # Create meryl database from HiFi reads
#     meryl_db = os.path.join(output_dir, "hifi_meryl_db")
#     meryl_cmd = [
#         "meryl",
#         "k=21",
#         "threads=" + str(threads),
#         "count",
#         "output", meryl_db,
#         hifi_reads
#     ]
    
#     run_command(meryl_cmd, "Building Meryl database")
    
#     # Run Merqury for QV statistics
#     merqury_out_prefix = os.path.join(output_dir, "merqury")
#     merqury_cmd = [
#         "merqury/merqury.sh",
#         meryl_db,
#         assembly_file,
#         merqury_out_prefix
#     ]
    
#     run_command(merqury_cmd, "Merqury QV analysis")
    
#     # Parse QV results
#     qv_results = {}
#     qv_file = f"{merqury_out_prefix}.qv"
    
#     if os.path.exists(qv_file):
#         with open(qv_file, 'r') as f:
#             for line in f:
#                 if not line.startswith("#"):
#                     parts = line.strip().split()
#                     if len(parts) >= 2:
#                         qv_results["qv_score"] = parts[1]
    
#     return qv_results

# def run_inspector(assembly_file, hifi_reads, ont_reads, output_dir, threads):
#     """Run Inspector to identify potential misassemblies"""
#     inspector_cmd = [
#         "Inspector/inspector.py",
#         "--assembly", assembly_file,
#         "--pacbio", hifi_reads,
#         "--ont", ont_reads,
#         "--output", output_dir,
#         "--threads", str(threads)
#     ]
    
#     run_command(inspector_cmd, "Inspector misassembly analysis")
    
#     # Parse results
#     inspector_results = {}
#     summary_file = os.path.join(output_dir, "summary.txt")
    
#     if os.path.exists(summary_file):
#         with open(summary_file, 'r') as f:
#             for line in f:
#                 line = line.strip()
#                 if "Total number of misassemblies" in line:
#                     inspector_results["misassemblies"] = line.split(":")[-1].strip()
    
#     return inspector_results

# def generate_report(results, output_file):
#     """Generate a comprehensive report with all evaluation results"""
#     with open(output_file, 'w') as f:
#         f.write("# Scincus mitranus Genome Assembly Evaluation Report\n\n")
        
#         # Basic assembly metrics
#         f.write("## QUAST Basic Assembly Metrics\n\n")
#         f.write(f"- Total Length: {results.get('quast', {}).get('total_length', 'N/A')}\n")
#         f.write(f"- Number of Contigs: {results.get('quast', {}).get('num_contigs', 'N/A')}\n")
#         f.write(f"- Largest Contig: {results.get('quast', {}).get('largest_contig', 'N/A')}\n")
#         f.write(f"- GC Content: {results.get('quast', {}).get('gc_content', 'N/A')}\n")
#         f.write(f"- N50: {results.get('quast', {}).get('n50', 'N/A')}\n")
#         f.write(f"- N90: {results.get('quast', {}).get('n90', 'N/A')}\n")
#         f.write(f"- L50: {results.get('quast', {}).get('l50', 'N/A')}\n\n")
        
#         # BUSCO results
#         f.write("## BUSCO Gene Completeness\n\n")
#         f.write(f"- Complete BUSCOs: {results.get('busco', {}).get('complete', 'N/A')}\n")
#         f.write(f"- Complete Single-Copy BUSCOs: {results.get('busco', {}).get('single_copy', 'N/A')}\n")
#         f.write(f"- Complete Duplicated BUSCOs: {results.get('busco', {}).get('duplicated', 'N/A')}\n")
#         f.write(f"- Fragmented BUSCOs: {results.get('busco', {}).get('fragmented', 'N/A')}\n")
#         f.write(f"- Missing BUSCOs: {results.get('busco', {}).get('missing', 'N/A')}\n\n")
        
#         # Merqury results
#         f.write("## Merqury QV Score\n\n")
#         f.write(f"- QV Score: {results.get('merqury', {}).get('qv_score', 'N/A')}\n\n")
        
#         # Inspector results
#         f.write("## Inspector Misassembly Analysis\n\n")
#         f.write(f"- Total Misassemblies: {results.get('inspector', {}).get('misassemblies', 'N/A')}\n\n")
        
#         # Overall assessment
#         f.write("## Overall Assessment\n\n")
        
#         # Assess if this is a high-quality assembly
#         qv_score = float(results.get('merqury', {}).get('qv_score', '0')) if results.get('merqury', {}).get('qv_score', '0').replace('.', '', 1).isdigit() else 0
#         complete_buscos = int(results.get('busco', {}).get('complete', '0')) if results.get('busco', {}).get('complete', '0').isdigit() else 0
#         misassemblies = int(results.get('inspector', {}).get('misassemblies', '0')) if results.get('inspector', {}).get('misassemblies', '0').isdigit() else 0
        
#         if qv_score >= 40 and complete_buscos >= 90 and misassemblies < 100:
#             f.write("This assembly is of HIGH QUALITY based on the evaluation metrics.\n")
#         elif qv_score >= 30 and complete_buscos >= 80 and misassemblies < 200:
#             f.write("This assembly is of MODERATE QUALITY based on the evaluation metrics.\n")
#         else:
#             f.write("This assembly needs improvement based on the evaluation metrics.\n")
        
#         f.write("\nPotential improvements:\n")
        
#         if qv_score < 40:
#             f.write("- Improve base accuracy through additional sequencing or better error correction\n")
#         if complete_buscos < 90:
#             f.write("- Improve gene completeness by addressing missing regions or fragmented genes\n")
#         if misassemblies >= 100:
#             f.write("- Address structural issues to reduce the number of misassemblies\n")

# def main():
#     parser = argparse.ArgumentParser(description='Genome Assembly and Evaluation Pipeline')
#     parser.add_argument('--outdir', default=DEFAULT_OUTDIR, help=f'Output directory (default: {DEFAULT_OUTDIR})')
#     parser.add_argument('--threads', type=int, default=DEFAULT_THREADS, help=f'Number of threads (default: {DEFAULT_THREADS})')
#     parser.add_argument('--subsample', type=float, help='Subsample reads to this fraction (e.g., 0.1 for 10%)')
#     parser.add_argument('--skip-assembly', action='store_true', help='Skip assembly step (use existing assembly)')
#     parser.add_argument('--input-dir', default=".", help='Directory containing input files')
#     parser.add_argument('--assembly-path', help='Path to existing assembly (if skipping assembly step)')
    
#     args = parser.parse_args()
    
#     # Make sure output directory exists
#     os.makedirs(args.outdir, exist_ok=True)
#     os.makedirs(os.path.join(args.outdir, "data"), exist_ok=True)
#     os.makedirs(os.path.join(args.outdir, "assembly"), exist_ok=True)
#     os.makedirs(os.path.join(args.outdir, "evaluation"), exist_ok=True)
    
#     # Check dependencies
#     check_dependencies()
    
#     # Prepare input files (with optional subsampling)
#     input_paths = {}
#     for key, filename in INPUT_FILES.items():
#         input_path = os.path.join(args.input_dir, filename)
        
#         if not os.path.exists(input_path):
#             logger.error(f"Input file not found: {input_path}")
#             sys.exit(1)
        
#         if args.subsample and key in ["hifi", "ont"]:
#             output_file = os.path.join(args.outdir, "data", f"subsampled_{filename.replace('.gz', '')}")
#             input_paths[key] = subsample_reads(input_path, output_file, args.subsample)
#         else:
#             input_paths[key] = input_path
    
#     assembly_files = {}
    
#     # Run assembly pipeline unless --skip-assembly is specified
#     if not args.skip_assembly:
#         logger.info("Starting genome assembly with Hifiasm...")
        
#         assembly_files = run_hifiasm(
#             input_paths["hifi"],
#             input_paths["ont"],
#             input_paths["hic_r1"],
#             input_paths["hic_r2"],
#             os.path.join(args.outdir, "assembly", "sandfish"),
#             args.threads
#         )
#     else:
#         if not args.assembly_path:
#             logger.error("--assembly-path must be specified when using --skip-assembly")
#             sys.exit(1)
        
#         assembly_files["primary"] = args.assembly_path
#         logger.info(f"Skipping assembly step, using existing assembly: {assembly_files['primary']}")
    
#     # For evaluation, we'll use the primary assembly
#     assembly_to_evaluate = assembly_files["primary"]
    
#     # Create evaluation directories
#     quast_dir = os.path.join(args.outdir, "evaluation", "quast")
#     busco_dir = os.path.join(args.outdir, "evaluation", "busco")
#     merqury_dir = os.path.join(args.outdir, "evaluation", "merqury")
#     inspector_dir = os.path.join(args.outdir, "evaluation", "inspector")
    
#     os.makedirs(quast_dir, exist_ok=True)
#     os.makedirs(busco_dir, exist_ok=True)
#     os.makedirs(merqury_dir, exist_ok=True)
#     os.makedirs(inspector_dir, exist_ok=True)
    
#     # Run evaluation tools
#     logger.info("Starting assembly evaluation...")
    
#     results = {}
    
#     # QUAST
#     logger.info("Running QUAST for basic metrics...")
#     results["quast"] = run_quast(assembly_to_evaluate, quast_dir, args.threads)
    
#     # BUSCO
#     logger.info("Running BUSCO for gene completeness...")
#     results["busco"] = run_busco(assembly_to_evaluate, busco_dir, "vertebrata_odb10", args.threads)
    
#     # Merqury
#     logger.info("Running Merqury for k-mer distribution and QV score...")
#     results["merqury"] = run_merqury(assembly_to_evaluate, input_paths["hifi"], merqury_dir, args.threads)
    
#     # Inspector
#     logger.info("Running Inspector to identify misassemblies...")
#     results["inspector"] = run_inspector(
#         assembly_to_evaluate, 
#         input_paths["hifi"],
#         input_paths["ont"],
#         inspector_dir,
#         args.threads
#     )
    
#     # Generate comprehensive report
#     report_file = os.path.join(args.outdir, "assembly_evaluation_report.md")
#     generate_report(results, report_file)
#     logger.info(f"Assembly evaluation complete. Report written to {report_file}")
    
#     # Output summary
#     logger.info("\n=== Assembly Evaluation Summary ===")
#     logger.info(f"Assembly size: {results.get('quast', {}).get('total_length', 'N/A')}")
#     logger.info(f"Number of contigs: {results.get('quast', {}).get('num_contigs', 'N/A')}")
#     logger.info(f"N50: {results.get('quast', {}).get('n50', 'N/A')}")
#     logger.info(f"BUSCO completeness: {results.get('busco', {}).get('complete', 'N/A')}")
#     logger.info(f"QV score: {results.get('merqury', {}).get('qv_score', 'N/A')}")
#     logger.info(f"Detected misassemblies: {results.get('inspector', {}).get('misassemblies', 'N/A')}")
#     logger.info(f"Full report available at: {report_file}")

# if __name__ == "__main__":
#     main()