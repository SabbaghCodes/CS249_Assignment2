# olc_assembly/main.py
"""
Main script for Overlap-Layout-Consensus assembly.
"""
import argparse
import sys
import os

from fastq_parser import FastqParser
from overlap_detector import OverlapDetector
from overlap_graph import OverlapGraph
from layout_generator import LayoutGenerator
from consensus_builder import ConsensusBuilder
from contig_generator import ContigGenerator
from fasta_writer import FastaWriter
from assembly_metrics import AssemblyMetrics
from visualize import Visualizer

def run_assembly(fastq_file, min_overlap, output_fasta, max_error_rate=0.05, 
                min_contig_length=None, use_quality=False, verbose=False):
    """
    Run the Overlap-Layout-Consensus assembly pipeline.
    
    Args:
        fastq_file (str): Path to the input FASTQ file.
        min_overlap (int): Minimum overlap length.
        output_fasta (str): Path to the output FASTA file.
        max_error_rate (float): Maximum error rate in overlaps.
        min_contig_length (int, optional): Minimum contig length to keep.
        use_quality (bool): Whether to use quality scores in consensus building.
        verbose (bool): Whether to print verbose output.
        
    Returns:
        tuple: (contigs, metrics) tuple.
    """
    if verbose:
        print(f"Starting OLC assembly with min overlap={min_overlap}, max error rate={max_error_rate}")
        print(f"Reading sequences from {fastq_file}...")
    
    # Parse FASTQ file
    parser = FastqParser(fastq_file)
    reads = parser.parse_as_dict()
    
    if verbose:
        print(f"Read {len(reads)} sequences")
    
    # Compute read overlaps
    if verbose:
        print("Computing read overlaps...")
    
    overlap_detector = OverlapDetector(min_overlap, max_error_rate)
    overlaps = overlap_detector.compute_all_overlaps(reads)
    
    if verbose:
        print(f"Found {len(overlaps)} overlaps")
    
    # Filter contained reads
    if verbose:
        print("Filtering contained reads...")
    
    contained_reads = overlap_detector.filter_contained_reads(reads, overlaps)
    
    if verbose:
        print(f"Filtered out {len(contained_reads)} contained reads")
    
    # Get dovetail overlaps
    dovetail_overlaps = overlap_detector.get_dovetail_overlaps(reads, overlaps, contained_reads)
    
    if verbose:
        print(f"Found {len(dovetail_overlaps)} dovetail overlaps")
    
    # Build overlap graph
    if verbose:
        print("Building overlap graph...")
    
    graph = OverlapGraph()
    # Filter reads to exclude contained reads
    filtered_reads = {read_id: read for read_id, read in reads.items() 
                     if read_id not in contained_reads}
    graph.build_from_overlaps(filtered_reads, dovetail_overlaps)
    
    if verbose:
        print(f"Graph built with {len(graph.graph.nodes())} nodes and {len(graph.graph.edges())} edges")
    
    # Remove transitive edges
    if verbose:
        print("Removing transitive edges...")
    
    graph.remove_transitive_edges()
    
    if verbose:
        print(f"Graph after transitive reduction: {len(graph.graph.nodes())} nodes and {len(graph.graph.edges())} edges")
    
    # Find paths through the graph
    if verbose:
        print("Finding paths through the graph...")
    
    paths = graph.find_nonbranching_paths()
    
    if verbose:
        print(f"Found {len(paths)} paths")
    
    # Generate contigs
    if verbose:
        print("Generating contigs...")
    
    contig_generator = ContigGenerator(filtered_reads, dovetail_overlaps)
    
    if use_quality:
        contigs = contig_generator.generate_contigs_with_quality(paths)
    else:
        contigs = contig_generator.generate_contigs(paths)
    
    # Filter contigs by length if requested
    if min_contig_length:
        original_count = len(contigs)
        contigs = contig_generator.filter_contigs(contigs, min_contig_length)
        if verbose:
            print(f"Filtered contigs by minimum length {min_contig_length}: {original_count} -> {len(contigs)}")
    
    # Write contigs to FASTA file
    if verbose:
        print(f"Writing {len(contigs)} contigs to {output_fasta}...")
    
    fasta_writer = FastaWriter(output_fasta)
    fasta_writer.write_contigs(contigs)
    
    # Calculate assembly metrics
    if verbose:
        print("Calculating assembly metrics...")
    
    metrics = AssemblyMetrics.calculate_basic_metrics(contigs)
    
    if verbose:
        print("\nAssembly Metrics:")
        for metric, value in metrics.items():
            print(f"  {metric}: {value}")
    
    return contigs, metrics


def main():
    """
    Main entry point for the command-line application.
    """
    parser = argparse.ArgumentParser(description="Overlap-Layout-Consensus Genome Assembly")
    
    parser.add_argument("--fastq", "-f", required=True, help="Input FASTQ file")
    parser.add_argument("--min-overlap", "-m", type=int, required=True, help="Minimum overlap length")
    parser.add_argument("--output", "-o", required=True, help="Output FASTA file")
    parser.add_argument("--error-rate", "-e", type=float, default=0.05, help="Maximum error rate in overlaps")
    parser.add_argument("--min-length", "-l", type=int, help="Minimum contig length")
    parser.add_argument("--use-quality", "-q", action="store_true", help="Use quality scores in consensus building")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    try:
        run_assembly(
            args.fastq,
            args.min_overlap,
            args.output,
            args.error_rate,
            args.min_length,
            args.use_quality,
            args.verbose
        )
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())