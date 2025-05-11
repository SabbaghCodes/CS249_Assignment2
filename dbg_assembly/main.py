# dbg_assembly/main.py
"""
Main script for De Bruijn Graph assembly.
"""
import argparse
import sys
import os

from fastq_parser import FastqParser
from kmer_generator import KmerGenerator
from de_bruijn_graph import DeBruijnGraph
from path_finder import PathFinder
from contig_generator import ContigGenerator
from fasta_writer import FastaWriter
from assembly_metrics import AssemblyMetrics
from gfa_export import GfaExporter
from visualize import Visualizer

def run_assembly(fastq_file, k, output_fasta, output_gfa=None, min_contig_length=None, verbose=False):
    """
    Run the De Bruijn Graph assembly pipeline.
    
    Args:
        fastq_file (str): Path to the input FASTQ file.
        k (int): K-mer length.
        output_fasta (str): Path to the output FASTA file.
        output_gfa (str, optional): Path to the output GFA file.
        min_contig_length (int, optional): Minimum contig length to keep.
        verbose (bool): Whether to print verbose output.
        
    Returns:
        tuple: (contigs, metrics) tuple.
    """
    if verbose:
        print(f"Starting De Bruijn Graph assembly with k={k}")
        print(f"Reading sequences from {fastq_file}...")
    
    # Parse FASTQ file
    parser = FastqParser(fastq_file)
    sequences = list(parser.parse())
    
    if verbose:
        print(f"Read {len(sequences)} sequences")
    
    # Create and build the De Bruijn Graph
    if verbose:
        print("Building De Bruijn Graph...")
    
    dbg = DeBruijnGraph(k)
    for seq in sequences:
        dbg.add_sequence(seq)
    
    if verbose:
        print(f"De Bruijn Graph built with {len(dbg.graph.nodes())} nodes and {len(dbg.graph.edges())} edges")
    
    # Find contigs in the graph
    if verbose:
        print("Finding contigs...")
    
    contig_paths = PathFinder.find_contigs(dbg)
    
    if verbose:
        print(f"Found {len(contig_paths)} contig paths")
    
    # Generate contig sequences from paths
    if verbose:
        print("Generating contig sequences...")
    
    contig_generator = ContigGenerator()
    contigs = contig_generator.generate_contigs_from_paths(contig_paths, k)
    
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
    
    # Export graph to GFA if requested
    if output_gfa:
        if verbose:
            print(f"Exporting graph to GFA format: {output_gfa}")
        
        gfa_exporter = GfaExporter(output_gfa)
        gfa_exporter.export_graph(dbg)
    
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
    parser = argparse.ArgumentParser(description="De Bruijn Graph Genome Assembly")
    
    parser.add_argument("--fastq", "-f", required=True, help="Input FASTQ file")
    parser.add_argument("--kmer", "-k", type=int, required=True, help="K-mer length")
    parser.add_argument("--output", "-o", required=True, help="Output FASTA file")
    parser.add_argument("--gfa", "-g", help="Output GFA file (optional)")
    parser.add_argument("--min-length", "-m", type=int, help="Minimum contig length")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    try:
        run_assembly(
            args.fastq,
            args.kmer,
            args.output,
            args.gfa,
            args.min_length,
            args.verbose
        )
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())