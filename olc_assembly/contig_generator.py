# olc_assembly/contig_generator.py
"""
Generate contigs from consensus sequences.
"""
from typing import Dict, List, Tuple
from fastq_parser import Read
from layout_generator import LayoutGenerator
from consensus_builder import ConsensusBuilder

class ContigGenerator:
    """Class for generating contigs from paths in the graph."""
    
    def __init__(self, reads: Dict[str, Read], overlaps: Dict[Tuple[str, str], Tuple[int, int]]):
        """
        Initialize the contig generator.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            overlaps (dict): Dictionary mapping read ID pairs to (overlap_length, errors).
        """
        self.reads = reads
        self.overlaps = overlaps
        self.layout_generator = LayoutGenerator(reads, overlaps)
        self.consensus_builder = ConsensusBuilder(reads)
    
    def generate_contigs(self, paths: List[List[str]]) -> List[str]:
        """
        Generate contigs from paths.
        
        Args:
            paths (list): List of paths, where each path is a list of read IDs.
            
        Returns:
            list: List of contig sequences.
        """
        contigs = []
        
        for path in paths:
            # Skip paths that are too short
            if len(path) < 2:
                continue
                
            # Generate layout
            layout = self.layout_generator.generate_layout(path)
            
            # Build consensus
            consensus = self.consensus_builder.build_consensus(layout)
            
            contigs.append(consensus)
        
        return contigs
    
    def generate_contigs_with_quality(self, paths: List[List[str]]) -> List[str]:
        """
        Generate contigs from paths using quality scores.
        
        Args:
            paths (list): List of paths, where each path is a list of read IDs.
            
        Returns:
            list: List of contig sequences.
        """
        contigs = []
        
        for path in paths:
            # Skip paths that are too short
            if len(path) < 2:
                continue
                
            # Generate layout
            layout = self.layout_generator.generate_layout(path)
            
            # Build consensus with quality
            consensus = self.consensus_builder.build_consensus_with_quality(layout)
            
            contigs.append(consensus)
        
        return contigs
    
    def filter_contigs(self, contigs: List[str], min_length: int = 0) -> List[str]:
        """
        Filter contigs based on minimum length.
        
        Args:
            contigs (list): List of contig sequences.
            min_length (int): Minimum contig length.
            
        Returns:
            list: Filtered list of contigs.
        """
        return [contig for contig in contigs if len(contig) >= min_length]