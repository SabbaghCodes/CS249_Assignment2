# olc_assembly/consensus_builder.py
"""
Build consensus sequences from layouts.
"""
from typing import Dict, List, Tuple
from collections import Counter
from fastq_parser import Read

class ConsensusBuilder:
    """Class for building consensus sequences from layouts."""
    
    def __init__(self, reads: Dict[str, Read]):
        """
        Initialize the consensus builder.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
        """
        self.reads = reads
    
    def build_consensus(self, layout: List[Tuple[str, int]]) -> str:
        """
        Build a consensus sequence from a layout.
        
        Args:
            layout (list): List of (read_id, position) tuples.
            
        Returns:
            str: Consensus sequence.
        """
        # Find the total length of the contig
        max_end = 0
        for read_id, pos in layout:
            read_end = pos + len(self.reads[read_id])
            max_end = max(max_end, read_end)
        
        # Initialize an array to hold base counts at each position
        base_counts = [{} for _ in range(max_end)]
        
        # Add bases from each read to the counts
        for read_id, pos in layout:
            read_seq = self.reads[read_id].sequence
            for i, base in enumerate(read_seq):
                position = pos + i
                if position < max_end:
                    counts = base_counts[position]
                    counts[base] = counts.get(base, 0) + 1
        
        # Build the consensus sequence
        consensus = []
        for counts in base_counts:
            if not counts:
                # No coverage at this position
                consensus.append('N')
            else:
                # Get the most common base
                most_common_base = max(counts.items(), key=lambda x: x[1])[0]
                consensus.append(most_common_base)
        
        return ''.join(consensus)
    
    def build_consensus_with_quality(self, layout: List[Tuple[str, int]]) -> str:
        """
        Build a consensus sequence using quality scores.
        
        Args:
            layout (list): List of (read_id, position) tuples.
            
        Returns:
            str: Consensus sequence.
        """
        # Find the total length of the contig
        max_end = 0
        for read_id, pos in layout:
            read_end = pos + len(self.reads[read_id])
            max_end = max(max_end, read_end)
        
        # Initialize an array to hold weighted base counts at each position
        weighted_counts = [{} for _ in range(max_end)]
        
        # Add weighted bases from each read to the counts
        for read_id, pos in layout:
            read = self.reads[read_id]
            for i, base in enumerate(read.sequence):
                position = pos + i
                if position < max_end:
                    # Convert quality score to error probability
                    qual_char = read.quality[i]
                    quality = ord(qual_char) - 33  # PHRED score
                    error_prob = 10 ** (-quality / 10)
                    weight = 1 - error_prob
                    
                    counts = weighted_counts[position]
                    counts[base] = counts.get(base, 0) + weight
        
        # Build the consensus sequence
        consensus = []
        for counts in weighted_counts:
            if not counts:
                # No coverage at this position
                consensus.append('N')
            else:
                # Get the base with the highest weighted count
                best_base = max(counts.items(), key=lambda x: x[1])[0]
                consensus.append(best_base)
        
        return ''.join(consensus)