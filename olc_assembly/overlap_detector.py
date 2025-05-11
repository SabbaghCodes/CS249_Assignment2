# olc_assembly/overlap_detector.py
"""
Functions for detecting overlaps between reads.
"""
from typing import Tuple, List, Dict, Set, Optional
import numpy as np
from fastq_parser import Read

class OverlapDetector:
    """Class for detecting overlaps between reads."""
    
    def __init__(self, min_overlap_length: int, max_error_rate: float = 0.05):
        """
        Initialize the overlap detector.
        
        Args:
            min_overlap_length (int): Minimum length of overlap to consider.
            max_error_rate (float): Maximum allowed error rate in overlaps.
        """
        self.min_overlap_length = min_overlap_length
        self.max_error_rate = max_error_rate
    
    def find_overlap(self, read1: Read, read2: Read) -> Tuple[int, int]:
        """
        Find the best suffix-prefix overlap between two reads.
        
        Args:
            read1 (Read): First read (suffix).
            read2 (Read): Second read (prefix).
            
        Returns:
            tuple: (overlap_length, errors) or (0, 0) if no significant overlap found.
        """
        seq1 = read1.sequence
        seq2 = read2.sequence
        
        max_overlap = min(len(seq1), len(seq2))
        best_overlap = 0
        best_errors = 0
        
        # Try different overlap lengths
        for overlap_len in range(min(max_overlap, len(seq1)), self.min_overlap_length - 1, -1):
            # Get the suffix of seq1 and prefix of seq2
            suffix = seq1[-overlap_len:]
            prefix = seq2[:overlap_len]
            
            # Count mismatches
            errors = sum(s != p for s, p in zip(suffix, prefix))
            error_rate = errors / overlap_len
            
            # Check if this overlap is valid
            if error_rate <= self.max_error_rate:
                best_overlap = overlap_len
                best_errors = errors
                break
        
        return best_overlap, best_errors
    
    def compute_all_overlaps(self, reads: Dict[str, Read]) -> Dict[Tuple[str, str], Tuple[int, int]]:
        """
        Compute all pairwise overlaps between reads.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            
        Returns:
            dict: Dictionary mapping read ID pairs to (overlap_length, errors).
        """
        overlaps = {}
        read_ids = list(reads.keys())
        
        for i in range(len(read_ids)):
            for j in range(len(read_ids)):
                if i != j:  # Don't compute overlap of a read with itself
                    read1 = reads[read_ids[i]]
                    read2 = reads[read_ids[j]]
                    
                    overlap, errors = self.find_overlap(read1, read2)
                    
                    if overlap >= self.min_overlap_length:
                        overlaps[(read_ids[i], read_ids[j])] = (overlap, errors)
        
        return overlaps
    
    def filter_contained_reads(self, reads: Dict[str, Read], 
                               overlaps: Dict[Tuple[str, str], Tuple[int, int]]) -> Set[str]:
        """
        Identify reads that are contained within other reads.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            overlaps (dict): Dictionary mapping read ID pairs to (overlap_length, errors).
            
        Returns:
            set: Set of read IDs that are contained within other reads.
        """
        contained_reads = set()
        
        for (read_id1, read_id2), (overlap, _) in overlaps.items():
            read1 = reads[read_id1]
            read2 = reads[read_id2]
            
            # If the overlap is the entire length of read1, then read1 is contained in read2
            if overlap == len(read1) and len(read1) < len(read2):
                contained_reads.add(read_id1)
            
            # If the overlap is the entire length of read2, then read2 is contained in read1
            if overlap == len(read2) and len(read2) < len(read1):
                contained_reads.add(read_id2)
        
        return contained_reads
    
    def get_dovetail_overlaps(self, reads: Dict[str, Read], 
                             overlaps: Dict[Tuple[str, str], Tuple[int, int]],
                             contained_reads: Set[str]) -> Dict[Tuple[str, str], Tuple[int, int]]:
        """
        Filter overlaps to keep only dovetail overlaps (suffix-prefix).
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            overlaps (dict): Dictionary mapping read ID pairs to (overlap_length, errors).
            contained_reads (set): Set of read IDs that are contained within other reads.
            
        Returns:
            dict: Dictionary mapping read ID pairs to (overlap_length, errors) for dovetail overlaps.
        """
        dovetail_overlaps = {}
        
        for (read_id1, read_id2), (overlap, errors) in overlaps.items():
            # Skip contained reads
            if read_id1 in contained_reads or read_id2 in contained_reads:
                continue
            
            read1 = reads[read_id1]
            read2 = reads[read_id2]
            
            # A dovetail overlap is one where:
            # 1. The overlap is less than the length of both reads
            # 2. The overlap is a suffix of read1 and a prefix of read2
            if overlap < len(read1) and overlap < len(read2):
                dovetail_overlaps[(read_id1, read_id2)] = (overlap, errors)
        
        return dovetail_overlaps