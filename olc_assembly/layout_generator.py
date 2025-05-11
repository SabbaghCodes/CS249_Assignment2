# olc_assembly/layout_generator.py
"""
Generate layouts from paths in the overlap graph.
"""
from typing import Dict, List, Tuple
from fastq_parser import Read

class LayoutGenerator:
    """Class for generating layouts from paths in the overlap graph."""
    
    def __init__(self, reads: Dict[str, Read], overlaps: Dict[Tuple[str, str], Tuple[int, int]]):
        """
        Initialize the layout generator.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            overlaps (dict): Dictionary mapping read ID pairs to (overlap_length, errors).
        """
        self.reads = reads
        self.overlaps = overlaps
    
    def generate_layout(self, path: List[str]) -> List[Tuple[str, int]]:
        """
        Generate a layout from a path of read IDs.
        
        A layout specifies the position of each read relative to the start of the contig.
        
        Args:
            path (list): List of read IDs representing a path through the graph.
            
        Returns:
            list: List of (read_id, position) tuples.
        """
        layout = []
        current_pos = 0
        
        # Add the first read at position 0
        layout.append((path[0], current_pos))
        
        # Process the rest of the reads in the path
        for i in range(1, len(path)):
            prev_read_id = path[i-1]
            curr_read_id = path[i]
            
            # Get the overlap between the previous and current read
            overlap_info = self.overlaps.get((prev_read_id, curr_read_id))
            if overlap_info:
                overlap_len, _ = overlap_info
            else:
                # If no direct overlap, use a default small overlap
                overlap_len = 1
            
            # Calculate the position of the current read
            current_pos += len(self.reads[prev_read_id]) - overlap_len
            layout.append((curr_read_id, current_pos))
        
        return layout