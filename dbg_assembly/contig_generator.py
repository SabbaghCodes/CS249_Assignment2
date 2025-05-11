# dbg_assembly/contig_generator.py
"""
Generate contigs from paths in the De Bruijn Graph.
"""

class ContigGenerator:
    """Generate contigs from paths in the De Bruijn Graph."""
    
    @staticmethod
    def path_to_contig(path, k):
        """
        Convert a path of nodes to a contig sequence.
        
        Args:
            path (list): List of nodes in the path.
            k (int): K-mer length.
            
        Returns:
            str: Contig sequence.
        """
        if not path:
            return ""
        
        # First node contributes all its characters
        contig = path[0]
        
        # Subsequent nodes contribute only their last character
        for node in path[1:]:
            contig += node[-1]
        
        return contig
    
    @staticmethod
    def generate_contigs_from_paths(paths, k):
        """
        Generate contigs from a list of paths.
        
        Args:
            paths (list): List of paths (each path is a list of nodes).
            k (int): K-mer length.
            
        Returns:
            list: List of contigs (strings).
        """
        return [ContigGenerator.path_to_contig(path, k) for path in paths]
    
    @staticmethod
    def filter_contigs(contigs, min_length=None):
        """
        Filter contigs based on minimum length.
        
        Args:
            contigs (list): List of contigs.
            min_length (int, optional): Minimum contig length.
            
        Returns:
            list: Filtered list of contigs.
        """
        if min_length is None:
            return contigs
            
        return [contig for contig in contigs if len(contig) >= min_length]