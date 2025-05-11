# dbg_assembly/gfa_export.py
"""
Export De Bruijn Graph to GFA format for visualization.
"""

class GfaExporter:
    """Export De Bruijn Graph to GFA format."""
    
    def __init__(self, file_path):
        """
        Initialize the GFA exporter.
        
        Args:
            file_path (str): Path to the output GFA file.
        """
        self.file_path = file_path
    
    def export_graph(self, dbg):
        """
        Export the De Bruijn Graph to GFA format.
        
        Args:
            dbg (DeBruijnGraph): De Bruijn Graph to export.
        """
        with open(self.file_path, 'w') as f:
            # Write header
            f.write("H\tVN:Z:1.0\n")
            
            # Write segments (nodes)
            for node in dbg.graph.nodes():
                # In GFA, each segment needs a unique identifier
                f.write(f"S\t{node}\t{node}\n")
            
            # Write links (edges)
            for u, v, data in dbg.graph.edges(data=True):
                # In the DBG, edges represent k-mers
                kmer = data.get('kmer', '')
                weight = data.get('weight', 1)
                
                # GFA link format: L <segment1> <+/-> <segment2> <+/-> <CIGAR>
                f.write(f"L\t{u}\t+\t{v}\t+\t{len(u)-1}M\tRC:i:{weight}\n")
    
    def export_contigs(self, contigs, k):
        """
        Export contigs to GFA format.
        
        Args:
            contigs (list): List of contig sequences.
            k (int): K-mer length used in assembly.
        """
        with open(self.file_path, 'w') as f:
            # Write header
            f.write("H\tVN:Z:1.0\n")
            
            # Write segments (contigs)
            for i, contig in enumerate(contigs):
                segment_id = f"contig_{i+1}"
                f.write(f"S\t{segment_id}\t{contig}\n")