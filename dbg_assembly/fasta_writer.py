# dbg_assembly/fasta_writer.py
"""
Write contigs to a FASTA file.
"""

class FastaWriter:
    """Write contigs to a FASTA file."""
    
    def __init__(self, file_path):
        """
        Initialize the FASTA writer.
        
        Args:
            file_path (str): Path to the output FASTA file.
        """
        self.file_path = file_path
    
    def write_contigs(self, contigs):
        """
        Write contigs to the FASTA file.
        
        Args:
            contigs (list): List of contig sequences.
        """
        with open(self.file_path, 'w') as f:
            for i, contig in enumerate(contigs):
                header = f">contig_{i+1} length={len(contig)}"
                f.write(f"{header}\n")
                
                # Write sequence in lines of 60 characters
                for j in range(0, len(contig), 60):
                    f.write(f"{contig[j:j+60]}\n")
    
    def write_contigs_with_headers(self, contigs_with_headers):
        """
        Write contigs with custom headers to the FASTA file.
        
        Args:
            contigs_with_headers (list): List of (header, sequence) tuples.
        """
        with open(self.file_path, 'w') as f:
            for header, sequence in contigs_with_headers:
                f.write(f">{header}\n")
                
                # Write sequence in lines of 60 characters
                for j in range(0, len(sequence), 60):
                    f.write(f"{sequence[j:j+60]}\n")