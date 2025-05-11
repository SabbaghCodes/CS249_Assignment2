"""
FASTQ file parser for reading sequence data.
"""

class FastqParser:
    """Parser for FASTQ formatted files."""
    
    def __init__(self, file_path):
        """
        Initialize the FASTQ parser.
        
        Args:
            file_path (str): Path to the FASTQ file.
        """
        self.file_path = file_path
    
    def parse(self):
        """
        Parse the FASTQ file and yield sequences.
        
        Yields:
            str: DNA sequence from the FASTQ file.
        """
        with open(self.file_path, 'r') as file:
            line_num = 0
            for line in file:
                line_num += 1
                # Sequence lines are the second line of each 4-line group
                if line_num % 4 == 2:
                    yield line.strip()
    
    def parse_with_quality(self):
        """
        Parse the FASTQ file and yield sequences with their quality scores.
        
        Yields:
            tuple: (sequence_id, sequence, quality_scores)
        """
        with open(self.file_path, 'r') as file:
            while True:
                header = file.readline().strip()
                if not header:
                    break
                    
                sequence = file.readline().strip()
                file.readline()  # Skip the '+' line
                quality = file.readline().strip()
                
                if header.startswith('@'):
                    seq_id = header[1:]  # Remove the '@' prefix
                else:
                    seq_id = header
                    
                yield (seq_id, sequence, quality)