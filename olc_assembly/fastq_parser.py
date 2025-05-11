# olc_assembly/fastq_parser.py
"""
FASTQ file parser for reading sequence data.
"""
from dataclasses import dataclass

@dataclass
class Read:
    """Class for storing a read with its ID, sequence, and quality scores."""
    id: str
    sequence: str
    quality: str
    
    def __len__(self):
        """Return the length of the read sequence."""
        return len(self.sequence)

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
        Parse the FASTQ file and yield Read objects.
        
        Yields:
            Read: Read object containing ID, sequence, and quality scores.
        """
        with open(self.file_path, 'r') as file:
            while True:
                # Read the four lines of a FASTQ entry
                header = file.readline().strip()
                if not header:
                    break  # End of file
                
                sequence = file.readline().strip()
                file.readline()  # Skip the '+' line
                quality = file.readline().strip()
                
                if header.startswith('@'):
                    read_id = header[1:]  # Remove the '@' prefix
                else:
                    read_id = header
                
                yield Read(id=read_id, sequence=sequence, quality=quality)
    
    def parse_as_dict(self):
        """
        Parse the FASTQ file and return reads as a dictionary.
        
        Returns:
            dict: Dictionary mapping read IDs to Read objects.
        """
        reads = {}
        for read in self.parse():
            reads[read.id] = read
        
        return reads