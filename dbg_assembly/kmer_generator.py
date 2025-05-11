# dbg_assembly/kmer_generator.py
"""
Generate k-mers from DNA sequences.
"""

class KmerGenerator:
    """Generate k-mers from DNA sequences."""
    
    @staticmethod
    def generate_kmers(sequence, k):
        """
        Generate k-mers from a DNA sequence.
        
        Args:
            sequence (str): Input DNA sequence.
            k (int): K-mer length.
            
        Returns:
            list: List of k-mers.
        """
        return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    
    @staticmethod
    def get_prefix(kmer, k):
        """
        Get the prefix of length k-1 from a k-mer.
        
        Args:
            kmer (str): K-mer.
            k (int): K-mer length.
            
        Returns:
            str: Prefix of length k-1.
        """
        return kmer[:-1]
    
    @staticmethod
    def get_suffix(kmer, k):
        """
        Get the suffix of length k-1 from a k-mer.
        
        Args:
            kmer (str): K-mer.
            k (int): K-mer length.
            
        Returns:
            str: Suffix of length k-1.
        """
        return kmer[1:]