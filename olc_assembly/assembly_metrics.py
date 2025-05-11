# olc_assembly/assembly_metrics.py
"""
Calculate metrics for genome assembly.
"""
import numpy as np

class AssemblyMetrics:
    """Calculate metrics for genome assembly."""
    
    @staticmethod
    def calculate_basic_metrics(contigs):
        """
        Calculate basic assembly metrics.
        
        Args:
            contigs (list): List of contig sequences.
            
        Returns:
            dict: Dictionary of basic metrics.
        """
        # Calculate lengths
        lengths = [len(contig) for contig in contigs]
        total_length = sum(lengths)
        
        if not lengths:
            return {
                'total_length': 0,
                'num_contigs': 0,
                'mean_length': 0,
                'median_length': 0,
                'max_length': 0,
                'min_length': 0,
                'n50': 0,
                'n90': 0,
                'l50': 0,
                'gc_content': 0
            }
        
        # Sort lengths in descending order
        lengths.sort(reverse=True)
        
        # Calculate N50 and N90
        cumulative_length = 0
        n50 = None
        n90 = None
        l50 = None
        
        for i, length in enumerate(lengths):
            cumulative_length += length
            
            if n50 is None and cumulative_length >= total_length * 0.5:
                n50 = length
                l50 = i + 1
                
            if n90 is None and cumulative_length >= total_length * 0.9:
                n90 = length
                break
        
        # Calculate GC content
        gc_count = 0
        for contig in contigs:
            gc_count += contig.count('G') + contig.count('g') + contig.count('C') + contig.count('c')
        
        gc_content = gc_count / total_length if total_length > 0 else 0
        
        return {
            'total_length': total_length,
            'num_contigs': len(contigs),
            'mean_length': total_length / len(contigs) if contigs else 0,
            'median_length': np.median(lengths),
            'max_length': max(lengths) if lengths else 0,
            'min_length': min(lengths) if lengths else 0,
            'n50': n50 or 0,
            'n90': n90 or 0,
            'l50': l50 or 0,
            'gc_content': gc_content
        }
    
    @staticmethod
    def compare_to_reference(contigs, reference):
        """
        Compare assembly to a reference sequence.
        
        Args:
            contigs (list): List of contig sequences.
            reference (str): Reference sequence.
            
        Returns:
            dict: Dictionary of comparison metrics.
        """
        # This would typically involve alignment, but let's implement a simpler version
        # In a real implementation, you'd use a tool like BLAST or minimap2
        
        # For now, let's just check how many k-mers from the reference are in our contigs
        k = 31  # Typical k-mer size for comparison
        
        # Generate k-mers from reference
        ref_kmers = set()
        for i in range(len(reference) - k + 1):
            ref_kmers.add(reference[i:i+k])
        
        # Check how many reference k-mers are in contigs
        found_kmers = set()
        for contig in contigs:
            for i in range(len(contig) - k + 1):
                kmer = contig[i:i+k]
                if kmer in ref_kmers:
                    found_kmers.add(kmer)
        
        # Calculate metrics
        kmer_recall = len(found_kmers) / len(ref_kmers) if ref_kmers else 0
        
        return {
            'reference_length': len(reference),
            'kmer_recall': kmer_recall,
            'estimated_genome_fraction': kmer_recall  # Very simplistic estimate
        }