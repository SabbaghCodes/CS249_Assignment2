# dbg_assembly/visualize.py
"""
Visualization functions for De Bruijn Graph assembly.
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class Visualizer:
    """Visualization for De Bruijn Graph assembly."""
    
    @staticmethod
    def plot_graph(dbg, figsize=(12, 10), node_size=500, width=1.5):
        """
        Plot the De Bruijn Graph.
        
        Args:
            dbg (DeBruijnGraph): De Bruijn Graph to plot.
            figsize (tuple): Figure size.
            node_size (int): Size of nodes in the plot.
            width (float): Width of edges in the plot.
        """
        plt.figure(figsize=figsize)
        
        # Create a position layout for the graph
        pos = nx.spring_layout(dbg.graph, seed=42)
        
        # Draw the graph
        nx.draw(dbg.graph, pos, with_labels=True, node_color='lightblue',
                edge_color='gray', node_size=node_size, width=width,
                font_size=10, arrowsize=15)
        
        plt.title(f"De Bruijn Graph (k={dbg.k})")
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_contig_length_distribution(contigs, bins=30, figsize=(10, 6)):
        """
        Plot the contig length distribution.
        
        Args:
            contigs (list): List of contig sequences.
            bins (int): Number of bins in the histogram.
            figsize (tuple): Figure size.
        """
        lengths = [len(contig) for contig in contigs]
        
        plt.figure(figsize=figsize)
        plt.hist(lengths, bins=bins, color='skyblue', edgecolor='black')
        plt.title('Contig Length Distribution')
        plt.xlabel('Contig Length (bp)')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_kmer_spectrum(dbg, figsize=(10, 6), max_count=100):
        """
        Plot the k-mer spectrum.
        
        Args:
            dbg (DeBruijnGraph): De Bruijn Graph.
            figsize (tuple): Figure size.
            max_count (int): Maximum k-mer count to include in plot.
        """
        # Get k-mer counts
        counts = list(dbg.kmer_counts.values())
        counts = [c for c in counts if c <= max_count]  # Filter extreme values
        
        plt.figure(figsize=figsize)
        plt.hist(counts, bins=min(50, max(10, len(set(counts)))), 
                 color='lightgreen', edgecolor='black')
        plt.title(f'K-mer Spectrum (k={dbg.k})')
        plt.xlabel('K-mer Count')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def plot_metrics_comparison(metrics_list, labels, figsize=(12, 8)):
        """
        Plot a comparison of assembly metrics.
        
        Args:
            metrics_list (list): List of metrics dictionaries.
            labels (list): Labels for each set of metrics.
            figsize (tuple): Figure size.
        """
        metrics_to_plot = ['num_contigs', 'total_length', 'n50', 'l50', 'gc_content']
        
        fig, axes = plt.subplots(len(metrics_to_plot), 1, figsize=figsize)
        
        for i, metric in enumerate(metrics_to_plot):
            values = [metrics.get(metric, 0) for metrics in metrics_list]
            axes[i].bar(labels, values, color='skyblue')
            axes[i].set_title(f'{metric.replace("_", " ").title()}')
            axes[i].set_ylabel('Value')
            axes[i].grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.show()