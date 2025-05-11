# olc_assembly/visualize.py
"""
Visualization functions for Overlap-Layout-Consensus assembly.
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class Visualizer:
    """Visualization for OLC assembly."""
    
    @staticmethod
    def plot_overlap_graph(graph, figsize=(12, 10), node_size=500, width=1.5):
        """
        Plot the overlap graph.
        
        Args:
            graph (nx.DiGraph): Overlap graph.
            figsize (tuple): Figure size.
            node_size (int): Size of nodes in the plot.
            width (float): Width of edges in the plot.
        """
        plt.figure(figsize=figsize)
        
        # Create a position layout for the graph
        pos = nx.spring_layout(graph, seed=42)
        
        # Draw the graph
        nx.draw(graph, pos, with_labels=True, node_color='lightblue',
                edge_color='gray', node_size=node_size, width=width,
                font_size=10, arrowsize=15)
        
        # Draw edge labels for overlaps
        edge_labels = {(u, v): f"{d['overlap']}" for u, v, d in graph.edges(data=True)}
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_size=8)
        
        plt.title("Overlap Graph")
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
    def plot_layout(layout, reads, figsize=(15, 8)):
        """
        Visualize the layout of reads in a contig.
        
        Args:
            layout (list): List of (read_id, position) tuples.
            reads (dict): Dictionary mapping read IDs to Read objects.
            figsize (tuple): Figure size.
        """
        plt.figure(figsize=figsize)
        
        # Sort layout by position
        sorted_layout = sorted(layout, key=lambda x: x[1])
        
        # Find max position
        max_pos = 0
        for read_id, pos in sorted_layout:
            max_pos = max(max_pos, pos + len(reads[read_id]))
        
        # Plot each read as a horizontal bar
        for i, (read_id, pos) in enumerate(sorted_layout):
            read_len = len(reads[read_id])
            plt.barh(i, read_len, left=pos, height=0.8, color='skyblue', alpha=0.7, 
                    edgecolor='black')
            plt.text(pos + read_len/2, i, read_id, ha='center', va='center', fontsize=8)
        
        plt.yticks(range(len(sorted_layout)), [f"Read {i+1}" for i in range(len(sorted_layout))])
        plt.xlabel('Position (bp)')
        plt.title('Read Layout in Contig')
        plt.xlim(0, max_pos + 10)
        plt.grid(axis='x', alpha=0.3)
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