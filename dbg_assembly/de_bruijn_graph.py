# dbg_assembly/de_bruijn_graph.py
"""
Implementation of De Bruijn Graph for genome assembly.
"""
import networkx as nx
from collections import defaultdict, Counter

class DeBruijnGraph:
    """De Bruijn Graph implementation for genome assembly."""
    
    def __init__(self, k):
        """
        Initialize the De Bruijn Graph.
        
        Args:
            k (int): K-mer length.
        """
        self.k = k
        self.graph = nx.MultiDiGraph()
        self.kmer_counts = Counter()
    
    def add_kmer(self, kmer):
        """
        Add a k-mer to the graph.
        
        Args:
            kmer (str): K-mer to add.
        """
        self.kmer_counts[kmer] += 1
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        # Add nodes if they don't exist
        if not self.graph.has_node(prefix):
            self.graph.add_node(prefix)
        if not self.graph.has_node(suffix):
            self.graph.add_node(suffix)
        
        # Add the edge from prefix to suffix
        if self.graph.has_edge(prefix, suffix):
            # Increment the weight if edge already exists
            self.graph[prefix][suffix][0]['weight'] += 1
        else:
            # Create a new edge with weight 1
            self.graph.add_edge(prefix, suffix, weight=1, kmer=kmer)
    
    def add_sequence(self, sequence):
        """
        Add all k-mers from a sequence to the graph.
        
        Args:
            sequence (str): DNA sequence.
        """
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            self.add_kmer(kmer)
    
    def get_node_degrees(self):
        """
        Get in-degree and out-degree for each node.
        
        Returns:
            dict: Dictionary mapping node to (in_degree, out_degree) tuple.
        """
        degrees = {}
        for node in self.graph.nodes():
            in_degree = self.graph.in_degree(node)
            out_degree = self.graph.out_degree(node)
            degrees[node] = (in_degree, out_degree)
        return degrees
    
    def get_source_nodes(self):
        """
        Get nodes with in-degree 0 (potential starting points).
        
        Returns:
            list: List of source nodes.
        """
        return [n for n in self.graph.nodes() if self.graph.in_degree(n) == 0]
    
    def get_sink_nodes(self):
        """
        Get nodes with out-degree 0 (potential ending points).
        
        Returns:
            list: List of sink nodes.
        """
        return [n for n in self.graph.nodes() if self.graph.out_degree(n) == 0]
    
    def get_unbalanced_nodes(self):
        """
        Get nodes where in-degree != out-degree.
        
        Returns:
            dict: Dictionary mapping node to (in_degree, out_degree) tuple.
        """
        degrees = self.get_node_degrees()
        return {node: degree for node, degree in degrees.items() 
                if degree[0] != degree[1]}
    
    def is_eulerian(self):
        """
        Check if the graph is Eulerian (has an Eulerian path).
        
        Returns:
            bool: True if the graph is Eulerian, False otherwise.
        """
        # First, check if the graph is connected
        if not nx.is_strongly_connected(self.graph):
            return False
        
        # Check if all nodes have equal in and out degrees
        for node in self.graph.nodes():
            if self.graph.in_degree(node) != self.graph.out_degree(node):
                return False
        
        return True
    
    def is_semi_eulerian(self):
        """
        Check if the graph has an Eulerian path (but not necessarily a cycle).
        
        Returns:
            bool: True if the graph has an Eulerian path, False otherwise.
        """
        # Count unbalanced nodes
        unbalanced = self.get_unbalanced_nodes()
        
        # For a graph to have an Eulerian path, it must have at most 2 unbalanced nodes
        # and if there are 2, one must have in-degree = out-degree + 1 (end node)
        # and the other must have out-degree = in-degree + 1 (start node)
        if len(unbalanced) > 2:
            return False
        
        if len(unbalanced) == 2:
            # Check if one node has in_degree = out_degree + 1 and
            # the other has out_degree = in_degree + 1
            start_found = False
            end_found = False
            
            for node, (in_degree, out_degree) in unbalanced.items():
                if in_degree == out_degree + 1:
                    end_found = True
                elif out_degree == in_degree + 1:
                    start_found = True
            
            return start_found and end_found
        
        # If there are 0 or 1 unbalanced nodes, the graph has an Eulerian path
        return True