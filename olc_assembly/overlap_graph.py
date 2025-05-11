# olc_assembly/overlap_graph.py
"""
Construct and manipulate overlap graphs.
"""
from typing import Dict, Tuple, List, Set, Optional
import networkx as nx
from fastq_parser import Read

class OverlapGraph:
    """Class for constructing and manipulating overlap graphs."""
    
    def __init__(self):
        """Initialize the overlap graph."""
        self.graph = nx.DiGraph()
    
    def build_from_overlaps(self, reads: Dict[str, Read], 
                           overlaps: Dict[Tuple[str, str], Tuple[int, int]]):
        """
        Build the overlap graph from reads and their overlaps.
        
        Args:
            reads (dict): Dictionary mapping read IDs to Read objects.
            overlaps (dict): Dictionary mapping read ID pairs to (overlap_length, errors).
        """
        # Add nodes for each read
        for read_id, read in reads.items():
            self.graph.add_node(read_id, sequence=read.sequence, length=len(read))
        
        # Add edges for overlaps
        for (read_id1, read_id2), (overlap, errors) in overlaps.items():
            self.graph.add_edge(read_id1, read_id2, overlap=overlap, errors=errors, 
                               weight=-overlap)  # Negative weight for path finding
    
    def remove_transitive_edges(self):
        """
        Remove transitive edges from the graph.
        
        A transitive edge is one that can be replaced by a path through other nodes.
        """
        transitive_edges = []
        
        # For each node, examine all pairs of neighbors
        for node in self.graph.nodes():
            # Get all successors (outgoing edges)
            successors = list(self.graph.successors(node))
            
            # For each pair of successors
            for i in range(len(successors)):
                for j in range(len(successors)):
                    if i != j:
                        succ_i = successors[i]
                        succ_j = successors[j]
                        
                        # If there is a path from succ_i to succ_j, the edge from node to succ_j is transitive
                        if self.graph.has_edge(succ_i, succ_j):
                            # Check if the transitive path has similar or better overlap
                            direct_overlap = self.graph[node][succ_j]['overlap']
                            path_overlap = min(self.graph[node][succ_i]['overlap'], 
                                              self.graph[succ_i][succ_j]['overlap'])
                            
                            if path_overlap >= direct_overlap * 0.8:  # 80% threshold
                                transitive_edges.append((node, succ_j))
        
        # Remove the transitive edges
        for u, v in transitive_edges:
            if self.graph.has_edge(u, v):  # Check if edge still exists
                self.graph.remove_edge(u, v)
    
    def get_contigs_paths(self) -> List[List[str]]:
        """
        Find paths through the graph that correspond to contigs.
        
        Returns:
            list: List of paths, where each path is a list of read IDs.
        """
        # Find sources (nodes with no incoming edges)
        sources = [node for node in self.graph.nodes() if self.graph.in_degree(node) == 0]
        
        # Find sinks (nodes with no outgoing edges)
        sinks = [node for node in self.graph.nodes() if self.graph.out_degree(node) == 0]
        
        # Find all paths from sources to sinks
        all_paths = []
        for source in sources:
            for sink in sinks:
                try:
                    # Use simple_paths to find all simple paths
                    paths = list(nx.all_simple_paths(self.graph, source, sink))
                    all_paths.extend(paths)
                except nx.NetworkXNoPath:
                    continue
        
        # If no source-to-sink paths, try to find paths starting from any node
        if not all_paths:
            visited = set()
            for node in self.graph.nodes():
                if node not in visited:
                    # Find the longest path starting from this node
                    path = self._find_longest_path_from_node(node)
                    if path:
                        all_paths.append(path)
                        visited.update(path)
        
        return all_paths
    
    def _find_longest_path_from_node(self, start_node: str) -> List[str]:
        """
        Find the longest path starting from a given node.
        
        Args:
            start_node (str): Starting node.
            
        Returns:
            list: Longest path as a list of node IDs.
        """
        # Use a modified BFS to find the longest path
        queue = [(start_node, [start_node])]
        longest_path = [start_node]
        
        while queue:
            node, path = queue.pop(0)
            
            # If this path is longer than the current longest, update
            if len(path) > len(longest_path):
                longest_path = path
            
            # Add all successors to the queue
            for succ in self.graph.successors(node):
                if succ not in path:  # Avoid cycles
                    queue.append((succ, path + [succ]))
        
        return longest_path
    
    def find_nonbranching_paths(self) -> List[List[str]]:
        """
        Find maximal non-branching paths in the graph.
        
        Returns:
            list: List of non-branching paths.
        """
        paths = []
        visited_edges = set()
        
        # Process each node
        for node in self.graph.nodes():
            # If node has indegree and outdegree 1, it's part of a non-branching path
            if self.graph.in_degree(node) == 1 and self.graph.out_degree(node) == 1:
                continue  # Skip for now - we'll process these when we encounter a branching node
            
            # If node has outdegree > 0, start a new path
            if self.graph.out_degree(node) > 0:
                for successor in self.graph.successors(node):
                    if (node, successor) not in visited_edges:
                        # Start a new path
                        path = [node, successor]
                        visited_edges.add((node, successor))
                        
                        # Extend path until we hit a branching node
                        current = successor
                        while (self.graph.in_degree(current) == 1 and 
                               self.graph.out_degree(current) == 1):
                            next_node = list(self.graph.successors(current))[0]
                            if (current, next_node) in visited_edges:
                                break
                            
                            path.append(next_node)
                            visited_edges.add((current, next_node))
                            current = next_node
                        
                        paths.append(path)
        
        # Process isolated cycles
        for node in self.graph.nodes():
            if (self.graph.in_degree(node) == 1 and 
                self.graph.out_degree(node) == 1):
                
                predecessor = list(self.graph.predecessors(node))[0]
                successor = list(self.graph.successors(node))[0]
                
                if (predecessor, node) not in visited_edges and node == successor:
                    # This is a self-loop
                    paths.append([node, node])
                    visited_edges.add((node, node))
                
                elif ((predecessor, node) not in visited_edges and 
                      (node, successor) not in visited_edges):
                    # Start of an unvisited cycle
                    path = [node]
                    current = node
                    
                    while True:
                        next_node = list(self.graph.successors(current))[0]
                        if (current, next_node) in visited_edges:
                            break
                        
                        visited_edges.add((current, next_node))
                        if next_node == node:
                            break
                        
                        path.append(next_node)
                        current = next_node
                    
                    path.append(node)  # Close the cycle
                    paths.append(path)
        
        return paths