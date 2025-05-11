# dbg_assembly/path_finder.py
"""
Find Eulerian paths in a De Bruijn Graph for genome assembly.
"""
import networkx as nx
import copy

class PathFinder:
    """Find paths in a De Bruijn Graph for assembly."""
    
    @staticmethod
    def find_eulerian_path(graph):
        """
        Find an Eulerian path in the graph.
        
        Args:
            graph (nx.MultiDiGraph): De Bruijn graph.
            
        Returns:
            list: Eulerian path as a list of nodes.
        """
        # Make a copy of the graph to avoid modifying the original
        g = graph.copy()
        
        # If graph is empty, return empty path
        if len(g.nodes()) == 0:
            return []
        
        # Find a start node (in_degree < out_degree if exists, otherwise any node)
        start_node = None
        for node in g.nodes():
            if g.in_degree(node) < g.out_degree(node):
                start_node = node
                break
        
        # If no such node exists, pick any node with outgoing edges
        if start_node is None:
            for node in g.nodes():
                if g.out_degree(node) > 0:
                    start_node = node
                    break
        
        # If we still don't have a start node, the graph has no edges
        if start_node is None:
            return list(g.nodes())
        
        # Find the Eulerian path using Hierholzer's algorithm
        path = []
        stack = [start_node]
        
        while stack:
            current = stack[-1]
            
            # If there are outgoing edges, follow one
            if g.out_degree(current) > 0:
                # Get the next edge and remove it
                neighbor = list(g.successors(current))[0]
                edge_data = g.get_edge_data(current, neighbor)[0]
                g.remove_edge(current, neighbor)
                
                # Add the neighbor to the stack
                stack.append(neighbor)
            else:
                # No more outgoing edges, add to path and backtrack
                path.append(stack.pop())
        
        # Reverse the path to get the correct order
        return path[::-1]
    
    @staticmethod
    def find_contigs(dbg):
        """
        Find contigs in the De Bruijn Graph by identifying non-branching paths.
        
        Args:
            dbg (DeBruijnGraph): De Bruijn graph.
            
        Returns:
            list: List of contigs (node paths).
        """
        graph = dbg.graph
        contigs = []
        
        # Find nodes that are branching points (in or out degree > 1)
        branching_nodes = set()
        for node in graph.nodes():
            if graph.in_degree(node) > 1 or graph.out_degree(node) > 1:
                branching_nodes.add(node)
        
        # Add source and sink nodes to branching points
        for node in graph.nodes():
            if graph.in_degree(node) == 0 or graph.out_degree(node) == 0:
                branching_nodes.add(node)
        
        # Find non-branching paths (maximal paths without branching)
        visited_edges = set()
        
        for node in branching_nodes:
            if graph.out_degree(node) > 0:
                for successor in graph.successors(node):
                    # Check if edge has been visited
                    edge_key = (node, successor, 0)  # Using key 0 for MultiDiGraph
                    if edge_key in visited_edges:
                        continue
                    
                    # Start a new path
                    path = [node, successor]
                    visited_edges.add(edge_key)
                    
                    # Extend path until we hit a branching node
                    current = successor
                    while (graph.in_degree(current) == 1 and 
                           graph.out_degree(current) == 1 and
                           current not in branching_nodes):
                        next_node = list(graph.successors(current))[0]
                        path.append(next_node)
                        visited_edges.add((current, next_node, 0))
                        current = next_node
                    
                    contigs.append(path)
        
        # Handle loops (cycles without branching)
        for node in graph.nodes():
            if (node not in branching_nodes and 
                graph.in_degree(node) == 1 and 
                graph.out_degree(node) == 1):
                
                # Check if node is part of an unvisited cycle
                successor = list(graph.successors(node))[0]
                edge_key = (node, successor, 0)
                if edge_key not in visited_edges:
                    # Start a new path
                    path = [node]
                    current = node
                    
                    # Follow path until we return to the start
                    while True:
                        next_node = list(graph.successors(current))[0]
                        visited_edges.add((current, next_node, 0))
                        
                        if next_node == node:
                            break
                            
                        path.append(next_node)
                        current = next_node
                    
                    path.append(node)  # Close the cycle
                    contigs.append(path)
        
        return contigs