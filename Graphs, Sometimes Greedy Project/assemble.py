import graph
import random

log = False

def overlap(s1, s2):
    """The length of the longest match between the a suffix of s1 and a prefix of s2"""
    for i in range(min(len(s1), len(s2)), 0, -1):
        if s1.endswith(s2[:i]):
            return i
    return 0

class OverlapGraph(graph.ComponentTrackingGraph,
                   graph.VertexLabeledDirectedGraph,
                   graph.EdgeWeightedDirectedGraph,                   
                   graph.AdjacencyListDirectedGraph):
    """A graph for representing reads and read overlaps for fragment assembly."""
    
    def superstring(self):
        """Returns the superstring represented this graph, assuming it is connected"""
        assert(self.num_components() == 1)
        substrings = []
        last_overlap = 0
        i = self.source()
        while True:
            substrings.append(self.vertex_label(i)[last_overlap:])
            if self.outdegree(i) > 0:
                j = self.out_edges(i)[0][1]
                last_overlap = self.edge_weight(i, j)
                i = j
            else:
                break
        return "".join(substrings)
    
    def source(self):
        """Returns the first vertex in the graph that has indegree == 0"""
        for i in range(self.num_vertices()):
            if self.indegree(i) == 0:
                return i
        else:
            return None

def sometimes_greedy_assemble(reads, p=1.0):
    # Construct initial read graph without any edges
    g = OverlapGraph(len(reads))
    for i, read in enumerate(reads):
        g.set_vertex_label(i, read)
        
    # Construct all possible edges
    print("Computing overlaps...")
    q = []
    for i in range(len(reads)):
        for j in range(len(reads)):
            if i != j:
                overlap_length = overlap(reads[i], reads[j])
                q.append((i, j, overlap_length))
    
    # Construct a queue of edges
    # Sort edges by decreasing weight, with lexicographical comparison of
    # read strings of edge vertices as tiebreaker              
    q.sort(key=lambda edge: (-edge[2], reads[edge[0]], reads[edge[1]]),
           reverse=True)
    
    print("Running greedy algorithm...")
    h = []
    while g.num_components() > 1:
        #if log: print(q)
        e = q.pop()
        i, j, weight = e
        if (g.outdegree(i) == 0 and
            g.indegree(j) == 0 and
            not g.are_in_same_component(i, j)):
            x = random.random()
            if x < p:
                if log: print("adding edge")
                g.add_weighted_edge(i, j, weight)
                while h:
                    q.append(h.pop())
            else:
                if log: print("skipping edge")
                h.append(e)
        else:
            if log: print("invalid edge")
        if not q:
            if log: print("q is empty")
            while h:
                q.append(h.pop())
                
    # return the superstring represented by the graph
    return g.superstring()
