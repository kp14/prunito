"""Subclass of DiGraph to represent taxonomic trees."""
from networkx.classes.digraph import DiGraph
from networkx.readwrite import json_graph
import json
from collections import Counter

class TaxGraph(DiGraph):
    """Subclass of DiGraph to represent taxonomic trees.

    To be used with taxonomic node data taken from UniProt entries.
    The Biopython parser returns the OC line data from UniProt entries
    as a list whose order is the same as in the entry, i.e. it reflects
    the hierarchy.

    Defines a new method, add_edges_from_ordered_nodes(), to pair up
    adjacent nodes as edges which can then be used for construction of
    a DiGraph.

    Graphs can be exported to JSON.

    """

    def __init__(self):
        """Initialize as an empty DiGraph by calling the DiGraph initializer."""
        super(TaxGraph, self).__init__()
        self.add_node("root")

    def add_edges_from_ordered_nodes(self, node_list):
        """Pairs up successive nodes in a list of nodes to form edge pairs.
        This assumes that the order of the list reflects a hierarchy.

        Parameter:
        node_list: list containing an hierachically ordered number of nodes
                    (mandatory)
        """
        for node in node_list:
            idx = node_list.index(node)
            if idx == 0:
                self.add_edge("root", node_list[idx], weight=1)
            else:
                self.add_edge(node_list[idx-1], node_list[idx], weight=1)

    #overrides add_edge() from DiGraph
    def add_edge(self, u, v, attr_dict=None, **attr):
        """Add an edge between u and v.

        attr_dict is implemented as a Counter.

        The nodes u and v will be automatically added if they are
        not already in the graph.

        Edge attributes can be specified with keywords or by providing
        a dictionary with key/value pairs.  See examples below.

        Parameters
        ----------
        u,v : nodes
            Nodes can be, for example, strings or numbers.
            Nodes must be hashable (and not None) Python objects.
        attr_dict : dictionary, optional (default= no attributes)
            Dictionary of edge attributes.  Key/value pairs will
            update existing data associated with the edge.
        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        See Also
        --------
        add_edge() in class DiGraph

        """
        # set up attribute dict
        if attr_dict is None:
            attr_dict=Counter(attr)
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(\
                    "The attr_dict argument must be a dictionary.")
        # add nodes
        if u not in self.succ:
            self.succ[u]={}
            self.pred[u]={}
            self.node[u] = {}
        if v not in self.succ:
            self.succ[v]={}
            self.pred[v]={}
            self.node[v] = {}
        # add the edge
        datadict=self.adj[u].get(v,{})
        datadict.update(attr_dict)
        self.succ[u][v]=datadict
        self.pred[v][u]=datadict

    def tree_data2json(self, save=False, path2file=None):
        """Return taxonmic TaxGraph in json format for use elsewhere.

        Parameter:
        save: boolean, optional, dafault=False
                    False: pretty print json to stdout
                    True: save json to file given in path2file
        path2file: path, optional, default=None
                    Path for saving json data to.

        """
        data = json_graph.tree_data(self, root="root")
        if save == False:
            print( json.dumps(data, indent = 2))
        elif save == True and path2file is not None:
            with open(path2file, "w") as f:
                json.dump(data, f)
        else:
            print( "Please check parameters for tree_data2json().")


if __name__ == '__main__':
    nodes = ["a", "b", "c"]
    tg = TaxGraph()
    tg.add_edges_from_ordered_nodes(nodes)
    print( tg.degree(weight="weight"))
