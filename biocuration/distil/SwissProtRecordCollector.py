from collections import namedtuple
import itertools
import json
import re

from IPython.display import HTML, display
from jinja2 import Template
import networkx as nx

import pointcloud
from .UniProtGraph import TaxGraph


class SwissProtRecordCollector(object):
    """Class to collect and count annotations from SwissProt entries.

    Considered types of annotations are:
                * accessions (AC line)
                * descriptions (DE lines)
                * gene names (GN line)
                * organelles (OG line)
                * organism classification (taxonomy) (OC lines)
                * comments (CC lines):
                  specified in class CommentTypes
                * cross references (DR lines):
                  specified in the class CrossRef
                * keywords (KW lines)
                * features (FT lines)

    """

    def __init__(self):
        """Initialize with all fields empty lists or dicts, apart from taxGraph,
        which is an empty DiGraph

        """
        self.accessions = []
        self.taxGraph = TaxGraph()
        self.graph = nx.Graph()
        self.current_entry = None
        self.context = None
        self.sequences = []
        self.template = Template("""
<head>
    <title>Test</title>
</head>
<body>
    <table>
    {% for item in data %}
        <tr>
        <td>{{ item[0] }}</td>
        <td>{{ item[1] }} </td>
        <td>{{ item[2] }}</td>
        <td>{{ item[3] }}</td>
        </tr>
    {% endfor %}
    </table>
</body>
""")

    def collect_record(self, sp_record):
        """The main work function that collects the parsed annotation from
        each SwissProt record as it comes in.

        """
        self.current_entry = sp_record.accessions[0]
        self.graph.add_node(self.current_entry, {'class': 'ac', 'typ': None, 'group': 1})
        self.__collect_accessions(sp_record.accessions)
        self.__collect_descriptions(sp_record.description)
        self.__collect_gene_names(sp_record.gene_name)
        self.__collect_organelles(sp_record.organelle)
        self.__collect_organism_classifications(sp_record.organism_classification)
        self.__collect_comments(sp_record.comments)
        self.__collect_cross_references(sp_record.cross_references)
        self.__collect_keywords(sp_record.keywords)
        #self.__collect_features(sp_record.features)
        self.__collect_sequences(sp_record)

    def __collect_accessions(self, acc_list):
        """Appends new accessions to list."""
        self.accessions.append(acc_list[0])

    def __collect_comments(self, list_of_comments):
        """Counts individual sentences in comments, grouped by comment type."""
        for comment in list_of_comments:
            for comment_type in CommentTypes.InterestingTypes:
                if comment.startswith(comment_type):
                    self.context = comment_type
                    type_specific_string_list = StringHelper.split_and_strip(comment, splitter=". ")
                    for specific_string in type_specific_string_list:
                        cleaned_string = StringHelper.remove_prefix(comment_type, specific_string)
                        cleaned_string = SPQualifiers.remove_qualifiers(cleaned_string)
                        if cleaned_string:
                            self.add_or_count_key(cleaned_string, 'comment', typ=self.context, group=2)

    def __collect_gene_names(self, gene_name):
        """Counts gene names."""
        self.context = 'gene'
        self.add_or_count_key(gene_name, 'name', typ=self.context, group=3)

    def __collect_organelles(self, organelle):
        """Counts organelles."""
        self.context = 'organelle'
        if organelle != "":
            self.add_or_count_key(organelle, 'organelle', group=10)

    def __collect_organism_classifications(self, a_list):
        """Counts individual taxonomic nodes and adds those nodes to a DiGraph
        in which the hierarchy, i.e. node relations, are maintained.

        """
        for item in a_list:
            self.add_or_count_key(item, 'taxon', typ=None, group=4)
        self.taxGraph.add_edges_from_ordered_nodes(a_list)
        for item in a_list:
            self.taxGraph.add_edge(self.current_entry, item)

    def __collect_descriptions(self, a_string):
        a_list = StringHelper.split_and_strip(a_string, splitter=";")
        for item in a_list[:-1]:
            try:
                prefix, name = item.split('=')
            except ValueError:
                print( "Raised ValueError: {}".format(item))
            if prefix in ["RecName: Full", "AltName: Full", "EC"]:
                self.context = prefix[:7]
                self.add_or_count_key(name, 'name', typ=self.context, group=2)
                words_in_name = name.split()
                if len(words_in_name) > 1:
                    for word in words_in_name:
                        self.add_or_count_key(word.lower(), 'name', typ=self.context, group=2)

    def __collect_cross_references(self, list_of_tuples):
        for item in list_of_tuples:
            if item[0] in CrossRefs.InterestingXRef:
                self.context = item[0]
                self.add_or_count_key(item[1], 'xref', typ=item[0], group=4)

    def __collect_keywords(self, a_list):  # TODO: Have a whitelist or blacklist
        self.context = 'keyword'
        for item in a_list:
            self.add_or_count_key(item, 'keyword', group=5)

    def __collect_features(self, list_of_tuples):  # TODO: Have a whitelist or blacklist
        Feature = namedtuple("Feature", ["key", "description"])
        for item in list_of_tuples:
            f = Feature(key=item[0], description=SPQualifiers.remove_qualifiers(item[3]).rstrip(" ."))
            self.add_or_count_key("{0}-{1}".format(f.key, f.description), 'feature', typ=f.key, group=6)

    def __collect_sequences(self, sp_record):
        self.sequences.append(SeqHelper.to_fasta(sp_record))

    def add_or_count_key(self, a_key, _class, typ=None, group=1):
        if self.graph.has_node(a_key):
            self.graph.node[a_key]['freq'] += 1
        else:
            self.graph.add_node(a_key, {'freq': 1, 'class': _class, 'typ': typ, 'group': group})
        self.graph.add_edge(self.current_entry, a_key)

    def summarize_all(self, cutoff=0.0, save_it=False, plot_it=None):
        print( "\n{0} entries were analyzed\n".format(self.get_number_of_entries()))
        # Let's exclude accession nodes
        try:
            data_nodes = [n for n in self.graph.nodes() if not self.graph.node[n]['class'] == 'ac']
        except KeyError:
            print( "We have a problem extracting the nodes names from the graph.")

        report = []

        for n in data_nodes:
            actual_ratio = self.graph.node[n]['freq'] / float(self.get_number_of_entries())
            if actual_ratio > cutoff:
                report.append([self.graph.node[n]['class'],
                               self.graph.node[n]['typ'],
                               self._construct_url(n, self.graph.node[n]),
                               str(self.graph.node[n]['freq'])])
        # sort on first item in lists, i.e. node class
        report_sorted = sorted(report)
        display(HTML(self.template.render(data=report_sorted)))
        #for item in report_sorted:
            #print( "{0}\t{1}:\t{2} ---> {3}".format(*item))

        if save_it:
            target = "d3_example/force/force.json"
            dmp = nx.readwrite.json_graph.node_link_data(self.graph)
            json.dump(dmp, open(target, 'w'))
            print( "Wrote graph data to json in {}".format(target))
        if plot_it and plot_it == '3d':
            self._make_3D_point_cloud()
        if plot_it and plot_it == '2d':
            self._make_2D_point_cloud()

    def _construct_url(self, node_text, node):
        url_mapping = {"comment": None,
                       "name": None,
                       "feature": None,
                       "organelle": None,
                       "taxon": "http://www.uniprot.org/taxonomy/?query=",
                       "keyword": "http://www.uniprot.org/keywords/?query=",
                       "interpro": "https://www.ebi.ac.uk/interpro/search?q=",
                       "go": "http://www.ebi.ac.uk/QuickGO/GTerm?id="
                       }
        node_class = node["class"]
        if node_class == "xref":
            if node["typ"] == "GO":
                node_class = "go"
            else:
                node_class = "interpro"
        if url_mapping[node_class]:
            url = self._wrap_link(node_text, url_mapping[node_class])
        else:
            url = node_text
        return url

    def _wrap_link(self, node_text, node_class):
        target = "".join([node_class, node_text])
        return "<a href=\"{0}\">{1}<a>".format(target, node_text)


    def summarize_web(self, cutoff=0.0, save_it=False, plot_it=None):
        #print "\n{0} entries were analyzed\n".format(self.get_number_of_entries())
        # Let's exclude accession nodes
        report = {}
        try:
            data_nodes = [n for n in self.graph.nodes() if not self.graph.node[n]['class'] == 'ac']
        except KeyError:
            report['result'] = "We have a problem extracting the nodes names from the graph."
            return report

        result = []

        for n in data_nodes:
            actual_ratio = self.graph.node[n]['freq'] / float(self.get_number_of_entries())
            if actual_ratio > cutoff:
                result.append([self.graph.node[n]['class'],
                               self.graph.node[n]['typ'],
                               n,
                               str(self.graph.node[n]['freq'])])
        # sort on first item in lists, i.e. node class
        result_sorted = sorted(result)
        text = ""
        for item in result_sorted:
            text = text + "\n{0}\t{1}:\t{2} ---> {3}".format(*item)

        report['result'] = text
        return report

        # if save_it:
        #     target = "d3_example/force/force.json"
        #     dmp = nx.readwrite.json_graph.node_link_data(self.graph)
        #     json.dump(dmp, open(target, 'w'))
        #     print "Wrote graph data to json in {}".format(target)
        # if plot_it and plot_it == '3d':
        #     self._make_3D_point_cloud()
        # if plot_it and plot_it == '2d':
        #     self._make_2D_point_cloud()

    def get_number_of_entries(self):
        return len(self.accessions)

    def get_graph(self):
        print( 'Starting to plot...')
        import matplotlib.pyplot as plt
        #composed_graph = nx.compose(self.graph, self.taxGraph)
        composed_graph = self.taxGraph
        nx.draw(composed_graph, pos=nx.spring_layout(composed_graph))
        plt.savefig('graph.pdf')

    def _make_3D_point_cloud(self):
        """Visualize compared entries as points in 3D space.

        Args:
            self, nbunch

        Returns:
            matplotlib display
        """
        try:
            data_nodes = [n for n in self.graph.nodes() if not self.graph.node[n]['class'] == 'ac']
        except KeyError:
            print( "We have a problem extracting the nodes names from the graph.")
        #number of data nodes/vectors
        number_of_vectors = len(data_nodes)

        #calculate one vector for each node
        cart_vectors = pointcloud.get_points_equiangularly_distanced_on_sphere(numberOfPoints=number_of_vectors)

        spher_vectors = []

        for cvec in cart_vectors:
            spher_vectors.append(pointcloud.cartesian_to_spherical(cvec))

        maximum = self.get_number_of_entries()

        for name, coord in itertools.izip(data_nodes, spher_vectors):
            #coord[0] *= pointcloud.remap(self.graph.node[name]['freq'], maximum, slope=0.3)
            coord[0] *= self.graph.node[name]['freq']
            self.graph.node[name]['vector'] = pointcloud.spherical_to_cartesian(coord)

        #get accessions
        acc_nodes = [n for n in self.graph.nodes() if self.graph.node[n]['class'] == 'ac']

        final_coord = []
        #get neighbors of each acc
        for acc in acc_nodes:
            neighbors = self.graph.neighbors(acc)
            coord_product = [0, 0, 0]
            for neighbor in neighbors:
                coord_product = pointcloud.scalar_vector_product(coord_product, self.graph.node[neighbor]['vector'])
            self.graph.node[acc]['vector'] = coord_product
            final_coord.append(coord_product)

        # for c in final_coord:
        #     print c

        xs = [n[0] for n in final_coord]
        ys = [n[1] for n in final_coord]
        zs = [n[2] for n in final_coord]

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xs, ys, zs)
        fig.show()

    def _make_2D_point_cloud(self):
        """Visualize compared entries as points in 2D space.

        Args:
            self

        Returns:
            matplotlib display
        """
        try:
            data_nodes = [n for n in self.graph.nodes() if not self.graph.node[n]['class'] == 'ac']
        except KeyError:
            print( "We have a problem extracting the nodes names from the graph.")
        #number of data nodes/vectors
        number_of_vectors = len(data_nodes)





class CrossRefs(object):
    """
    Class listing the types of cross references to be considered.
    """
    InterestingXRef = {"GO",
                       "HAMAP",
                       "InterPro",
                       "Gene3D",
                       "SUPFAM",
                       "PANTHER",
                       "Pfam",
                       "PIRSF",
                       "PRINTS",
                       "ProDom",
                       "SMART",
                       "TIGRFAMs",
                       "PROSITE"}


class CommentTypes(object):
    """
    Class listing the types of comments to be considered.
    """
    InterestingTypes = {"FUNCTION",
                        "CATALYTIC ACTIVITY",
                        "COFACTOR",
                        "ENZYME REGULATION",
                        "SUBUNIT",
                        "PATHWAY",
                        "SUBCELLULAR LOCATION",
                        "INDUCTION",
                        "DOMAIN",
                        "SIMILARITY"}


class FeatureTypes(object):
    """
    Class listing the types of features to be considered.
    """
    InterestingTypes = {"FUNCTION",
                        "CATALYTIC ACTIVITY",
                        "COFACTOR",
                        "ENZYME REGULATION",
                        "SUBUNIT",
                        "PATHWAY",
                        "SUBCELLULAR LOCATION",
                        "INDUCTION",
                        "DOMAIN",
                        "SIMILARITY"}


class SPQualifiers(object):
    """
    Utility class providing both qualifiers used in Swiss-Prot
    and a static method to remove them from an input string.
    """
    Qualifiers = {"(By similarity)",
                  "(Potential)",
                  "(Probable)"}

    @staticmethod
    def remove_qualifiers(a_string):
#        for qualifier in SPQualifiers.Qualifiers:
        for qualifier in ["(By similarity)", "(Potential)", "(Probable)"]:
            if qualifier in a_string:
                a_string = a_string.replace(qualifier, "").rstrip()
        # we add removal of tags here for now; should refactored at some point.
        no_tags = re.sub('\{ECO.*\}|\(PubMed:.*\)', '', a_string)
        return no_tags


class StringHelper(object):
    """Utility class for string manipulation relevant to Swiss-Prot."""

    @staticmethod
    def remove_empty_strings(a_list):
        """ Deletes empty string from list.
        """
        while '' in a_list:
            a_list.remove('')
        return a_list

    @staticmethod
    def remove_prefix(prefix, a_string):
        """ Removes specified prefix from string.
        """
        if prefix.endswith("="):
            return a_string.replace(prefix, "")
        else:
            prefix += ": " # Comment types are followed by a colon and a space
            return a_string.replace(prefix, "")

    @staticmethod
    def split_and_strip(a_string, splitter=";"):
        """ Splits string at specified token and strips white space and full stops.
        """
        a_list = a_string.split(splitter)
        a_list = [item.strip(' .') for item in a_list]
        return StringHelper.remove_empty_strings(a_list)


class SeqHelper(object):
    """Utility class for sequence manipulation."""

    @staticmethod
    def to_fasta(sp_record):
        """Returns sequence in fasta format as a string."""
        return ">{0}\n{1}".format(sp_record.accessions[0], sp_record.sequence)
