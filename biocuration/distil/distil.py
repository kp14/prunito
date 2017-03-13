import sys
import argparse

# partial imports
from IPython.display import clear_output, display
from ipywidgets import widgets

# import from own modules/packages
from .SwissProtRecordCollector import SwissProtRecordCollector
from ..uniprot import parse_txt_compatible
from ..uniprot import search_all, search_reviewed
from ..utils import Borg


class SPAnalysis(Borg):
    """A class modelling a commonprot analysis."""

    def __init__(self):
        Borg.__init__(self)
        if 'cache' in self.__dict__:
            pass
        else:
            self.cache = {}

    def run_analysis(self, query, reviewed, ratio):
        """Distill data from hits of a UniProtKB search.

        If the combination of query and reviewed has been seen before, results
        are reused and analyzed using the possibly new ratio.

        @query: a UniProtKB query; str
        @reviewed: search Swiss-Prot only; boolean
        @ratio: frequency of reported data points; float
        """
        collector = self._check_whether_cached(query, reviewed)

        if collector:
            self._re_run(collector, ratio)

        else:
            self._run_first_time(query, reviewed, ratio)

    def _check_whether_cached(self, query, reviewed):
        """Check whether a given analysis has been run before."""
        try:
            return self.cache[(query, reviewed)]
        except KeyError:
            return None

    def _re_run(self, collector, ratio):
        """Re-analyze the collected results from a previous analysis
        using a different ratio.
        """
        clear_output()
        print("Re-running analysis with new ratio using cached data.")
        collector.summarize_all(cutoff=ratio, save_it=False, plot_it=None)

    def _run_first_time(self, query, reviewed, ratio):
        """Run an analysis for the first time, i.e., results are not cached."""
        if reviewed:
            results = search_reviewed(query, frmt='txt', file=True)
        else:
            results = search_all(query, frmt='txt', file=True)

        #query_final = query.encode('ascii', 'ignore')

        #results = search_mode(query, frmt='txt', file=True)
        entries = parse_txt_compatible(results)
        if not entries:
            sys.exit('No entries that could be parsed were retrieved. Check your query.')
        collector = SwissProtRecordCollector()
        for entry in entries:
            collector.collect_record(entry)

        clear_output()

        self.cache[(query, reviewed)] = collector

        collector.summarize_all(cutoff=ratio, save_it=False, plot_it=None)


def run_notebook():
    """Create the GUI controls and any control logic."""

    #Create Borg for others to communicate with
    analysis = SPAnalysis()

    # text box for entering UniProt queries
    query = widgets.Text(description="Query")

    # checkbox to select Swiss-Prot only searches
    reviewed_check = widgets.Checkbox(description="Search reviewed",
                                            value=True)

    # slider for the ratio
    ratio_slider = widgets.FloatSlider(description="Ratio",
                                             min=0.0,
                                             max=1.0,
                                             step=0.05,
                                             value=0.9)
    # bundle them up...
    controls = [query,
                reviewed_check,
                ratio_slider,
                ]
    #... and pass them to the container for display
    interface = widgets.Box(children=controls)
    display(interface)

    # Create the buttons...
    run_button = widgets.Button(description="Run", button_style='primary')
    clear_cache_button = widgets.Button(description="Clear cache", button_style='warning')
    button_container = widgets.Box(children=[run_button, clear_cache_button])

    #... and display them; CSS changes have to be run after the display() function
    display(button_container)
    #button_container.remove_class('vbox')
    #button_container.add_class('hbox')

    def on_run_button_clicked(b):
        #This Borg knows what the outer Borg, created
        #when starting the interface, knows. I.e. it
        #knows about previous searches -> caching!
        temp_analysis = SPAnalysis()

        querystr = query.value
        reviewed = reviewed_check.value
        ratio = ratio_slider.value
        temp_analysis.run_analysis(querystr, reviewed, ratio)

    def on_clear_cache_button_clicked(b):
        #This Borg knows what the outer Borg, created
        #when starting the interface, knows.
        #Deleting its cache means deleting the entire Borg cache
        temp_analysis = SPAnalysis()
        clear_output()
        try:
            temp_analysis.cache = {}
            print("Cached searches have been cleared.")
        except AttributeError:
            print("Cache is empty already.")

    run_button.on_click(on_run_button_clicked)
    clear_cache_button.on_click(on_clear_cache_button_clicked)
