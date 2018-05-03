.. _uniprot_searching:

UniProt REST services
=====================

Using |up| REST services is about retrieving data, for example based on a query.
By carefully choosing their format, retrieved data might constitute the true end result,
they might be input for another procedure (e.g. FASTA sequences feeding into a similarity search),
but often results have to processed further.
This is called *parsing*.
Another section of the help will take a closer look at this.

|up| vs. |papi|
---------------

Services via |up|_ allow the following:

* :ref:`searching-uniprot`
* :ref:`mapping-identifiers`
* :ref:`retrieving-batches`
* :ref:`converting-formats`

The |papi|_ provides access data (sets) that often go beyond what you would find in |upkb|.
Many have not been included in ``prunito`` yet.

* Protein and isoform sequences
* Residue-specific annotations (also called *features*)
* Variation data
* Proteomics data
* Antigen binding sites
* Proteome data
* :ref:`getting-taxonomy`
* Genomic coordinates
* Access to UniParc

Common result formats for the REST services are:

* TXT, also called flat file version. The classic format. Contains the most important information. Hard to parse though.
* TAB, tabular format. Data fields (columns) can be specified.
* GFF, for residue-specific annotations, i.e., anything that can be linked to a (range of) amino acid(s).
* XML
* LIST, just accession numbers.
* FASTA, sequences only, in FASTA format.
* JSON

.. note::
    Note that while the classic |up| REST API does not provide data in JSON format the |papi|_ does.
    On the other hand, |papi| does not provide GFF or tabular output.

.. _searching-uniprot:

Searching UniProtKB
-------------------

The |up|_ website has powerful search functionality supporting both free-text and field-based queries.
To get the most concise and relevant data set, field-based queries (also called advanced search) should be used.
As UniProt's advanced search is RESTful anyway, searching via ``prunito`` can be done using queries copied directly
from the website.
This means that complicated queries can be developed and tested on the website first which might come in handy,
especially when learning.
There is a `help page <https://www.uniprot.org/help/advanced_search>`_ listing all the possible fields for the
advanced search.

There is really only one function for searching, :py:func:`prunito.uniprot.rest.uniprotapi.search`.
A query string is the only mandatory parameter.
In that case full |upkb| entries in text (flat file) format are retrieved.
The default limit is 2000 entries; this can be changed, of course.

For convenience, methods limiting results to reviewed (i.e. Swiss-Prot) or unreviewed (i.e. TrEMBL) entries can be
used.
Unsurprisingly, these are :py:func:`prunito.uniprot.rest.uniprotapi.search_reviewed` and
:py:func:`prunito.uniprot.rest.uniprotapi.search_unreviewed`.

A few example queries:

.. code-block:: python

    from prunito import uniprot as up

    result = up.search('name:laccase AND reviewed:yes')
    # or using the convenience function
    result = up.search_reviewed('name:laccase')
    # getting sequences only
    result = up.search_reviewed('name:laccase', frmt='fasta')
    # As laccases are enzyme, get relevant enzyme data
    result = up.search_reviewed('name:laccase', frmt='tab', columns='id,entry name,comment(FUNCTION),ec')


.. _mapping-identifiers:

Mapping identifiers
--------------------

A paper might contain a list of identifiers for 3D protein structures.
A repository for such structure is PDB and their IDs look like *1ABC*.
Say we wanted to map those PDB IDs to Ensembl ones--this is what the mapping does.
As mappings always have to include |up| accessions as either source or target,
mapping from PDB to Ensembl is s two-step process.

.. code-block:: python

    up_from_pdb = up.map_to_or_from_uniprot(['1abc', '5ukz', '3tui'], 'PDB_ID', 'ACC')
    ensembl = up.map_to_or_from_uniprot(list(up_from_pdb), 'ACC', 'ENSEMBL_ID')

A full list of sources, targets and their abbreviations can be found `here <https://www.uniprot.org/help/api_idmapping>`_.

.. _converting-formats:

Converting between different UniProt formats
--------------------------------------------

I don't think this is used much.
One could, for example, convert the text version of a |up| entry into XML.
The text entry would have to be without any errors though for this to work.

.. _retrieving-batches:

Retrieving batches of entries
-----------------------------

If one already has a list of |up| accessions these can be retrieved using the batch functionality.

.. code-block:: python
    result = up.retrieve_batch(['P12345', 'P12344'], frmt='txt')

.. _getting-taxonomy:

Retrieving taxonomy data
------------------------
