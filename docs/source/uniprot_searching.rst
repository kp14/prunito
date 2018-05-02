.. _uniprot_searching:

UniProt REST services
=====================
Searching, retrieving and mapping data
--------------------------------------

Working with |up|_ data means getting the data first and then extracting the relevant bits.
This section deals more with the former.
Services via `uniprot.org <https://www.uniprot.org>`_ allow the following:

* :ref:`searching-uniprot`
* :ref:`mapping-identifiers`
* :ref:`retrieving-batches`
* :ref:`converting-formats`

Common formats for the above REST services are:

* TXT, also called flat file version. The classic format. Contains the most important information. Hard to parse though.
* TAB, tabular format. Data fields (columns) can be specified.
* GFF, for residue-specific annotations, i.e., anything that can be linked to a (range of) amino acid(s).
* XML
* LIST, just accession numbers.
* FASTA, sequences only, in FASTA format.

.. note::
    Note that while the classic REST API does not provide |upkb|_ entries in JSON format
    the |papi|_ does.

.. _searching-uniprot:

Searching UniProtKB
-------------------

The |up|_ website has powerful search functionality supporting both free-text and field-based queries.
To get the most concise and relevant data set, field-based queries (also called advanced search) should be used.
As UniProt's advanced search is RESTful anyway, searching via prunito can be done using queries copied directly
from the website.
This means that complicated queries can be developed and tested on the website first which might come in handy,
especially when learning.
There is a `help page <https://www.uniprot.org/help/advanced_search>`_ listing all the possible fields.

There is really only one function for searching, :py:func:`prunito.uniprot.rest.uniprotapi.search`.

.. code-block:: python

    from prunito import uniprot as up

    result = up.search('name:laccase AND reviewed:yes', frmt='txt)

    for entry in result:
        print(entry.primary_accession, entry.recommended_full_name)

For convenience, methods limiting results to reviewed (i.e. Swiss-Prot) or unreviewed (i.e. TrEMBL) entries can be
used.
Unsurprisingly, these are :py:func:`prunito.uniprot.rest.uniprotapi.search_reviewed` and
:py:func:`prunito.uniprot.rest.uniprotapi.search_unreviewed`.



.. _mapping-identifiers:

Mapping identifiers
--------------------

.. _converting-formats:

Converting between different UniProt formats
--------------------------------------------

.. _retrieving-batches:

Retrieving batches of entries
-----------------------------

, i.e., without searching but with a list of identifiers.