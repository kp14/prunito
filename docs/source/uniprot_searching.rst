.. uniprot_searching:

UniProt REST services
=====================
Searching, retrieving and mapping data
--------------------------------------

Working with UniProt data means getting the data first and then extracting the relevant bits.
This section deals more with the former.
Services via `uniprot.org <https://www.uniprot.org>`_ allow the following:

* Searching the Uniprot Knowledgebase (UniProtKB)
* Retrieving batches of entries, i.e., without searching but with a list of identifiers.
* Converting between database identifiers.
* Converting between different UniProt formats.

Common formats for the above REST services are:

* TXT, also called flat file version. The classic format. Contains the most important information. Hard to parse though.
* TAB, tabular format. Data fields (columns) can be specified.
* GFF, for residue-specific annotations, i.e., anything that can be linked to a (range of) amino acid(s).
* XML
* LIST, just accession numbers.
* FASTA, sequences only, in FASTA format.

Searching UniProtKB
-------------------


