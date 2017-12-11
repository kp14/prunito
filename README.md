# Prunito

A package providing tools for working with of proteins sequences and associated data.
The focus is on data from UniProtKB and this is reflected in the package name which is an anagram of UniProt.
In addition to UniProt data, a few tools for accessing and/or working data from [EuropePMC](https://europepmc.org/),
[InterPro](https://www.ebi.ac.uk/interpro/) and [ENA](https://www.ebi.ac.uk/ena) are provided.

## UniProt

This includes a parser for the UniProtKB text format. [Biopython](http://biopython.org/) provides a parser as well but writing my own was a way
of learning as much as being able to customize certain aspects. The returned Record objects, for example, provide many
convenience methods for accessing fields like the recommended name, the primary accession etc. 

Also included is a parser for the [UniRule](http://www.uniprot.org/help/unirule) XML format used by some of the automatic annotation pipelines feeding into especially the unreviewed
part of UniProtKB (aka UniProtKB/TrEMBL).

The [Proteins API](https://www.ebi.ac.uk/proteins/api/doc/), which serves UniProt-associated data but is hosted under an EBI URL, is used for retrieving data
on taxonomy nodes.

## Eureopean Nucleotide Archive (ENA)

Retrieve CDS in FASTA format based on protein ID (pid).

## InterPro

Provides a function of comparing UniProtKB coverage of [InterPro](https://www.ebi.ac.uk/interpro/) signatures. Currently broken.