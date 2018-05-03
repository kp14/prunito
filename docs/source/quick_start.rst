.. _quick_start:

Quick start
==========

First, import the relevant package:

.. code:: python

    from prunito import uniprot as up

A simple search for all reviewed UniProtKB entries with *tax* in their names.

.. code:: python

    result = up.search_reviewed('name:tax')

Check how many hits this search retrieved.
As results are sequence-like, they support ``__len__``.

.. code:: python

    result.size()
    # or
    #len(result)

What would happen if a given query had many hits?
By default, the maximum number of hits retrieved is 2000.
This can be changed using the parameter ``limit``.

Let's re-run the previous search but this time not just for reviewed entries but all of UniProtKB:

.. code:: python

    huge = up.search('name:tax', limit=1000)


.. parsed-literal::

    Partial dataset retrieved. Size: 1890. Retrieved: 1000.
    Consider increasing the limit and/or using offset.


If the set limit is lower than the actual number of search hits, the above hint is printed.
One could then set ``limit=2000``.

Back to the initial result.
What does ``size`` mean?
UniProtKB entries.
And these entries we would like to parse.
As prunito provides functionality for both searching and parsing UniProt, one can
directly iterate over the entries in a search result for convenience:

.. code:: python

    entries = list(result)
    # or
    # for entry in result: ...

Iterate over the entries, printing out primary accessions and recommended full names.
Both fields are provided for convenience.

.. code:: python

    for entry in entries:
        print(entry.primary_accession, entry.recommended_full_name)


.. parsed-literal::

    Q06507 Cyclic AMP-dependent transcription factor ATF-4
    P18848 Cyclic AMP-dependent transcription factor ATF-4
    Q9Y6D9 Mitotic spindle assembly checkpoint protein MAD1
    P47911 60S ribosomal protein L6
    P10070 Zinc finger protein GLI2
    Q0VGT2 Zinc finger protein GLI2
    ...


Which methods and fields are available on a Record object?
Basically, all the ones Biopython's REcord objects provide plus as few more for convenience.
See :py:class:`prunito.uniprot.parsers.parser_knowledgebase_txt.Record`.

Get isoforms for those entries that have them.
We use the presence of a keyword, *Alternative splicing*, as a filter here.

.. code:: python

    for e in entries:
        if 'Alternative splicing' in e.keywords:
            for i in e.isoforms():
                print(i)


.. parsed-literal::

    >sp|Q9Y6D9-2|MD1L1_HUMAN Isoform 2 of Mitotic spindle assembly checkpoint protein MAD1 OS=Homo sapiens (Human). OX=['9606']
    MLPARGCVRKRTVWPRLARVLIVTLLTLELSYAPLPCQLSGVPYNTGDPVGRWARPCIWP
    CPWHTTINALKGRISELQWSVMDQEMRVKRLESEKQELQEQLDLQHKKCQEANQKIQELQ
    ...
    >sp|P10070-1|GLI2_HUMAN Isoform 1 of Zinc finger protein GLI2 OS=Homo sapiens (Human). OX=['9606']
    MALTSINATPTQLSSSSNCLSDTNQNKQSSESAVSSTVNPVAIHKRSKVKTEPEGLRPAS
    PLALTQGQVSGHGSCGCALPLSQEQLADLKEDLDRDDCKQEAEVVIYETNCHWEDCTKEY
    ...
    >sp|P10070-2|GLI2_HUMAN Isoform 2 of Zinc finger protein GLI2 OS=Homo sapiens (Human). OX=['9606']
    MALTSINATPTQLSSSSNCLSDTNQNKQSSESAVSSTVNPVAIHKRSKVKTEPEGLRPAS
    PLALTQEQLADLKEDLDRDDCKQEAEVVIYETNCHWEDCTKEYDTQEQLVHHINNEHIHG
    ...


We would like to run a FASTA similarity search against Swiss-Prot for one of the sequences.
Let's take the canonical sequence of the first entry in *entries*.

Here we use the ``ebiwebservices`` module from ``prunito``.
The EBI web services require an email address to be set.

.. code:: python

    from prunito import ebiwebservices as ews

    ews.set_email('some@gmx.de')

    first_entry = entries[0]
    similar = ews.fasta_search(first_entry.as_fasta())

    print(similar.text[:600])

.. parsed-literal::

    # /nfs/public/release/wp-jdispatcher/latest/appbin/linux-x86_64/fasta-36.3.7b/fasta36 -l /nfs/public/ro/es/data/idata/latest/fastacfg/fasta3db -L -T 8 -p -m "F9 fasta-R20180501-155642-0060-16766253-p1m.m9" @:1- +uniprotkb_swissprot+
    FASTA searches a protein or DNA sequence data bank
     version 36.3.7b Jun, 2015(preload9)
    Please cite:
     W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

    Query: @
      1>>>sp|Q06507|ATF4_MOUSE Cyclic AMP-dependent transcription factor ATF-4 OS=Mus musculus (Mouse). OX=['10090'] - 349 aa
    Library: UniProtKB/Swiss-Prot
      199856860 residues in 557275 sequences

    Statistic...


How about using InterPro's HMMER search instead of FASTA?

.. code:: python

    from prunito import interpro as ip

    ip_similar = ip.search_phmmer(first_entry.as_fasta())
    print(ip_similar.summary())

.. parsed-literal::

    acc2	acc	desc	species	kg	evalue
    Q06507	ATF4_MOUSE	Cyclic AMP-dependent transcription factor ATF-4	Mus musculus	Eukaryota	1.0e-232
    Q9ES19	ATF4_RAT	Cyclic AMP-dependent transcription factor ATF-4	Rattus norvegicus	Eukaryota	2.6e-216
    P18848	ATF4_HUMAN	Cyclic AMP-dependent transcription factor ATF-4	Homo sapiens	Eukaryota	2.4e-195
    Q3ZCH6	ATF4_BOVIN	Cyclic AMP-dependent transcription factor ATF-4	Bos taurus	Eukaryota	1.9e-169
    Q6NW59	ATF4_DANRE	Cyclic AMP-dependent transcription factor ATF-4	Danio rerio	Eukaryota	5.0e-34
    Q9Y2D1	ATF5_HUMAN	Cyclic AMP-dependent transcription factor ATF-5	Homo sapiens	Eukaryota	6.3e-20
    Q6P788	ATF5_RAT	Cyclic AMP-dependent transcription factor ATF-5	Rattus norvegicus	Eukaryota	5.8e-18
    Q9GPH3	ATFC_BOMMO	Activating transcription factor of chaperone	Bombyx mori	Eukaryota	2.4e-16
    O70191	ATF5_MOUSE	Cyclic AMP-dependent transcription factor ATF-5	Mus musculus	Eukaryota	3.5e-13
    Q8TFF3	HAC1_HYPJE	Transcriptional activator hac1	Hypocrea jecorina (strain QM6a)	Eukaryota	5.4e-05

The result summary is also available as a dataframe if ``pandas`` is.

.. code:: python

    df_hmmer = ip_similar.as_dataframe()


Do some of the entries contain the same PubMed IDs?
Let's find the 5 most common ones.

.. code:: python

    from collections import Counter

    c = Counter()
    for e in entries:
        c.update(e.all_pubmed_ids)
    print(c.most_common(5))

.. parsed-literal::

    [('15489334', 24), ('20068231', 9), ('14702039', 8), ('23186163', 8), ('21269460', 7)]


Which are the accession numbers and species of those 24 entries containing the most common one (15489334)?

.. code:: python

    for e in entries:
        if '15489334' in e.all_pubmed_ids:
            print(e.primary_accession, e.organism)


.. parsed-literal::

    Q06507 Mus musculus (Mouse).
    P18848 Homo sapiens (Human).
    Q9Y6D9 Homo sapiens (Human).
    P47911 Mus musculus (Mouse).
    Q0VGT2 Mus musculus (Mouse).
    ...


So, which paper is hiding behind this PMID 15489334?
Here we use another module for accessing `EuropePMC <https://europepmc.org>` from ``prunito``.
EuropePMC returns data for example in JSON format.
We can iterate over the results.

.. code:: python

    from prunito import europepmc as epmc

    paper = epmc.get_pmid_metadata('15489334')
    for p in paper:
        print(p['title'])
        print(p['abstractText'])

.. parsed-literal::

    The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).

    "The National Institutes of Health's Mammalian Gene Collection (MGC) project was designed to generate and
    sequence a publicly accessible cDNA resource containing a complete open reading frame (ORF) for every human
    and mouse gene. The project initially used a random strategy to select clones from a large number of cDNA
    libraries from diverse tissues. Candidate clones were chosen based on 5'-EST sequences, and then fully sequenced
    to high accuracy and analyzed by algorithms developed for this project. Currently, more than 11,000 human and
    10,000 mouse genes are represented in MGC by at least one clone with a full ORF. The random selection approach
    is now reaching a saturation point, and a transition to protocols targeted at the missing transcripts is now
    required to complete the mouse and human collections. Comparison of the sequence of the MGC clones to reference
    genome sequences reveals that most cDNA clones are of very high sequence quality, although it is likely that some
    cDNAs may carry missense variants as a consequence of experimental artifact, such as PCR, cloning, or reverse
    transcriptase errors. Recently, a rat cDNA component was added to the project, and ongoing frog (Xenopus) and
    zebrafish (Danio) cDNA projects were expanded to take advantage of the high-throughput MGC pipeline."



The paper mentions the Mammalian Gene Collection.
Why not search EuropePMC for articles mentioning the collection in their abstracts?

.. code:: python

    mgc_papers = epmc.search('abstract:"Mammalian Gene Collection"')
    mgc_papers.size()
    #
    # len(mgc_papers)
    for idx, hit in enumerate(mgc_papers):
        print(idx, hit['title'])


.. parsed-literal::

    0 Identification of candidate transcription factor binding sites in the cattle genome.
    1 Selenoproteins in bladder cancer.
    2 NSrp70 is a novel nuclear speckle-related protein that modulates alternative pre-mRNA splicing in vivo.
    3 Generation of a genome scale lentiviral vector library for EF1Î± promoter-driven expression of human ORFs ...
    4 The completion of the Mammalian Gene Collection (MGC).
    5 A high-throughput platform for lentiviral overexpression screening of the human ORFeome.
    6 PRFdb: a database of computationally predicted eukaryotic programmed -1 ribosomal frameshift signals.
    7 Transcriptome analysis of a cDNA library from adult human epididymis.
    ...


Each hit/paper has many extra data fields including DOI, PubMed ID etc.
If the abstract is needed, ``resulttype='core'`` has to be specified as a search parameter.

.. code:: python

    for k, v in list(mgc_papers)[3].items():
        print(k + ':\t' + str(v))


.. parsed-literal::

    id:	23251614
    source:	MED
    pmid:	23251614
    pmcid:	PMC3520899
    doi:	10.1371/journal.pone.0051733
    title:	Generation of a genome scale lentiviral vector library for EF1Î± promoter-driven expression of human ORFs and identification of human genes affecting viral titer.
    authorString:	Å kalamera D, Dahmer M, Purdon AS, Wilson BM, Ranall MV, Blumenthal A, Gabrielli B, Gonda TJ.
    journalTitle:	PLoS One
    issue:	12
    journalVolume:	7
    pubYear:	2012
    journalIssn:	1932-6203
    pageInfo:	e51733
    pubType:	research support, non-u.s. gov't; research-article; journal article;
    isOpenAccess:	Y
    inEPMC:	Y
    inPMC:	Y
    hasPDF:	Y
    hasBook:	N
    hasSuppl:	Y
    citedByCount:	8
    hasReferences:	Y
    hasTextMinedTerms:	Y
    hasDbCrossReferences:	Y
    dbCrossReferenceList:	{'dbName': ['EMBL']}
    hasLabsLinks:	Y
    hasTMAccessionNumbers:	Y
    tmAccessionTypeList:	{'accessionType': ['gen']}
    firstPublicationDate:	2012-12-12