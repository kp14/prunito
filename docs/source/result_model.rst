.. _result_model:

Dealing with results from web services
======================================

Internally, ``prunito`` uses the wonderful `requests <http://docs.python-requests.org/en/master/>`_
package and that provides a sophisticated :class:`~requests.Response` object.
To avoid losing the power of that object, ``prunito`` provides a wrapper around it,
the :class:`~prunito.utils.WSResponse` class (where WS stands for web service).
Or course, format and content of results returned by calls to web services
depend on the service and so there are custom subclasses of :class:`~prunito.utils.WSResponse`
for various web services.

What does WSResponse (sub)classes do?
-------------------------------------

They provide helper methods.
For example, ``size()`` returns the number of hits if that particular web service
gives that information or it can be quickly determined.
:meth:`as_file_obejct() <prunito.utils.WSResponse>` wraps the result's text representation
in :class:`io.StringIO`.
As many results are sequences-like--we do expect one or more hits for a call after all--
special methods for sequences are provided wherever possible to allow iteration or
slicing.
Obviously, this depends on the service and kind of results.

Any attribute they don't recognize is passed on to the wrapped Response object.
This means that things like the the Response's attributes can be accessed via as
expected: ``text``, ``url``, ``status``.
Methods like ``json()`` are available, too.

Results from searching |upkb|
-----------------------------

In addition to what the parent class has, this one also:

*   Has a method, ``release()``, for getting the current release.
    Releases are specified as year_number, e.g. 2017_10.
*   Has a method, ``date()``, to get the date of the current release.
    This can be returned as a string or a :class:`datetime.datetime` obejct.
*   Allows iterating over results for some formats.

    When full |upkb| text entries are retrieved, they are parsed and the iterator
    returns the ``Record`` object. For tabular data, like *list*, *gff* or *tab*,
    lines are iterated over.
*   Tabular data can also be returned as a pandas dataframe.

See also :class:`~prunito.uniprot.rest.uniprotapi.WSResponseUniprot`.

Results from identifiers mappings
---------------------------------

Results for mapping calls are returned as a table with two columns, *From* and *To*.
This table can be accessed as text via ``<obj>.text``,
the lines can be iterated over
but, for convenience, a dictionary of the results is prepared as ``<obj>.as_dict()``.
A list of mapped-to target IDs is ``<obj>.target_ids()``.
Identifiers that could not be mapped will be silently ignored,
i.e., there won't be any mappings in the result set.

See also :class:`~prunito.uniprot.rest.uniprotapi.WSResponseUniprotMapping`.

Results from HMMER searches
---------------------------

*todo*
See also :class:`~prunito.utils.WSResponseHmmer`.

Results from searching EuropePMC
--------------------------------

*todo*
See also :class:`~prunito.utils.WSResponseEPMC`.

Taxonomy results
----------------

*todo*
See also :class:`~prunito.utils.WSResponseTax`.