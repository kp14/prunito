.. _result_model:

Dealing with results from web services
======================================

Internally, ``prunito`` uses the wonderful `requests <http://docs.python-requests.org/en/master/>`_
package and that provides a sophisticated ``Response`` object.
To avoid losing the power of that object, ``prunito`` provides a wrapper around it,
the ``WSResponse`` class (where WS stands for web service).
Or course, format and content of results returned by calls to web services
depend on the service and so there are custom subclasses of ``WSResponse``
for various web services.

What does WSResponse (sub)classes do?
-------------------------------------

They provide helper methods.
For example, ``size()`` returns the number of hits if that particular web service
gives that information or it can be quickly determined.
``as_file_obejct()`` wraps the result's text representation in :class:`io.StringIO`.
As many results are sequences-like--we do expect one or more hits for a call after all--
special methods for sequences are provided wherever possible to allow iteration or
slicing.
Obviously, this depends on the service and kind of results.

Example: WSResponseUniprot class
--------------------------------

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

