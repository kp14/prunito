"""Module for representing UniProtKB data as single annotations.

Many of the fields used in UniProtKB are free text and/or multi-labeled
with many evidence tags attached to one data point. Sometimes it would
be beneficial to be able to treat all data points as separate annotations,
i.e., a data point reported by two different papers would be two separate
annotations. I call these atomic annotations.
"""

import hashlib
import logging
import re
from collections import defaultdict
import peewee
from .parser_knowledgebase_txt import parse_txt
from ...utils import UNIPROT_EVIDENCE_REGEX, PUBMED_REGEX


ECO_PTRN = re.compile('ECO:[0-9]{7}')

db = peewee.SqliteDatabase(':memory:')

class MyBaseModel(peewee.Model):
    class Meta:
        database = db


class Entity(MyBaseModel):
    identifier = peewee.CharField(max_length=50, unique=True)


class Statement(MyBaseModel):
    type = peewee.CharField(max_length=50)
    value = peewee.BlobField()
    start = peewee.IntegerField(default=0)
    end = peewee.IntegerField(default=0)
    # checksum = peewee.CharField(max_length=100, unique=True)

    class Meta:
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('type', 'value', 'start', 'end'), True),
        )

    def __unicode__(self):
        return '{}-{}'.format(self.type, self.value)


class Evidence(MyBaseModel):
    code = peewee.CharField(max_length=11, default='ECO:0000000')
    source = peewee.CharField(max_length=50, default='')
    # checksum = peewee.CharField(max_length=100, unique=True, default=_checksum_ev)

    class Meta:
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('code', 'source'), True),
        )

    def is_sim(self):
        return self.code == 'ECO:0000250'

    def is_exp(self):
        return self.code == 'ECO:0000269'

    def _is_auto(self):
        return self.code in ['ECO:0000256', 'ECO:0000259', 'ECO:0000244']

    def __unicode__(self):
        return '{}-{}'.format(str(self.code), str(self.source))


class Annotation(MyBaseModel):
    entity = peewee.ForeignKeyField(Entity, backref='annotations')
    statement = peewee.ForeignKeyField(Statement)
    evidence = peewee.ForeignKeyField(Evidence)

    class Meta:
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('entity', 'statement', 'evidence'), True),
        )

    def value(self):
        return self.statement.value

    def type(self):
        return self.statement.field_type

    def evidence_code(self):
        return self.evidence.code

    def source(self):
        return self.evidence.source

    def __unicode__(self):
        return '{0}, {1}, {2}'.format(self.entity.identifier,
                                      self.statement.__unicode__(),
                                      self.evidence.__unicode__())


def parse_txt_atomic(source):
    """Parse atomic annotations from the UniProtKB flat file format.

    An atomic annotation consists of a statement associated with an
    entity and some evidence. All the line types currently found in
    UniProtKB entries are  parsed, including ** comments which are
    sometimes found.

    source: source containing one or more UniProtKB
            entries. Can be a file object, file name/path string,
            a pathlib.Path instance or a WSResponse object, i.e.,
            the result of a prunito.uniprot.search() call can be
            fed directly into the parser

    Returns:


    """
    db.connect()
    db.create_tables([Entity, Statement, Evidence, Annotation])
    for entry in parse_txt(source):
        entity, _ = Entity.get_or_create(identifier=entry.primary_accession)
        for comment in entry.comments:
            typ, value = comment.split(': ', maxsplit=1)
            parser2use = _get_specific_parser(comment)
            try:
                for annotation in parser2use(typ, value, entity):
                    with db.atomic() as transaction:
                        try:
                            annotation.save()
                        except peewee.IntegrityError as e:
                            logging.error(e)
                            logging.error(annotation)
                            transaction.rollback()
            except TypeError:
                pass
        for feature in entry.features:
            for annotation in _parse_feature(feature, entity):
                with db.atomic() as transaction:
                    try:
                        annotation.save()
                    except peewee.IntegrityError:
                        logging.error('Integrity error.')
                        logging.error(annotation)
                        transaction.rollback()


def _get_specific_parser(comment):
    typ, value = comment.split(': ', maxsplit=1)
    parser_func = mapper.get(typ[:6].lower(), _parse_freetext)
    return parser_func


def _parse_feature(feature, entity):
    typ, start, end, description, _ = feature
    try:
        start = int(start)
    except ValueError:
        start = re.sub("[^0-9]", "", start) # >1 etc found sometimes
    try:
        end = int(end)
    except ValueError:
        end = re.sub("[^0-9]", "", end)
    try:
        text, evs = description.split('. {')
    except ValueError:
        if description.startswith('{'):
            text = ''
            evs = description[1:]
        else:
            text = description
            evs = None
    value = typ + ' ' + text
    statement , _ = Statement.get_or_create(type=typ, value=value, start=start, end=end)
    if evs:
        evidences = []
        for token in evs.rstrip('}.').split(', '):
            try:
                code, source = token.split('|')
            except ValueError:
                code, source = token, ''
            ev, _ = Evidence.get_or_create(code=code, source=source)
            evidences.append(ev)
        for ev in evidences:
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)
    else:
        ev, _ = Evidence.get_or_create()
        yield Annotation(entity=entity,
                         statement=statement,
                         evidence=ev)


def _parse_freetext(typ, value, entity):
    """Extract Annotations from freetext comments.

    Args:
        typ (str): type of UniProt comment
        value (str): freetext body of a UniProt comment

    Return:
          Annotation instances
    """
    evidences = _extract_evidences(value)
    # There might not be any tags, we provide a token one
    if not evidences:
        ev, _ = Evidence.get_or_create()
        evidences.append(ev)
    # handling body of annotations
    body = value.split(' {')[0]
    statements = re.split('\. ', body)
    from_pubmed = []
    from_sim = []
    from_unclear = []
    for statement in statements:
            if 'PubMed:' in statement:
                from_pubmed.append(statement)
            elif '(By similarity)' in statement:
                from_sim.append(statement)
            elif '(Probable)' in statement:
                from_pubmed.append(statement)
            else:
                from_unclear.append(statement)
    for s in from_pubmed:
        text = s.split(' (PubMed:')[0]
        statement, _ = Statement.get_or_create(type=typ, value=text)
        for ev in evidences:
            if ev.source in s:
                yield Annotation(entity=entity,
                                 statement=statement,
                                 evidence=ev)
    for s in from_sim:
        text = s.split(' (By sim')[0]
        statement, _ = Statement.get_or_create(type=typ, value=text)
        sim_tags = [tag for tag in evidences if tag.is_sim()]
        if len(sim_tags) == 1:
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=sim_tags[0])
        else:
            ev, _ = Evidence.get_or_create(code='ECO:0000250', source='')
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)
    for s in from_unclear:
        text = s.rstrip('. ')
        statement, _ = Statement.get_or_create(type=typ, value=text)
        if not from_sim and not from_pubmed and len(evidences) == 1:# 1 tag applicable to all
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=evidences[0])
        elif len(evidences) > 1 and len(set([e.code for e in evidences])) == 1:# 1 type of tag applicable to all but unclear source
            this_code = evidences[0].code
            ev, _ = Evidence.get_or_create(code=this_code, source='')
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)
        elif len(evidences) > 1 and len(set([e.code for e in evidences])) > 1 and from_sim:
            code = list(set([e.code for e in evidences if not e.code == 'ECO:0000250']))[0]# a hack
            ev, _ = Evidence.get_or_create(code=code, source='')
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)
        else:
            ev, _ = Evidence.get_or_create()
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)


def _extract_evidences(blob):
    evidences = []
    for match in re.finditer(UNIPROT_EVIDENCE_REGEX, blob):
        ev_string = match.group()
        try:
            code, source = ev_string.split('|')
        except ValueError:
            code, source = ev_string, ''
        ev, _ = Evidence.get_or_create(code=code, source=source)
        evidences.append(ev)
    return evidences


def _automatic_tags(list_of_evidences):
    auto_only = [tag for tag in list_of_evidences if tag.is_auto()]
    return auto_only


def _sim_tags(list_of_evidences):
    sim_only = [tag for tag in list_of_evidences if tag.is_sim()]
    return len(sim_only)


def _exp_or_opi_tags(list_of_evidences):
    exp_only = [tag for tag in list_of_evidences if tag.is_exp()]
    return len(exp_only)


def _find_tag_for_pubmed(pubmed, list_of_evidences ):
    for tag in list_of_evidences:
        if tag.source == pubmed:
            return tag
            break

    # body_and_ev = value.split(' {')
    # try:
    #     body, ev = body_and_ev
    # except ValueError:
    #     print('Weird splitting pattern for comment: {} {}'.format(typ, value))
    # else:
    #     # handle evidences
    #     evidences = []
    #     for token in ev.rstrip('}.').split(', '):
    #         try:
    #             code, source = token.split('|')
    #         except ValueError:
    #             code, source = token, None
    #         evidences.append(Evidence(code=code, source=source))
    #     # handle statements
    #     stmts = re.split('\. ', body)
    #     for stmt in stmts:
    #         text = re.split('\(PubMed:', stmt, 1)[0]
    #         if len(evidences) == 1:
    #             anno = Annotation(self.entry.primary_accession,
    #                               Statement(text, typ),
    #                               evidence=ev)
    #             yield anno
    #         else:
    #             for ev in evidences:
    #                 if ev.source in stmt:
    #                     anno = Annotation(self.entry.primary_accession,
    #                                       Statement(text, typ),
    #                                       evidence=ev)
    #                     yield anno


def _parse_subcellular_location(typ, value, entity):
    """Extract Annotations from subcellular location comments.

    Args:
        typ (str): type of UniProt comment
        value (str): body of a UniProt comment

    Return:
          Annotation instances
    """
    note = '' #TODO: Handle the Note
    try:
        locations_all, note = value.split('. Note=')
    except ValueError: # there is no Note
        locations = value.split('. ')
    else:
        locations = locations_all.split('. ')
    for location in locations:
        try:
            loc, evs = location.split(' {')
        except ValueError:
            loc = location
            evs = None
        statement, _ = Statement.get_or_create(type=typ, value=loc)
        if evs:
            evidences = []
            for token in evs.rstrip('}.').split(', '):
                try:
                    code, source = token.split('|')
                except ValueError:
                    code, source = token, ''
                ev, _ = Evidence.get_or_create(code=code, source=source)
                evidences.append(ev)
            for ev in evidences:
                yield Annotation(entity=entity,
                                 statement=statement,
                                 evidence=ev)
        else:
            ev, _ = Evidence.get_or_create()
            yield Annotation(entity=entity,
                             statement=statement,
                             evidence=ev)
    if note:
        yield from _parse_freetext('SUBCELLULAR LOCATION', note, entity)


def _parse_cofactor(typ, value, entity):
    """Extract Annotations from cofactor comments.

    Args:
        typ (str): type of UniProt comment
        value (str): body of a UniProt comment

    Return:
          Annotation instances
    """
    pass


def _parse_interaction(typ, value, entity):
    """Extract Annotations from interaction comments.

        Args:
            typ (str): type of UniProt comment
            value (str): body of a UniProt comment

        Return:
              Annotation instances
        """
    pass

mapper = {'intera': _parse_interaction,
          'subcel': _parse_subcellular_location,
          'cofact': _parse_cofactor,
          'domain': _parse_freetext,
          'ptm': _parse_freetext,
          }