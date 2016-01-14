# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:09:22 2015

@author: kp14
"""
import logging

from lxml import etree
from lxml import objectify


logging.basicConfig(level=logging.WARN)


NS = '{http://uniprot.org/unirule-1.0}'
NS_uniprot = '{http://uniprot.org/uniprot}'


class UniRule():
    '''Class representing a rule as used by the UniRule system.'''

    def __init__(self):
        self.meta = {}
        self.main = BasicRule()
        self.cases = []

    def has_cases(self):
        '''Return whether rule has cases.'''
        return bool(self.cases)

    def has_ec(self):
        '''Return whether rule has EC numbers either in main or in cases.'''
        has_ec = False
        for annot in self.iter_annotations():
            has_ec = self._is_enzyme_number(annot)
            if has_ec:
                break
        return has_ec

    def is_rulebase(self):
        '''Return whether rule comes from the Rulebase pipeline.'''
        return self._is_from_pipeline('rulebase')

    def is_hamap(self):
        '''Return whether rule comes from the HAMAP pipeline.'''
        return self._is_from_pipeline('hamap')

    def is_pir(self):
        '''Return whether rule comes from PIR`s pipeline.'''
        return self._is_from_pipeline('pir')

    def get_ec(self):
        '''Return EC numbers in rule.

        Returns:
        list of strings
        '''
        ec_list = []
        for ann in self.iter_annotations():
            if self._is_enzyme_number(ann):
                ec_list.append(ann.value)
        return ec_list

    def iter_conditions(self):
        '''Iterate over all conditions in main and cases.'''
        yield from self.iter_main_conditions()
        yield from self.iter_case_conditions()

    def iter_annotations(self):
        '''Iterate over all annotations in main and cases.'''
        yield from self.iter_main_annotations()
        yield from self.iter_case_annotations()

    def iter_main_conditions(self):
        '''Iterate over conditions in main only.'''
        return self.main.iter_conditions()

    def iter_main_annotations(self):
        '''Iterate over annotations in main only.'''
        return self.main.iter_annotations()

    def iter_case_conditions(self):
        '''Iterate over conditions in cases only.'''
        for c in self.cases:
            yield from c.iter_conditions()

    def iter_case_annotations(self):
        '''Iterate over annotations in cases only.'''
        for c in self.cases:
            yield from c.annotations

    def _is_enzyme_number(self, annot):
        '''Helper method to determine whether an EC number is among the annotations.'''
        return annot.subtype == 'ecNumber'

    def _is_from_pipeline(self, ppln):
        '''Helper method to determine pipeline of origin.'''
        mapper = {'hamap': 'MF',
                  'pir': 'PIR',
                  'rulebase': 'RU'}
        ppln_low = ppln.lower()
        try:
            return self.meta['oldRuleNum'].startswith(mapper[ppln_low])
        except KeyError:
            print('Invalid pipeline. Choose hamap, pir or rulebase.')

    def __str__(self):
        template = ('Rule ID: {0}\n'
                    'Main:\n'
                    'Number of condition sets: {1}\n'
                    'Number of annotations {2}\n'
                    'Number of cases: {3}\n')
        string = template.format(self.meta['id'],
                                 len(self.main.conditions),
                                 len(self.main.annotations),
                                 len(self.cases)
                                 )
        return string


class BasicRule():
    '''Basic components of a rule: a set of condition(s) (sets) and annotations.'''

    def __init__(self):
        self.conditions = []
        self.annotations = []

    def iter_annotations(self):
        '''Iterate over annotations.'''
        yield from self.annotations

    def iter_conditions(self):
        '''Iterate over conditions.'''
        for condition_set in self.conditions:
            yield from condition_set

    def __str__(self):
        template = ('Number of condition sets: {0}\n'
                    'Number of annotations {1}\n')
        string = template.format(len(self.conditions),
                                 len(self.annotations))
        return string


class Condition():
    '''Represents conditions as used in UniRule.

    Conditions are prerequisites for a rule to be applied to a given
    UniProtKB entry. Currently, conditions used are InterPro signatures,
    taxonomic nodes, proteome properties, sequence flags and length.
    '''

    def __init__(self):
        self.type = None
        self.negative = False
        self.value = None

    def __str__(self):
        return '{0}, {1}, {2}'.format(self.type, self.value, self.negative)


class Annotation():
    '''Annotations as used in UniRule.

    Annotations are applied once conditions are met by a UniProtKB entry.
    Many types of annotations are used. In accordance with the UniProtKB
    flat file format, the annotations have a class equivalent to their
    line type, an optional type like FUNCTION and subtype like FullName.
    '''

    def __init__(self):
        self.class_ = None
        self.type = None
        self.subtype = None
        self.value = None

    def __str__(self):
        return self.__dict__.__str__()


def parse_rules(filename):
    '''Extract rule information from UniRules in a file.'''
    with open(filename, 'rb') as data:
        logging.info('Starting work on file: {}'.format(filename))
        xml = data.read()
        root = objectify.fromstring(xml)
        objectify.deannotate(root, cleanup_namespaces=True)
        for rule in root.unirule:
            rule_id = rule.attrib['id']
            logging.info('Parsing rule: {}'.format(rule_id))
            uni = UniRule()
            logging.info('Extracting meta from : {}'.format(rule_id))
            extract_meta(rule, uni)
            logging.info('Extracting conditions from : {}'.format(rule_id))
            extract_main_conditions(rule, uni)
            logging.info('Extracting annotations from : {}'.format(rule_id))
            extract_main_annotations(rule, uni)
            try:
                for case in rule.cases.case:
                    logging.info('Found a case.')
                    basic_rule = BasicRule()
                    uni.cases.append(basic_rule)
                    extract_case_conditions(case, uni)
                    extract_case_annotations(case, uni)
            except AttributeError:
                logging.info('Rule appears to have no cases: {}'.format(rule_id))
            yield uni

def extract_meta(rule, uni):
    uni.meta.update(rule.attrib)
    for info in rule.information.getchildren():
        try:
            key = info.tag.replace(NS, '')
            val = info.text
            uni.meta[key] = val
        except:
            print('Error in: {}'.format(info))


def extract_main_conditions(rule, uni):
    for c_set in rule.main.conditionSets.conditionSet:
        condition_list = _extract_conditions(c_set)
        uni.main.conditions.append(condition_list)


def extract_case_conditions(case, uni):
    for c_set in case.conditionSets.conditionSet:
        condition_list = _extract_conditions(c_set)
        logging.info('Extracted condition list from case: {}'.format(condition_list))
        uni.cases[-1].conditions.append(condition_list)


def _extract_conditions(rule_element):
    c_list = []
    for child in rule_element.getchildren():
        cond = Condition()
        cond.type = child.attrib['type']
        cond.negative = child.attrib['negative']
        try:
            cond.value = child.value.text
            if child.value.attrib:
                for key, val in child.value.attrib.items():
                    setattr(cond, key, val)
        except AttributeError:
            pass
        try:
            cond.range = child.range
            start = child.range.attrib['start']
            end = child.range.attrib['end']
            cond.start = start
            cond.end = end
            string_ = ': start:{0} end:{1}'.format(start, end)
            if cond.value:
                cond.value += string_
            else:
                cond.value = string_
        except AttributeError:
            pass
        c_list.append(cond)
    return c_list


def extract_main_annotations(rule, uni):
    try:
        for annot in rule.main.annotations.annotation:
            uni.main.annotations.extend(_extract_annotations(annot))
            logging.info('Extracting annotations from main in rule: {}'.format(rule.attrib['id']))
    except AttributeError:
        logging.warn('Rule appears to have no annotations in main: {}'.format(rule.attrib['id']))


def extract_case_annotations(case, uni):
    annotation_list = []
    try:
        for annot in case.annotations.annotation:
            annotation_list.extend(_extract_annotations(annot))
    except AttributeError:
        logging.warn('Case appears to have no annotations: {}'.format(rule.attrib['id']))
    uni.cases[-1].annotations.extend(annotation_list)


def _extract_annotations(annotation_element):
    annotation_list = []
    class_element = annotation_element.getchildren()[0] # Only one toplevel element
    class_ = class_element.tag.replace(NS, '')
    logging.info('Parsing class: {}'.format(class_))
    if class_ == 'comment':
        attribs = class_element.attrib
        if attribs['type'] not in ['subcellular location', 'cofactor']:
            attribs['value'] = class_element.getchildren()[0].text
            attribs['class'] = class_
            annotation_list.append(_create_annotation(attribs))
        elif attribs['type'] == 'subcellular location':
            for location in class_element.getchildren():
                for loc in location.getchildren():
                    attribs[loc.tag.replace(NS_uniprot, '')] = loc.text
                sub_cell_components = []
                for k in ['location', 'topology', 'orientation']:
                    try:
                        sub_cell_components.append(attribs[k])
                    except KeyError:
                        pass
                attribs['value'] = ' / '.join(sub_cell_components)
                attribs['class_'] = class_
                annotation_list.append(_create_annotation(attribs))
        elif attribs['type'] == 'cofactor':
            for cfac in class_element.getchildren():
                if 'cofactor' in cfac.tag:
                    for data in cfac.getchildren():
                        if 'name' in data.tag:
                            attribs['value'] = data.text
                        elif 'dbReference' in data.tag:
                            attribs['id'] = data.attrib['id']
                    attribs['class_'] = class_
                    annotation_list.append(_create_annotation(attribs))
                elif 'text' in cfac.tag:
                    attribs.clear()
                    attribs['class_'] = class_
                    attribs['type'] = 'cofactor'
                    attribs['note'] = cfac.text
                    annotation_list.append(_create_annotation(attribs))

    elif class_ == 'keyword':
        attribs = class_element.attrib
        attribs['value'] = class_element.text
        attribs['class_'] = class_
        annotation_list.append(_create_annotation(attribs))
    elif class_ == 'gene':
        gene_ele = class_element.getchildren()[0]
        attribs = gene_ele.attrib
        attribs['value'] = gene_ele.text
        attribs['class_'] = class_
        annotation_list.append(_create_annotation(attribs))
    elif class_ == 'protein':
        for typ in class_element.getchildren():
            typ_ = typ.tag.replace(NS, '')
            if typ_ == 'alternativeName':
                attribs = {}
                attribs['type'] = typ_
                attribs['subtype'] = 'fullName'
                attribs['class_'] = class_
                attribs['value'] = typ.fullName.text
                annotation_list.append(_create_annotation(attribs))
            elif typ_ =='flag':
                attribs = {}
                attribs['type'] = typ_
                attribs['value'] = typ.value.text
                attribs['class_'] = class_
                annotation_list.append(_create_annotation(attribs))
            elif typ_ == 'recommendedName':
                for subtyp in typ.getchildren():
                    attribs = {}
                    attribs['type'] = typ_
                    attribs['subtype'] = subtyp.tag.replace(NS, '')
                    attribs['value'] = subtyp.text
                    attribs['class_'] = class_
                    annotation_list.append(_create_annotation(attribs))
    return annotation_list


def _create_annotation(adict):
    annotation = Annotation()
    annotation.__dict__.update(**adict)
    return annotation
