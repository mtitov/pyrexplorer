#
# Copyright 2015 Mikhail Titov
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Authors:
# - Mikhail Titov, <mikhail.titov@cern.ch>, 2015
#

__all__ = ['SPADE']

from collections import defaultdict
from itertools import ifilter

from .element import Element, ElementPool


class SPADE(object):

    def __init__(self):
        """Initialization."""
        self._sequences = {}
        self._minimum_support = None
        self._frequent_elementpool = None

    def set(self, **kwargs):
        """
        Set initial data.

        @param kwargs: Input parameters.
        @type kwargs: dict

        @keyword sequences: Dictionary of sequences {sid: {eid: <itemset>}}.
        @keyword minimum_support: Minimum support (number of distinct sids).
        """
        if isinstance(kwargs.get('sequences'), dict):
            self._sequences = kwargs['sequences']
        self._minimum_support = kwargs.get('minimum_support')
        self._frequent_elementpool = ElementPool()

    def is_maximal(self, element):
        """
        Check that element's sequence is maximal frequent sequence.

        @param element: Element object.
        @type element: Element
        @return: Flag that element doesn't have sub-seq for any frequent seq.
        @rtype: bool
        """
        output = True

        for freq_element in self._frequent_elementpool.itervalues():
            if element.has_subsequence(freq_element):
                output = False
                break

        return output

    def add_elements(self, elements, top_number=None):
        """
        Add elements to the ElementPool of frequent sequences.

        @param elements: List of Element objects.
        @type elements: list
        @param top_number: The number of top longest output sequences.
        @type top_number: int/None
        """
        for element in elements:
            self._frequent_elementpool[element.sequence] |= element
        self.keep_top_elements(top_number=top_number)

    def keep_top_elements(self, top_number):
        """
        Keep only elements with top longest sequential patterns.

        @param top_number: The number of top longest output sequences.
        @type top_number: int
        """
        if top_number and len(self._frequent_elementpool) > top_number:

            top_freq_sequences = set(
                map(lambda x: x.sequence,
                    sorted(self._frequent_elementpool.values(),
                           key=lambda x: (len(x), x.num_itemsets),
                           reverse=True)[:top_number]))

            for seq in self._frequent_elementpool.keys():
                if seq not in top_freq_sequences:
                    del self._frequent_elementpool[seq]

    def generate_frequent_sequences(self, maximal=False, maximum_length=None):
        """
        Compute frequent 1-sequences and 2-sequences.

        @param maximal: Flag to get only maximal frequent sequences.
        @type maximal: bool
        @param maximum_length: The maximum length of sequential patterns.
        @type maximum_length: int/None
        @return: Two lists of frequent 1- and 2-sequences.
        @rtype: tuple(list, list)
        """
        if not self._sequences or not self._minimum_support:
            raise Exception('Initial sequences/support are not set')

        id_lists, sequences = {}, {}
        for sid in self._sequences:
            sequences.setdefault(sid, [])
            for eid in sorted(self._sequences[sid]):
                for item in sorted(set(self._sequences[sid][eid])):
                    id_lists.\
                        setdefault(item, {}).\
                        setdefault(sid, []).\
                        append(eid)
                    sequences[sid].append(item)

        freq_1_seq_elementpool = ElementPool()

        for item in id_lists:
            if len(id_lists[item]) < self._minimum_support:
                continue
            for sid in id_lists[item]:
                for eid in id_lists[item][sid]:
                    freq_1_seq_elementpool[item] |= \
                        Element(item, sid=sid, eid=eid)

        items_pair_frequency = defaultdict(int)
        for items in sequences.itervalues():
            items = filter(lambda x: x in freq_1_seq_elementpool, items)
            for idx_i in range(len(items)):
                for idx_j in range(idx_i + 1, len(items)):
                    items_pair_frequency[(items[idx_i], items[idx_j])] += 1

        freq_2_seq_elementpool = ElementPool()

        if maximum_length is None or maximum_length > 1:
            used_freq_items = set()

            for item_i, item_j in set(
                    map(lambda x: tuple(sorted(x[0])),
                        ifilter(lambda x: x[1] >= self._minimum_support,
                                items_pair_frequency.iteritems()))):
                for element in (freq_1_seq_elementpool[item_i] +
                                freq_1_seq_elementpool[item_j]):
                    if element.support < self._minimum_support:
                        continue
                    freq_2_seq_elementpool[element.sequence] |= element
                    used_freq_items.update([item_i, item_j])

            if maximal:
                for item in used_freq_items:
                    del freq_1_seq_elementpool[item]

        return freq_1_seq_elementpool.values(), freq_2_seq_elementpool.values()

    def enumerate_frequent_sequences(self, elements, maximal=False,
                                     maximum_length=None, top_number=None,
                                     _enh=False):
        """
        Compute frequent k-sequences (k > 2) with Depth-First Search.

        @param elements: List of Element objects.
        @type elements: list
        @param maximal: Flag to get only maximal frequent sequences.
        @type maximal: bool
        @param maximum_length: The maximum length of sequential patterns.
        @type maximum_length: int/None
        @param top_number: The number of top longest output sequences.
        @type top_number: int/None
        @param _enh: Performance enhancement (_internal usage_).
        @type _enh: bool
        """
        elements.sort(key=lambda x: x.sequence)
        while len(elements):

            # - get elements that belong to one equivalence class -
            idx, equiv_cls_elements = 0, [elements.pop(0)]
            while idx < len(elements):
                if equiv_cls_elements[0].\
                        has_equivalence_relation(elements[idx]):
                    equiv_cls_elements.append(elements.pop(idx))
                else:
                    idx += 1

            if not maximal:
                self.add_elements(elements=equiv_cls_elements)

            # - get joined elements (k+1 sequences) -
            frequent_inner_elementpool = ElementPool()

            for idx_i in range(len(equiv_cls_elements)):
                for idx_j in range(idx_i, len(equiv_cls_elements)):
                    for element in (equiv_cls_elements[idx_i] +
                                    equiv_cls_elements[idx_j]):
                        if element.support < self._minimum_support:
                            continue
                        frequent_inner_elementpool[element.sequence] |= element

            if maximal and self._frequent_elementpool:
                for seq in frequent_inner_elementpool.keys():
                    if not self.is_maximal(frequent_inner_elementpool[seq]):
                        del frequent_inner_elementpool[seq]

            # - used for performance enhancement (further research is needed) -
            if _enh and frequent_inner_elementpool:
                frequent_inner_elements = frequent_inner_elementpool.values()
                idx_i = 0
                while idx_i < len(elements):
                    set_next_idx = True
                    idx_j = 0
                    while idx_j < len(frequent_inner_elements):
                        if elements[idx_i].\
                                has_subsequence(frequent_inner_elements[idx_j]):
                            del elements[idx_i]
                            frequent_inner_elements.append(
                                frequent_inner_elements.pop(idx_j))
                            set_next_idx = False
                            break
                        else:
                            idx_j += 1
                    if set_next_idx:
                        idx_i += 1
            # - performance enhancement - end -

            if (len(frequent_inner_elementpool) > 1
                    and (len(equiv_cls_elements[0]) + 1) != maximum_length):

                # - move to the next level (to get k+2 sequences) -
                self.enumerate_frequent_sequences(
                    elements=frequent_inner_elementpool.values(),
                    maximal=maximal,
                    maximum_length=maximum_length,
                    top_number=top_number,
                    _enh=_enh
                )

            else:

                # - frequent_inner_elementpool contains 0 or 1 element -
                for element in frequent_inner_elementpool.itervalues():

                    if maximal:
                        for seq in self._frequent_elementpool.keys():
                            if self._frequent_elementpool[seq].\
                                    has_subsequence(element):
                                del self._frequent_elementpool[seq]

                    self.add_elements(elements=[element], top_number=top_number)

    def prune(self, sequence):
        raise NotImplementedError

    def execute(self, maximal=False, sort=False,
                maximum_length=None, top_number=None):
        """
        Execute SPADE algorithm for defined data with certain minimum support.

        @param maximal: Flag to get only maximal frequent sequences.
        @type maximal: bool
        @param sort: Flag to sort the output base on sequence length.
        @type sort: bool
        @param maximum_length: The maximum length of sequential patterns.
        @type maximum_length: int/None
        @param top_number: The number of top longest output sequences.
        @type top_number: int/None
        @return: List of frequent sequences (elements of type Element).
        @rtype: list
        """
        self._frequent_elementpool.clear()

        freq_1s_elements, freq_2s_elements = self.generate_frequent_sequences(
            maximal=maximal,
            maximum_length=maximum_length
        )

        self.enumerate_frequent_sequences(
            elements=freq_2s_elements,
            maximal=maximal,
            maximum_length=maximum_length,
            top_number=top_number
        )

        self.add_elements(elements=(freq_1s_elements + freq_2s_elements),
                          top_number=top_number)

        frequent_elements = self._frequent_elementpool.values()
        if sort and not top_number:
            frequent_elements.sort(key=lambda x: len(x))

        return frequent_elements
