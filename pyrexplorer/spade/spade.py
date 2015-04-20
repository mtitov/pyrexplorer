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
# - Shahin Saneinejad, <shahin.saneinejad@gmail.com, github.com/shahin>, 2012
# - Mikhail Titov, <mikhail.titov@cern.ch>, 2015
#

__all__ = ['SPADE']

from collections import defaultdict

from .element import Element, ElementPool


class SPADE(object):

    def __init__(self):
        """Initialization."""
        self._sequences = {}
        self._minimum_support = None

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

    def generate_frequent_sequences(self):
        """
        Get initial data and compute frequent 1-sequences and 2-sequences.

        @return: Two lists of frequent 1- and 2-sequences.
        @rtype: tuple(list, list)
        """
        if not self._sequences or not self._minimum_support:
            raise Exception('Initial sequences/support are not set')

        id_lists, sequences = {}, {}
        for sid in self._sequences:
            sequences.setdefault(sid, [])
            for eid in sorted(self._sequences[sid]):
                for item in sorted(self._sequences[sid][eid]):
                    id_lists.\
                        setdefault(item, {}).\
                        setdefault(sid, []).\
                        append(eid)
                    sequences[sid].append(item)

        freq_1_seq_elements = ElementPool()
        for item in id_lists:
            if self._minimum_support > len(id_lists[item]):
                continue
            for sid in id_lists[item]:
                for eid in id_lists[item][sid]:
                    freq_1_seq_elements[item] |= Element(item, sid=sid, eid=eid)

        items_pair_frequency = defaultdict(int)
        for items in sequences.itervalues():
            items = filter(lambda x: x in freq_1_seq_elements, items)
            for idx_i in range(len(items)):
                for idx_j in range(idx_i + 1, len(items)):
                    items_pair_frequency[(items[idx_i], items[idx_j])] += 1

        freq_2_seq_elements = ElementPool()
        for item_i, item_j in set(
                map(lambda x: tuple(sorted(x[0])),
                    filter(lambda x: x[1] >= self._minimum_support,
                           items_pair_frequency.iteritems()))):
            for element in (freq_1_seq_elements[item_i] +
                            freq_1_seq_elements[item_j]):
                if element.support >= self._minimum_support:
                    freq_2_seq_elements[element.sequence] |= element

        return freq_1_seq_elements.values(), freq_2_seq_elements.values()

    def enumerate_frequent_sequences(self, elements, k=None, top_number=None):
        """
        Compute frequent k-sequences (k > 2) with Depth-First Search.

        @param elements: List of Element objects.
        @type elements: list
        @param k: The maximum lenght of output sequences.
        @type k: int/None
        @param top_number: The number of top longest sequential patterns.
        @type top_number: int/None
        @return: List of frequent sequences (elements of type Element).
        @rtype: list
        """
        if k and len(elements) and len(elements[0]) == k:
            return []

        freq_seq_elements = ElementPool()
        for idx_i in range(len(elements)):

            freq_seq_inner_elements = ElementPool()
            for idx_j in range(idx_i + 1, len(elements)):

                for element in (elements[idx_i] + elements[idx_j]):
                    if element.support >= self._minimum_support:
                        freq_seq_inner_elements[element.sequence] |= element

            for element in freq_seq_inner_elements.itervalues():
                freq_seq_elements[element.sequence] |= element

            for element in self.enumerate_frequent_sequences(
                    freq_seq_inner_elements.values(), k=k):
                freq_seq_elements[element.sequence] |= element

        return self.get_top_elements(freq_seq_elements.values(), top_number)

    def prune(self, sequence):
        raise NotImplementedError

    @staticmethod
    def get_top_elements(elements, top_number=None):
        """
        Get top longest sequential patterns.

        @param elements: List of sequences (elements of type Element).
        @type: list
        @param top_number: The number of top longest sequential patterns.
        @type top_number: int/None
        @return: Top-n longest sequential patterns.
        @rtype: list
        """
        if not top_number:
            return elements
        return sorted(elements,
                      key=lambda x: (x.sequence_length, len(x)),
                      reverse=True)[:top_number]

    def execute(self, sort=False, k=None, top_number=None):
        """
        Execute SPADE algorithm for defined data with certain minimum support.

        @param sort: Flag to sort the output base on sequence length.
        @type sort: bool
        @param k: The lenght of output sequences.
        @type k: int/None
        @param top_number: The number of top longest sequential patterns.
        @type top_number: int/None
        @return: List of frequent sequences (elements of type Element).
        @rtype: list
        """
        output = []

        freq_1_seq, freq_2_seq = self.generate_frequent_sequences()

        if not k or k == 1:
            output.extend(freq_1_seq)
        if not k or k == 2:
            output.extend(freq_2_seq)
        if not k or k > 2:
            freq_seq = self.enumerate_frequent_sequences(freq_2_seq,
                                                         k=k,
                                                         top_number=top_number)
            if not top_number:
                if k:
                    freq_seq = filter(lambda x: len(x) == k, freq_seq)
                elif sort:
                    freq_seq = sorted(freq_seq, key=lambda x: len(x))
            output.extend(freq_seq)

        return self.get_top_elements(elements=output, top_number=top_number)
