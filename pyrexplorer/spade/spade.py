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
        self._id_lists = {}
        self._minimum_support = None

    def set(self, **kwargs):
        """
        Set initial data.

        @param kwargs: Input parameters.
        @type kwargs: dict

        @keyword sequences: Dictionary of sequences {sid: {eid: <itemset>}}.
        @keyword minimum_support: Minimum support (number of distinct sids).
        """
        if kwargs.get('sequences'):
            self._id_lists.clear()
            for sid in kwargs['sequences']:
                for eid in kwargs['sequences'][sid]:
                    for item in set(kwargs['sequences'][sid][eid]):
                        self._id_lists.\
                            setdefault(item, {}).\
                            setdefault(sid, []).\
                            append(eid)
        self._minimum_support = kwargs.get('minimum_support')

    def generate_frequent_sequences(self):
        """
        Get initial data and compute frequent 1-sequences and 2-sequences.

        @return: Two lists of frequent 1 and 2-sequences.
        @rtype: tuple(list, list)
        """
        if not self._id_lists or not self._minimum_support:
            raise Exception('Initial sequences/support are not defined')

        sequences = {}
        freq_1_seq_elements = ElementPool()
        for item in sorted(self._id_lists):

            if self._minimum_support > len(self._id_lists[item]):
                continue

            for sid in self._id_lists[item]:
                sequences.setdefault(sid, [])
                for eid in sorted(self._id_lists[item][sid]):
                    freq_1_seq_elements[item] |= Element(item, sid=sid, eid=eid)
                    sequences[sid].append((item, eid))

        items_pair_frequency = defaultdict(int)
        for pairs_list in sequences.itervalues():
            for idx_i in range(len(pairs_list)):
                item_i, eid_i = pairs_list[idx_i]
                for idx_j in range(idx_i + 1, len(pairs_list)):
                    item_j, eid_j = pairs_list[idx_j]
                    if eid_i > eid_j:
                        items_pair = (item_j, item_i)
                    else:
                        items_pair = (item_i, item_j)
                    items_pair_frequency[items_pair] += 1

        freq_2_seq_elements = ElementPool()
        for item_i, item_j in set(map(
                lambda x: tuple(sorted(x[0])),
                filter(lambda x: x[1] >= self._minimum_support,
                       items_pair_frequency.iteritems()))):

            for element in (freq_1_seq_elements[item_i] +
                            freq_1_seq_elements[item_j]):
                if element.support >= self._minimum_support:
                    freq_2_seq_elements[element.sequence] |= element

        return freq_1_seq_elements.values(), freq_2_seq_elements.values()

    def enumerate_frequent_sequences(self, elements, k=None):
        """
        Compute frequent k-sequences (k > 2).

        @param elements: List of Element objects.
        @type elements: list
        @param k: The lenght of output sequences.
        @type k: int/None
        @return: List of frequent sequences (elements of type Element).
        @rtype: list
        """
        freq_seq_elements = ElementPool()
        for idx_i in range(len(elements)):

            if k and len(elements[idx_i]) == k:
                break

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

        return freq_seq_elements.values()

    def prune(self, sequence):
        raise NotImplementedError

    def execute(self, sort=False, k=None):
        """
        Execute SPADE algorithm for defined data with certain minimum support.

        @param sort: Flag to sort the output base on sequence length.
        @type sort: bool
        @param k: The lenght of output sequences.
        @type k: int/None
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
            freq_seq = self.enumerate_frequent_sequences(freq_2_seq, k=k)
            if k:
                freq_seq = filter(lambda x: len(x) == k, freq_seq)
            elif sort:
                freq_seq = sorted(freq_seq, key=lambda x: len(x))
            output.extend(freq_seq)

        return output
