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

__all__ = ['SPADEm']

from collections import defaultdict, deque

from .element import Element, ElementPool


class SPADEm(object):

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
        @param top_number: The number of top longest output frequent sequences.
        @type top_number: int/None
        """
        for element in elements:
            self._frequent_elementpool[element.sequence] |= element
        self.keep_top_elements(top_number=top_number)

    def keep_top_elements(self, top_number):
        """
        Keep only elements with top longest sequential patterns.

        @param top_number: The number of top longest output frequent sequences.
        @type top_number: int
        """
        if top_number and len(self._frequent_elementpool) > top_number:

            top_frequent_sequences = set([
                element.sequence for element in sorted(
                    self._frequent_elementpool.values(),
                    key=lambda x: (len(x), x.num_itemsets),
                    reverse=True
                )[:top_number]
            ])

            for seq in self._frequent_elementpool.keys():
                if seq not in top_frequent_sequences:
                    del self._frequent_elementpool[seq]

    @staticmethod
    def sorted_elements(elements):
        """
        Get sorted Element objects.

        @param elements: List of Element objects.
        @type elements: list
        @return: Sorted (new) list of Element objects.
        @rtype: list
        """
        output = []

        elements_by_initial_event = {}
        for _ in xrange(len(elements)):
            element = elements.pop()
            min_eid = min(set([event.eid for event in element.id_list]))
            elements_by_initial_event.setdefault(min_eid, []).append(element)

        for eid in sorted(elements_by_initial_event.keys()):
            output.extend(sorted(elements_by_initial_event[eid],
                                 key=lambda x: (x.num_itemsets, x.sequence)))

        return output

    def generate_frequent_sequences(self, max_length=None):
        """
        Compute frequent 1-sequences and 2-sequences.

        @param max_length: The maximum length of sequential patterns.
        @type max_length: int/None
        @return: Two lists of frequent 1- and 2-sequences respectively.
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

        freq_1s_elementpool = ElementPool()

        for item in id_lists:
            if len(id_lists[item]) < self._minimum_support:
                continue
            for sid in id_lists[item]:
                for eid in id_lists[item][sid]:
                    freq_1s_elementpool[item] |= Element(item, sid=sid, eid=eid)

        itemspair_frequency = defaultdict(int)
        for items in sequences.itervalues():
            f_items = [x for x in items if x in freq_1s_elementpool]
            for idx_i in xrange(len(f_items)):
                for idx_j in xrange(idx_i + 1, len(f_items)):
                    itemspair_frequency[(f_items[idx_i], f_items[idx_j])] += 1

        freq_2s_elementpool = ElementPool()

        if max_length is None or max_length > 1:
            used_freq_items = set()

            for item_i, item_j in set([
                    tuple(sorted(k)) for k, v in itemspair_frequency.iteritems()
                    if v >= self._minimum_support]):

                joined_elements = freq_1s_elementpool[item_i].\
                    join(freq_1s_elementpool[item_j])

                if joined_elements is not None:
                    for element in joined_elements:
                        if element.support < self._minimum_support:
                            continue
                        freq_2s_elementpool[element.sequence] |= element
                        used_freq_items.update([item_i, item_j])

            for item in used_freq_items:
                del freq_1s_elementpool[item]

        return freq_1s_elementpool.values(), freq_2s_elementpool.values()

    def enumerate_frequent_sequences(self, elements,
                                     max_length=None, top_number=None):
        """
        Compute frequent k-sequences (k > 2) with Depth-First Search.

        @param elements: List of Element objects.
        @type elements: list
        @param max_length: The maximum length of sequential patterns.
        @type max_length: int/None
        @param top_number: The number of top longest output frequent sequences.
        @type top_number: int/None
        """
        elements_ = deque(self.sorted_elements(elements=elements))
        while len(elements_):

            frequent_elementpool_ = ElementPool()

            master_element = elements_.popleft()
            for _ in xrange(len(elements_)):
                joined_elements = master_element.join(elements_[0])
                if joined_elements is not None:  # same equivalence class
                    for element in joined_elements:
                        if element.support >= self._minimum_support:
                            frequent_elementpool_[element.sequence] |= element
                    elements_.popleft()
                else:
                    elements_.rotate(-1)

            if (len(frequent_elementpool_) > 1
                    and (len(master_element) + 1) != max_length):

                self.enumerate_frequent_sequences(
                    elements=frequent_elementpool_.values(),
                    max_length=max_length,
                    top_number=top_number
                )

            elif frequent_elementpool_:

                freq_elements = []
                for element in frequent_elementpool_.itervalues():

                    if not self.is_maximal(element):
                        continue

                    for seq in self._frequent_elementpool.keys():
                        if self._frequent_elementpool[seq].\
                                has_subsequence(element):
                            del self._frequent_elementpool[seq]

                    freq_elements.append(element)
                self.add_elements(elements=freq_elements, top_number=top_number)

            elif self.is_maximal(master_element):
                self.add_elements(elements=[master_element])

    def execute(self, sort=False, max_length=None, top_number=None):
        """
        Execute SPADE algorithm for defined data with certain minimum support.

        @param sort: Flag to sort the output base on sequence length.
        @type sort: bool
        @param max_length: The maximum length of sequential patterns.
        @type max_length: int/None
        @param top_number: The number of top longest output sequences.
        @type top_number: int/None
        @return: List of frequent sequences (elements of type Element).
        @rtype: list
        """
        self._frequent_elementpool.clear()

        freq_1s_elements, freq_2s_elements = \
            self.generate_frequent_sequences(max_length=max_length)

        if max_length is None or max_length > 2:
            self.enumerate_frequent_sequences(elements=freq_2s_elements,
                                              max_length=max_length,
                                              top_number=top_number)

        self.add_elements(elements=(freq_1s_elements + freq_2s_elements),
                          top_number=top_number)

        frequent_elements = self._frequent_elementpool.values()
        if sort:
            frequent_elements.\
                sort(key=lambda x: (len(x), x.num_itemsets, x.sequence))

        return frequent_elements
