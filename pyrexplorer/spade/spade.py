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

from .element import Element, ElementDict, EVENT_ATOM_TYPE, SEQUENCE_ATOM_TYPE


def is_subsequence(sequence_i, sequence_j, level=0):
    """
    Check if sequence_i is sub-sequence for sequence_j.

    @param sequence_i: First sequence (i.e. possible sub-sequence).
    @type sequence_i: tuple
    @param sequence_j: Second sequence (i.e. master sequence)
    @type sequence_j: tuple
    @param level: Number of sub-tuples (possible values: 0, 1).
    @type level: int
    @return: Flag that the 1st sequence is sub-sequence for the 2nd one.
    @rtype: bool
    """
    counter, idx_start_j = 0, 0
    for idx_i in xrange(len(sequence_i)):
        for idx_j in xrange(idx_start_j, len(sequence_j)):

            if not level:
                item_j = sequence_j[idx_j] if idx_i else abs(sequence_j[idx_j])
                if sequence_i[idx_i] != item_j:
                    continue

            elif level == 1:
                master_itemset = set(sequence_j[idx_j])

                is_subitemset = True
                for x in sequence_i[idx_i]:
                    if x not in master_itemset:
                        is_subitemset = False
                        break

                if not is_subitemset:
                    continue

            else:
                break

            counter += 1

            idx_start_j = idx_j + 1
            break

        if (idx_i + 1) != counter:
            break

    return counter == len(sequence_i)


class SPADEm(object):

    def __init__(self):
        """Initialization."""
        self._sequences = {}
        self._minimum_support = None

        self._cmap = {}
        self._frequent_elementdict = ElementDict()

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

    def is_maximal_sequence(self, element_sequence):
        """
        Check that element's sequence is maximal frequent sequence.

        @param element_sequence: Element's sequence.
        @type element_sequence: tuple of tuples
        @return: Flag that element's seq is not sub-seq for any frequent seq.
        @rtype: bool
        """
        output = True

        for sequence, _ in self._frequent_elementdict.items():
            if is_subsequence(element_sequence, sequence, level=1):
                output = False
                break

        return output

    def add_elements(self, elementdict, top_number=None):
        """
        Add elements to the ElementDict of frequent sequences.

        @param elementdict: Set of Element objects grouped into ElementDict.
        @type elementdict: ElementDict
        @param top_number: The number of top longest output frequent sequences.
        @type top_number: int/None
        """
        for sequence, element in elementdict.items():
            self._frequent_elementdict.update(key=sequence, element=element)

        if top_number and len(self._frequent_elementdict) > top_number:

            sequences = []
            for sequence, element in self._frequent_elementdict.items():
                sequences.append((element.sequence_length,
                                  element.sequence_size,
                                  sequence))
            sequences.sort(reverse=True)
            for sequence in set([x[2] for x in sequences[top_number:]]):
                self._frequent_elementdict.remove(key=sequence)

    def grouped(self, elements):
        """
        Get filtered and grouped Element objects (and sorted groups).

        @param elements: List of Element objects.
        @type elements: list
        @return: Element objects grouped by equivalence class.
        @rtype: collections.deque
        """
        output = []

        prefixes, grouped_elements = {}, {}
        for _ in xrange(len(elements)):
            element = elements.pop()

            min_eids = {}
            for event in element.id_list:
                min_eid = min_eids.setdefault(event.sid, event.eid)
                if event.eid < min_eid:
                    min_eids[event.sid] = event.eid

            prefix_min_eids = prefixes.setdefault(element.prefix, {})
            for sid, eid in min_eids.iteritems():
                min_eid = prefix_min_eids.setdefault(sid, eid)
                if eid < min_eid:
                    prefix_min_eids[sid] = eid

            grouped_elements.setdefault(element.prefix, []).append(
                (element, tuple(sorted(min_eids.items()))))

        grouped_maximum_sequences = deque([])
        for prefix, prefix_eids in sorted(prefixes.items(),
                                          key=lambda e: (sum(e[1].values()),
                                                         min(e[1].values()),
                                                         len(e[0]),
                                                         e[0])):
            # filter prefixes
            if len(prefixes) > 1:

                sid_parameters = {}
                for element, eids in grouped_elements[prefix]:
                    for event in element.id_list:
                        sid_parameters.setdefault(event.sid, []).append(
                            (event.eid, element.conn_type, element.key_item))

                sequences_per_sid = {}
                for sid in sid_parameters.keys():
                    sequences_per_sid[sid] = tuple(
                        [prefix[-1][-1]] +
                        [(-x[2] if x[1] == SEQUENCE_ATOM_TYPE else x[2])
                         for x in sorted(sid_parameters[sid])]
                    )
                    del sid_parameters[sid]

                for _ in xrange(len(grouped_maximum_sequences)):
                    grouped_sequences = grouped_maximum_sequences[0]

                    skip_sequences_set = False
                    for sid in grouped_sequences:
                        if not is_subsequence(grouped_sequences[sid],
                                              sequences_per_sid.get(sid, ())):
                            skip_sequences_set = True
                            break

                    if skip_sequences_set:
                        grouped_maximum_sequences.rotate(-1)
                    else:
                        grouped_maximum_sequences.popleft()

                skip_prefix = False
                for grouped_sequences in grouped_maximum_sequences:

                    skip_sequences_set = False
                    for sid in sequences_per_sid:
                        if not is_subsequence(sequences_per_sid[sid],
                                              grouped_sequences.get(sid, ())):
                            skip_sequences_set = True
                            break

                    if not skip_sequences_set:
                        skip_prefix = True
                        break

                if skip_prefix:
                    continue

                grouped_maximum_sequences.append(sequences_per_sid)

            if len(grouped_elements[prefix]) == 1:
                output.append({
                    'idx': [0],
                    'elements': deque([grouped_elements[prefix][0][0]])
                })
                continue

            grouped_elements[prefix].sort(key=lambda e: (e[0].conn_type,
                                                         e[1],
                                                         e[0].key_item))

            if grouped_elements[prefix][0][0].conn_type == EVENT_ATOM_TYPE:

                idx_list = []
                prev_min_eids = None
                for idx, p in enumerate(grouped_elements[prefix]):

                    if p[0].conn_type != EVENT_ATOM_TYPE:
                        break

                    elif p[1] != prev_min_eids:
                        idx_list.append(idx)
                        prev_min_eids = p[1]

                output.append({
                    'idx': idx_list,
                    'elements': deque([x[0] for x in grouped_elements[prefix]])
                })

            else:

                elements_ = [x[0] for x in sorted(grouped_elements[prefix],
                                                  key=lambda e: e[0].key_item)]

                idx_dict = {}
                for idx_i in xrange(len(elements_)):
                    item = elements_[idx_i].key_item

                    current_group = set()
                    for idx_j in xrange(idx_i + 1, len(elements_)):
                        if (elements_[idx_j].key_item
                                in self._cmap[item][EVENT_ATOM_TYPE]):
                            current_group.add(elements_[idx_j].key_item)
                    current_group.add(item)

                    skip_group = False
                    for idx in idx_dict:
                        counter = 0
                        for gr_item in current_group:
                            if gr_item in idx_dict[idx]:
                                counter += 1
                        if counter == len(current_group):
                            skip_group = True
                            break

                    if not skip_group:
                        idx_dict[idx_i] = current_group

                output.append({
                    'idx': idx_dict.keys(),
                    'elements': deque(elements_)
                })

        return deque(output)

    def update_cmap(self, element):
        """
        Update cmap (Co-occurrence Map) with data from 2-sequence element.

        @param element: Element object.
        @type element: Element
        """
        master_item, connected_item = element.prefix[0][0], element.key_item
        for item in [master_item, connected_item]:
            self._cmap.setdefault(item, {EVENT_ATOM_TYPE: set(),
                                         SEQUENCE_ATOM_TYPE: set()})
        self._cmap[master_item][element.conn_type].add(connected_item)

    def generate_frequent_sequences(self, max_length=None):
        """
        Compute frequent 1-sequences and 2-sequences.

        @param max_length: The maximum length of sequential patterns.
        @type max_length: int/None
        @return: Two ElementDicts of frequent 1- and 2-sequences respectively.
        @rtype: tuple(ElementDict, ElementDict)
        """
        if not self._sequences or not self._minimum_support:
            raise Exception('Initial sequences/support are not set')

        id_lists, sequences = {}, {}
        for sid in self._sequences:
            sequences.setdefault(sid, [])
            for e_idx, eid in enumerate(sorted(self._sequences[sid])):
                for item in sorted(set(self._sequences[sid][eid])):
                    # use index of itemset ("e_idx") instead of actual "eid"
                    id_lists.\
                        setdefault(item, {}).\
                        setdefault(sid, []).\
                        append(e_idx)
                    sequences[sid].append(item)

        freq_1s_elementdict = ElementDict()

        for item in id_lists:
            if len(id_lists[item]) < self._minimum_support:
                continue
            freq_1s_elementdict[item] = Element(item=item)
            for sid in id_lists[item]:
                for eid in id_lists[item][sid]:
                    freq_1s_elementdict.add_event(key=item, sid=sid, eid=eid)

        itemspair_frequency = defaultdict(int)
        for items in sequences.itervalues():
            f_items = [x for x in items if x in freq_1s_elementdict]
            for idx_i in xrange(len(f_items)):
                for idx_j in xrange(idx_i + 1, len(f_items)):
                    itemspair_frequency[(f_items[idx_i], f_items[idx_j])] += 1

        freq_2s_elementdict = ElementDict()

        if max_length is None or max_length > 1:
            used_freq_items = set()

            for item_i, item_j in set([
                    tuple(sorted(k)) for k, v in itemspair_frequency.iteritems()
                    if v >= self._minimum_support]):

                elementdict_ = Element.join(
                    element_i=freq_1s_elementdict[item_i],
                    element_j=freq_1s_elementdict[item_j]
                )

                counter = 0
                for sequence, element in elementdict_.items():
                    if element.support < self._minimum_support:
                        continue

                    if sequence not in freq_2s_elementdict:
                        self.update_cmap(element=element)

                    freq_2s_elementdict.update(key=sequence, element=element)
                    counter += 1

                if counter:
                    used_freq_items.update([item_i, item_j])

            for item in used_freq_items:
                freq_1s_elementdict.remove(key=item)

        return freq_1s_elementdict, freq_2s_elementdict

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
        frequent_inner_elementdict = ElementDict()
        frequent_master_elementdict = ElementDict()

        grouped_elements = self.grouped(elements=elements)
        while grouped_elements:

            data = grouped_elements.popleft()
            current_element_length = data['elements'][0].sequence_length

            for master_idx in data['idx']:
                master_element = data['elements'][master_idx]

                counter = 0
                for idx in xrange(len(data['elements'])):
                    current_element = data['elements'][idx]

                    if current_element is None:
                        continue

                    elementdict_ = Element.join(element_i=master_element,
                                                element_j=current_element,
                                                cmap=self._cmap)

                    for sequence, element in elementdict_.items():
                        if element.support < self._minimum_support:
                            continue

                        frequent_inner_elementdict.update(key=sequence,
                                                          element=element)
                        counter += 1

                if not counter:
                    sequence = master_element.sequence
                    if self.is_maximal_sequence(element_sequence=sequence):
                        frequent_master_elementdict[sequence] = master_element

                data['elements'][master_idx] = None

            if (len(frequent_inner_elementdict) > 1
                    and (current_element_length + 1) != max_length):

                new_grouped_elements = self.grouped(
                    elements=frequent_inner_elementdict.get_elements())
                new_grouped_elements.extend(grouped_elements)
                grouped_elements = new_grouped_elements
                # new_grouped_elements.reverse()  # python >= 2.7 (!)
                # grouped_elements.extendleft(new_grouped_elements)

            elif len(frequent_inner_elementdict):

                for sequence in frequent_inner_elementdict.get_keys():

                    if not self.is_maximal_sequence(element_sequence=sequence):
                        frequent_inner_elementdict.remove(key=sequence)
                        continue

                    for freq_sequence in self._frequent_elementdict.get_keys():
                        if is_subsequence(freq_sequence, sequence, level=1):
                            self._frequent_elementdict.remove(key=freq_sequence)

                self.add_elements(elementdict=frequent_inner_elementdict,
                                  top_number=top_number)

            if len(frequent_master_elementdict):
                self.add_elements(elementdict=frequent_master_elementdict)

            frequent_inner_elementdict.clear()
            frequent_master_elementdict.clear()

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
        self._frequent_elementdict.clear()
        self._cmap.clear()

        freq_1s_elementdict, freq_2s_elementdict = \
            self.generate_frequent_sequences(max_length=max_length)

        if len(freq_2s_elementdict) and (max_length is None or max_length > 2):
            self.enumerate_frequent_sequences(
                elements=freq_2s_elementdict.get_elements(),
                max_length=max_length,
                top_number=top_number
            )

        if len(freq_1s_elementdict):
            self.add_elements(elementdict=freq_1s_elementdict,
                              top_number=top_number)

        frequent_elements = self._frequent_elementdict.get_elements()
        if sort:
            frequent_elements.sort(key=lambda x: (x.sequence_length,
                                                  x.sequence_size,
                                                  x.prefix,
                                                  x.key_item))

        return frequent_elements
