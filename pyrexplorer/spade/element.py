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

__all__ = ['Element', 'ElementPool']

from collections import namedtuple, defaultdict


Event = namedtuple('Event', ['sid', 'eid'])


class Element(object):

    """Class represents element from atom set with corresponding id-list."""

    def __init__(self, sequence, sid=None, eid=None, *events):
        """
        Initialization.

        @param sequence: Sequence of itemsets (or single Item object).
        @type sequence: tuple of tuples/type(Item)
        @param sid: Sequence id.
        @type sid: int/None
        @param eid: Event id.
        @type eid: int/None
        @param events: List of Event objects (with parameters: sid and eid).
        @type events: list
        """
        if not isinstance(sequence, (list, tuple)):
            sequence = tuple([tuple([sequence])])

        self.sequence = tuple(sequence)
        self.id_list = set()

        if sid and eid:
            self.id_list.add(Event(sid=sid, eid=eid))

        for event in events:
            self.id_list.add(event)

    def __len__(self):
        """
        Get length of sequential pattern (k: k-sequence).

        @return: Sum of itemset lengths (items per event).
        @rtype: int
        """
        return sum([len(x) for x in self.sequence])

    @property
    def num_itemsets(self):
        """
        Get number of itemsets in Element's sequence.

        @return: Number of itemsets in sequence.
        @rtype: int
        """
        return len(self.sequence)

    @property
    def support(self):
        """
        Get Element's sequence support (frequency of sequence).

        @return: Number of distinct sids.
        @rtype: int
        """
        sids = set()
        for event in self.id_list:
            sids.add(event.sid)
        return len(sids)

    @staticmethod
    def get_sequences_diff(sequence_i, sequence_j):
        """
        Get difference between two sequences.

        @param sequence_i: First sequence.
        @type sequence_i: tuple of tuples
        @param sequence_j: Second sequence.
        @type sequence_j: tuple of tuples
        @return: 2 sequences of diffs between each element of input sequences.
        @rtype: tuple of 2 sequences (tuple of tuples)
        """
        diff_i, diff_j = [], []

        if sequence_i and sequence_j:
            for idx in xrange(max(len(sequence_i), len(sequence_j))):

                try:
                    iset_i = set(sequence_i[idx])
                except IndexError:
                    iset_i = set()

                try:
                    iset_j = set(sequence_j[idx])
                except IndexError:
                    iset_j = set()

                diff_i.append(tuple([x for x in iset_i if x not in iset_j]))
                diff_j.append(tuple([x for x in iset_j if x not in iset_i]))

        return tuple(diff_i), tuple(diff_j)

    def get_event_atom_union(self, item):
        """
        Get event atom (self sequence and item have the same eid).

        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        return (self.sequence[:-1] +
                tuple([tuple(sorted(set(self.sequence[-1] + tuple([item]))))]))

    def get_sequence_atom_union(self, item):
        """
        Get sequence atom (self sequence and item have different eids).

        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        return self.sequence + tuple([tuple([item])])

    def get_equivalence_relation_diff(self, other):
        """
        Get difference between two last itemsets of elements' sequences.

        @param other: Element object.
        @type other: Element
        @return: Difference between last itemsets in elements' sequences.
        @rtype: tuple(<item>/None, <item>/None)
        """
        output = [None, None]

        max_eq_idx = max(len(self.sequence), len(other.sequence)) - 2
        if self.sequence[:max_eq_idx] == other.sequence[:max_eq_idx]:

            seq_diff_i, seq_diff_j = self.get_sequences_diff(
                self.sequence[max_eq_idx:], other.sequence[max_eq_idx:])

            if len(seq_diff_i[-1]) == len(seq_diff_j[-1]) == 1:
                if (not [x for x in seq_diff_i[:-1] if x]
                        and not [x for x in seq_diff_j[:-1] if x]):
                    output[0] = seq_diff_i[-1][0]
                    output[1] = seq_diff_j[-1][0]

            elif len(seq_diff_i) > 1 and len(seq_diff_j) > 1:

                if (len(seq_diff_i[-1]) == len(seq_diff_j[-2]) == 1
                        and not seq_diff_i[-2] and not seq_diff_j[-1]
                        and (len(self.sequence) - len(other.sequence)) == 1):
                    output[0] = seq_diff_i[-1][0]

                elif (len(seq_diff_j[-1]) == len(seq_diff_i[-2]) == 1
                        and not seq_diff_j[-2] and not seq_diff_i[-1]
                        and (len(other.sequence) - len(self.sequence)) == 1):
                    output[1] = seq_diff_j[-1][0]

        return tuple(output)

    def join(self, other):
        """
        Temporal join of current element with other one (of the same prefix).

        @param other: Element object.
        @type other: Element
        @return: New Element objects.
        @rtype: list/None
        """
        last_item_i, last_item_j = self.get_equivalence_relation_diff(other)
        if last_item_i is None and last_item_j is None:
            return None

        elementpool = ElementPool()

        for pair_i in self.id_list:
            for pair_j in other.id_list:

                if pair_i.sid != pair_j.sid:
                    continue

                # - create sequence atom -
                if pair_i.eid < pair_j.eid and last_item_j:
                    atom = self.get_sequence_atom_union(last_item_j)
                    eid = pair_j.eid
                # - create sequence atom -
                elif pair_i.eid > pair_j.eid and last_item_i:
                    atom = other.get_sequence_atom_union(last_item_i)
                    eid = pair_i.eid
                # - create event atom -
                elif (pair_i.eid == pair_j.eid and last_item_i and last_item_j
                        and last_item_i != last_item_j):
                    atom = self.get_event_atom_union(last_item_j)
                    eid = pair_i.eid
                else:
                    continue

                elementpool[atom] |= Element(atom, sid=pair_i.sid, eid=eid)

        return elementpool.values()

    def has_subsequence(self, other):
        """
        Check if self's sequence is sub-sequence for other's sequence.

        @param other: Element object.
        @type other: Element
        @return: Flag that self has sub-sequence for other's sequence.
        @rtype: bool
        """
        counter, idx_start_j = 0, 0
        for idx_i in xrange(len(self.sequence)):
            for idx_j in xrange(idx_start_j, len(other.sequence)):

                master_itemset = set(other.sequence[idx_j])

                is_subitemset = True
                for x in self.sequence[idx_i]:
                    if x not in master_itemset:
                        is_subitemset = False
                        break

                if not is_subitemset:
                    continue

                counter += 1

                idx_start_j = idx_j + 1
                break

            if (idx_i + 1) != counter:
                break

        return counter == len(self.sequence)

    def __eq__(self, other):
        """
        Implements the comparison operator "==".

        @param other: Element object.
        @type other: Element
        @return: Flag that says two objects are equal or not.
        @rtype: bool
        """
        return self.sequence == other.sequence and self.id_list == other.id_list

    def __ior__(self, other):
        """
        Implements the assignment operator "|=".

        @param other: Element object.
        @type other: Element
        @return: Element with events as a union of the events of both Elements.
        @rtype: self
        """
        if self.sequence == other.sequence:
            self.id_list |= other.id_list
        return self

    def __add__(self, other):
        """
        Implements the assignment operator "+".

        @param other: Element object.
        @type other: Element
        @return: New Element objects.
        @rtype: list
        """
        return self.join(other)

    def __repr__(self):
        """
        String representation of an object.

        @return: Object description.
        @rtype: str
        """
        return self.__dict__.__repr__()


class ElementPool(defaultdict):

    """Dictionary class with pre-defined value type."""

    def __missing__(self, key):
        self[key] = value = Element(key)  # key can be either <item> or <seq>
        return value
