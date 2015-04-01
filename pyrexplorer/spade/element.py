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

__all__ = [
    'Element',
    'ElementPool'
]

from collections import namedtuple, defaultdict


class Element(object):

    """Class represents element from atom set with corresponding id-list."""

    Event = namedtuple('Event', ['sid', 'eid'])

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
            self.id_list.add(Element.Event(sid=sid, eid=eid))

        for event in events:
            self.id_list.add(event)

    def __len__(self):
        """
        Get lenght of Element's sequence (k: k-sequence).

        @return: Sum of itemset lengths (items per event).
        @rtype: int
        """
        return sum([len(x) for x in self.sequence])

    @property
    def support(self):
        """
        Get Element's sequence support (frequency of sequence).

        @return: Number of distinct sids.
        @rtype: int
        """
        return len(set([x.sid for x in self.id_list]))

    @staticmethod
    def get_difference_sequences(sequence_i, sequence_j):
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
            for idx in range(max(len(sequence_i), len(sequence_j))):

                try:
                    itemset_i = set(sequence_i[idx])
                except IndexError:
                    itemset_i = set()

                try:
                    itemset_j = set(sequence_j[idx])
                except IndexError:
                    itemset_j = set()

                diff_i.append(tuple(filter(
                    lambda x: x not in itemset_j, itemset_i)))
                diff_j.append(tuple(filter(
                    lambda x: x not in itemset_i, itemset_j)))

        return tuple(diff_i), tuple(diff_j)

    @staticmethod
    def event_atom_union(atom, item):
        """
        Get event atom (atom and item have the same eid).

        @param atom: Sequence of itemsets.
        @type atom: tuple of tuples
        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        _last_itemset = atom[-1]
        if item not in _last_itemset:
            _last_itemset += tuple([item])
        return atom[:-1] + tuple([tuple(sorted(_last_itemset))])

    @staticmethod
    def sequence_atom_union(atom, item):
        """
        Get sequence atom (atom and item have different eids).

        @param atom: Sequence of itemsets.
        @type atom: tuple of tuples
        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        return atom + tuple([tuple([item])])

    def get_event_atom_union(self, item):
        """
        Get event atom.

        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        return self.event_atom_union(atom=self.sequence, item=item)

    def get_sequence_atom_union(self, item):
        """
        Get sequence atom.

        @param item: Item object (str or int).
        @type item: type(Item)
        @return: United atom.
        @rtype: tuple of tuples
        """
        return self.sequence_atom_union(atom=self.sequence, item=item)

    def temporal_join(self, other):
        """
        Temporal join (of current element with other one).

        @param other: Element object.
        @type other: Element
        @return: New Element objects.
        @rtype: list
        """
        elements = ElementPool()

        seq_diff_i, seq_diff_j = self.get_difference_sequences(
            self.sequence, other.sequence)

        for pair_i in self.id_list:
            for pair_j in other.id_list:

                if pair_i.sid != pair_j.sid:
                    continue

                # - create sequence atom -
                if pair_i.eid < pair_j.eid:
                    itemset = seq_diff_j[-1] or other.sequence[-1]
                    atom = self.get_sequence_atom_union(itemset[-1])
                    eid = pair_j.eid
                # - create sequence atom -
                elif pair_i.eid > pair_j.eid:
                    itemset = seq_diff_i[-1] or self.sequence[-1]
                    atom = other.get_sequence_atom_union(itemset[-1])
                    eid = pair_i.eid
                # - create event atom -
                elif (seq_diff_i[-1] and seq_diff_j[-1]
                        and seq_diff_i[-1][0] != seq_diff_j[-1][0]):
                    atom = self.get_event_atom_union(seq_diff_j[-1][0])
                    eid = pair_i.eid
                else:
                    continue

                elements[atom] |= Element(atom, sid=pair_i.sid, eid=eid)

        return elements.values()

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
        return self.temporal_join(other)

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
