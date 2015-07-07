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

__all__ = ['Element', 'ElementDict', 'EVENT_ATOM_TYPE', 'SEQUENCE_ATOM_TYPE']

from collections import namedtuple

UNKNOWN_ATOM_TYPE = 0
EVENT_ATOM_TYPE = 1
SEQUENCE_ATOM_TYPE = 2  # EVENT_ATOM_TYPE < SEQUENCE_ATOM_TYPE

Event = namedtuple('Event', ['sid', 'eid'])


class Element(object):

    """Class represents atom set with corresponding id-list."""

    def __init__(self, item, prefix=None, conn_type=None, sid=None, eid=None):
        """
        Initialization.

        @param item: Atom's last item (i.e. key item).
        @type item: type(Item)
        @param prefix: Atom's prefix.
        @type prefix: tuple/None
        @param conn_type: Type of atom (how prefix connects to the last item).
        @type conn_type: int/None
        @param sid: Sequence id.
        @type sid: int/None
        @param eid: Event id.
        @type eid: int/None
        """
        self.key_item = item
        self.prefix = prefix
        self.conn_type = conn_type or UNKNOWN_ATOM_TYPE

        self.id_list = set()
        self.update_id_list(sid, eid)

    def update_id_list(self, sid=None, eid=None, id_list=None):
        """
        Update element's id-list with new event(s).

        @param sid: Sequence id.
        @type sid: int/None
        @param eid: Event id.
        @type eid: int/None
        @param id_list: List of Event objects (with parameters: sid and eid).
        @type id_list: list/None
        """
        if sid is not None and eid is not None:
            self.id_list.add(Event(sid=sid, eid=eid))
        self.id_list.update(id_list or [])

    @staticmethod
    def generate_sequence(last_item, prefix=None, is_sequence_atom=False):
        """
        Generate atom's sequence (prefix - conn_type - item).

        @param last_item: Atom's/sequence's last item (key item).
        @type last_item: type(Item)
        @param prefix: Atom's prefix.
        @type prefix: tuple/None
        @param is_sequence_atom: Flag that it is sequence atom (default=False).
        @type is_sequence_atom: bool
        @return: Sequence of itemsets (new atom).
        @rtype: tuple of tuples
        """
        prefix = prefix or ((),)
        if is_sequence_atom:
            output = prefix + tuple([tuple([last_item])])
        else:
            output = prefix[:-1] + tuple([prefix[-1] + tuple([last_item])])

        return output

    @property
    def sequence_length(self):
        """
        Get length of sequential pattern (k: k-sequence).

        @return: Sum of itemset lengths (items per event).
        @rtype: int
        """
        return sum([len(x) for x in self.prefix or ()]) + 1

    @property
    def sequence_size(self):
        """
        Get number of itemsets in Element's sequence.

        @return: Number of itemsets in sequence.
        @rtype: int
        """
        return (len(self.prefix or ((),)) +
                1 if self.conn_type == SEQUENCE_ATOM_TYPE else 0)

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

    @property
    def sequence(self):
        """
        Element's sequence.

        @return: Sequence (sequential pattern).
        @rtype: tuple of tuples
        """
        return self.generate_sequence(
            last_item=self.key_item,
            prefix=self.prefix,
            is_sequence_atom=(self.conn_type == SEQUENCE_ATOM_TYPE)
        )

    @classmethod
    def join(cls, element_i, element_j, cmap=None):
        """
        Temporal join of current element with other one (of the same prefix).

        @param element_i: Element object.
        @type element_i: Element
        @param element_j: Element object.
        @type element_j: Element
        @param cmap: Co-occurrence Map.
        @type cmap: dict/None
        @return: Set of Element objects grouped into ElementDict.
        @rtype: ElementDict
        """
        output = ElementDict()

        skip_by_cmap = False
        if cmap:

            if element_i.conn_type == element_j.conn_type == EVENT_ATOM_TYPE:

                if (element_i.key_item == element_j.key_item
                    or (element_i.key_item < element_j.key_item
                        and (element_j.key_item not in
                             cmap[element_i.key_item][EVENT_ATOM_TYPE]))
                    or (element_j.key_item < element_i.key_item
                        and (element_i.key_item not in
                             cmap[element_j.key_item][EVENT_ATOM_TYPE]))):
                    skip_by_cmap = True

            elif (element_i.conn_type == SEQUENCE_ATOM_TYPE
                    and element_j.conn_type == EVENT_ATOM_TYPE):

                if (element_i.key_item not in
                        cmap[element_j.key_item][SEQUENCE_ATOM_TYPE]):
                    skip_by_cmap = True

            elif (element_j.conn_type == SEQUENCE_ATOM_TYPE
                    and element_i.conn_type == EVENT_ATOM_TYPE):

                if (element_j.key_item not in
                        cmap[element_i.key_item][SEQUENCE_ATOM_TYPE]):
                    skip_by_cmap = True

        if not skip_by_cmap:

            for pair_i in element_i.id_list:
                for pair_j in element_j.id_list:

                    if pair_i.sid != pair_j.sid:
                        continue

                    # - create event atom -
                    if pair_i.eid == pair_j.eid:
                        if (element_i.conn_type != element_j.conn_type
                                or element_i.key_item == element_j.key_item):
                            continue
                        elif element_i.key_item < element_j.key_item:
                            key_item = element_j.key_item
                            prefix = element_i.sequence
                        elif element_i.key_item > element_j.key_item:
                            key_item = element_i.key_item
                            prefix = element_j.sequence

                        conn_type = EVENT_ATOM_TYPE
                        eid = pair_i.eid

                    # - create sequence atom -
                    else:
                        if pair_i.eid < pair_j.eid:
                            if element_j.conn_type == EVENT_ATOM_TYPE:
                                continue
                            key_item = element_j.key_item
                            prefix = element_i.sequence
                            eid = pair_j.eid

                        else:  # pair_i.eid > pair_j.eid
                            if element_i.conn_type == EVENT_ATOM_TYPE:
                                continue
                            key_item = element_i.key_item
                            prefix = element_j.sequence
                            eid = pair_i.eid

                        conn_type = SEQUENCE_ATOM_TYPE

                    sequence = cls.generate_sequence(
                        last_item=key_item,
                        prefix=prefix,
                        is_sequence_atom=(conn_type == SEQUENCE_ATOM_TYPE)
                    )

                    if sequence not in output:
                        output[sequence] = Element(
                            item=key_item, prefix=prefix, conn_type=conn_type)
                    output.add_event(key=sequence, sid=pair_i.sid, eid=eid)

        return output

    def __repr__(self):
        """
        String representation of an object.

        @return: Object description.
        @rtype: str
        """
        return '<Element: sequence=%s id_list=%s>' % (self.sequence,
                                                      self.id_list)


class ElementDict(object):

    """Class represents dictionary of elements."""

    def __init__(self):
        """Initialization."""
        self._data = {}

    def get(self, key):
        """
        Get element with defined key.

        @param key: Sequence or item.
        @type key: tuple/type(Item)
        @return: Element object.
        @rtype: Element/None
        """
        return self._data.get(key)

    def set(self, key, element):
        """
        Assign element to the defined key (sequence or item).

        @param key: Sequence or item.
        @type key: tuple/type(Item)
        @param element: Element object.
        @type element: Element
        """
        self._data[key] = element

    def remove(self, key):
        """
        Remove element with defined key from the dictionary.

        @param key: Sequence or item.
        @type key: tuple/type(Item)
        """
        if key in self:
            del self._data[key]

    def get_keys(self):
        """
        Get keys (sequences or items) of the dictionary.

        @return: List of keys.
        @rtype: list
        """
        return self._data.keys()

    def get_elements(self):
        """
        Get elements (without corresponding keys) from the dictionary.

        @return: List of Element objects.
        @rtype: list
        """
        return self._data.values()

    def update(self, key, element):
        """
        Update key's element with other element's id-list, if there is no key
        in dictionary then assign element to the defined key.

        @param key: Sequence or item.
        @type key: tuple/type(Item)
        @param element: Element object.
        @type element: Element
        """
        if key in self:
            self.get(key=key).update_id_list(id_list=element.id_list)
        else:
            self.set(key=key, element=element)

    def add_event(self, key, sid, eid):
        """
        Add event to the key's element.

        @param key: Sequence or item.
        @type key: tuple/type(Item)
        @param sid: Sequence id.
        @type sid: inte
        @param eid: Event id.
        @type eid: int
        """
        if key in self:
            self.get(key=key).update_id_list(sid=sid, eid=eid)

    def items(self):
        """
        Generator of (key, element) pairs.

        @return: Key-value pairs (key-element)
        @rtype: tuple
        """
        for key, element in self._data.iteritems():
            yield key, element

    def clear(self):
        """Clear dictionary data."""
        self._data.clear()

    def __getitem__(self, item):
        return self.get(key=item)

    def __setitem__(self, key, value):
        self.set(key=key, element=value)

    def __contains__(self, key):
        return key in self._data

    def __len__(self):
        return len(self._data)
