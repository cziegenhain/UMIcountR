"""
This code is powered by the pybktree project. I chose to copy the class & functions here to avoid installation issues when loading the module instead

BK-tree data structure to allow fast querying of "close" matches.
pybktree code is licensed under a permissive MIT license -- see LICENSE.txt.
The pybktree project lives on GitHub here:
https://github.com/benhoyt/pybktree
"""

#import numpy as np
#import pybktree
#import pandas as pd
from collections import deque
from operator import itemgetter

_getitem0 = itemgetter(0)


def hamming_distance(x, y):
    """Calculate the hamming distance (number of bits different) between the
    two integers given.
    >>> [hamming_distance(x, 15) for x in [0, 8, 10, 12, 14, 15]]
    [4, 3, 2, 2, 1, 0]
    """
    return bin(x ^ y).count('1')


class BKTree(object):
    """BK-tree data structure that allows fast querying of matches that are
    "close" given a function to calculate a distance metric (e.g., Hamming
    distance or Levenshtein distance).
    Each node in the tree (including the root node) is a two-tuple of
    (item, children_dict), where children_dict is a dict whose keys are
    non-negative distances of the child to the current item and whose values
    are nodes.
    """
    def __init__(self, distance_func, items=[]):
        """Initialize a BKTree instance with given distance function
        (which takes two items as parameters and returns a non-negative
        distance integer). "items" is an optional list of items to add
        on initialization.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree.distance_func is hamming_distance
        True
        >>> tree = BKTree(hamming_distance, [])
        >>> list(tree)
        []
        >>> tree = BKTree(hamming_distance, [0, 4, 5])
        >>> sorted(tree)
        [0, 4, 5]
        """
        self.distance_func = distance_func
        self.tree = None

        _add = self.add
        for item in items:
            _add(item)

    def add(self, item):
        """Add given item to this tree.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree.add(4)
        >>> sorted(tree)
        [4]
        >>> tree.add(15)
        >>> sorted(tree)
        [4, 15]
        """
        node = self.tree
        if node is None:
            self.tree = (item, {})
            return

        # Slight speed optimization -- avoid lookups inside the loop
        _distance_func = self.distance_func

        while True:
            parent, children = node
            distance = _distance_func(item, parent)
            node = children.get(distance)
            if node is None:
                children[distance] = (item, {})
                break

    def find(self, item, n):
        """Find items in this tree whose distance is less than or equal to n
        from given item, and return list of (distance, item) tuples ordered by
        distance.
        >>> tree = BKTree(hamming_distance)
        >>> tree.find(13, 1)
        []
        >>> tree.add(0)
        >>> tree.find(1, 1)
        [(1, 0)]
        >>> for item in [0, 4, 5, 14, 15]:
        ...     tree.add(item)
        >>> sorted(tree)
        [0, 0, 4, 5, 14, 15]
        >>> sorted(tree.find(13, 1))
        [(1, 5), (1, 15)]
        >>> sorted(tree.find(13, 2))
        [(1, 5), (1, 15), (2, 4), (2, 14)]
        >>> sorted(tree.find(0, 1000)) == [(hamming_distance(x, 0), x) for x in tree]
        True
        """
        if self.tree is None:
            return []

        candidates = deque([self.tree])
        found = []

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend
        _found_append = found.append
        _distance_func = self.distance_func

        while candidates:
            candidate, children = _candidates_popleft()
            distance = _distance_func(candidate, item)
            if distance <= n:
                _found_append((distance, candidate))

            if children:
                lower = distance - n
                upper = distance + n
                _candidates_extend(c for d, c in children.items() if lower <= d <= upper)

        found.sort(key=_getitem0)
        return found

    def __iter__(self):
        """Return iterator over all items in this tree; items are yielded in
        arbitrary order.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree = BKTree(hamming_distance, [1, 2, 3, 4, 5])
        >>> sorted(tree)
        [1, 2, 3, 4, 5]
        """
        if self.tree is None:
            return

        candidates = deque([self.tree])

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend

        while candidates:
            candidate, children = _candidates_popleft()
            yield candidate
            _candidates_extend(children.values())

    def __repr__(self):
        """Return a string representation of this BK-tree with a little bit of info.
        >>> BKTree(hamming_distance)
        <BKTree using hamming_distance with no top-level nodes>
        >>> BKTree(hamming_distance, [0, 4, 8, 14, 15])
        <BKTree using hamming_distance with 3 top-level nodes>
        """
        return '<{} using {} with {} top-level nodes>'.format(
            self.__class__.__name__,
            self.distance_func.__name__,
            len(self.tree[1]) if self.tree is not None else 'no',
        )




def makeDNAbits(instring):
    lookupdict = {'A': '00001', 'C': '00010', 'G': '00100', 'T': '01000', 'N': '10000'}
    outbits = '1'
    for c in instring:
        outbits += lookupdict[c]
    return int(outbits,2)

def makeDNAstring(inbits):
    lookupdict = {'00001': 'A', '00010': 'C', '00100': 'G', '01000': 'T', '10000': 'N'}
    outstring = ''
    inbitsbin = '{0:08b}'.format(inbits)
    inbitsbinlist = [inbitsbin[i:i+5] for i in range(1,len(inbitsbin),5)]
    for c in inbitsbinlist:
        outstring += lookupdict[c]
    return(outstring)

def splitter(seq, breakpoint):
    return seq[:breakpoint], seq[breakpoint:]

def return_dist_frame(uniqumis, hamdist, bp):
    bp = int(bp)
    umibits = [makeDNAbits(u) for u in uniqumis]
    trees = [{},{}]
    for dna, bits in zip(uniqumis, umibits):
        part1, part2 = splitter(dna, bp)
        if not part1 in trees[0]:
            trees[0][part1] = BKTree(hamming_distance)
        trees[0][part1].add(bits)
        if not part2 in trees[1]:
            trees[1][part2] = BKTree(hamming_distance)
        trees[1][part2].add(bits)
    #print('trees created...')
    close_hits = [(trees[0][splitter(dna, bp)[0]].find(umi, hamdist*2), trees[1][splitter(dna, bp)[1]].find(umi, hamdist*2)) for dna, umi in zip(uniqumis, umibits)]
    #print('closest found')
    outlist = []
    for bs in close_hits:
        for b in bs:
            b_0 = b[0][1]
            b_rest = b[1:]
            outlist.extend([(makeDNAstring(b_0), b_1[0]/2, makeDNAstring(b_1[1])) for b_1 in b_rest])
    return(outlist)
