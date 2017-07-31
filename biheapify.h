/*
 * biheapify.h
 *
 *  Created on: Jun 12, 2017
 *      Author: Matthew Gregory Krupa
 */
/* BIHEAP DEFINITION (QUICK):
We're given total_num_nodes elements indexed by 0, ..., total_num_nodes - 1.
These elements form a BiHeap if:
(1) the elements [0, ..., num_nodes_in_heap) form a min heap with the element
    at 0 being the min, and
(2) the elements [total_num_nodes - num_nodes_in_heap, ..., total_num_nodes)
   form a max heap with the element at total_num_nodes - 1 being the max,
 where num_nodes_in_heap is obtained from the function
GetNumNodesInHeapContainedInBiheap(), which is defined in biheap_common.h.
The quantity num_nodes_in_heap is a fundamentally important quantity associated
 with the biheap on total_num_nodes nodes. The formulas that define this
 quantity, although complicated, stem from the intuitive and natural
 graph-theoretic definition of a biheap given in the detailed definition below.
Given an iterator first, the function IsBiheap(), defined in biheap_common.h,
 may be used to check that the first total_num_nodes elements iterated by first
 (i.e. *first, *(first + 1), ..., *(first + total_num_nodes - 1))
 form a biheap.
It would likely be instructive to read the definition of
 IsBiheap(RAI first, size_type total_num_nodes) found in biheap_common.h
 before reading the more detailed definition of a biheap given below.
 */
/*
BIHEAP DEFINITION (DETAILED):
We define a new type of data structure, called a biheap, that
(a) is not excessively computationally costly to create, especially when
    there are an even number of elements being considered, and
(b) is not an excessively complicated operation, being not much more difficult
    to implement than a min heap and a max heap; indeed the BiHeapify()
    operation is primarily a combination of the four standard sifting
    operations, which are sifting elements up/down in a min/max heap operations.
(c) makes the minimum and maximum elements readily accessible.
(d) has a part of its data in a min heap, a part in a max heap, and a
    part that is in both simultaneously.
These properties make the biheap definition and the biheapify operation
easier to implement than double-ended heaps, its nearest counterpart, while
also having good performance. The author considers the definition of a biheap
to be a very natural data structure to define. It is also a definition that
appears to have gone unnoticed until now.

Besides making the minimum and maximum elements readily available in O(n) time
(based on experimentation), the author has noticed that, as expected, the
BiHeapify() functions tends to take randomly distributed elements (following
a uniform distribution) and reorders them so that many of the larger elements
are near the max and many of the smaller elements are near the min and the
value in middle of the data to be close to the median. Therefore, it is
reasonable to expect that applying BiHeapify to a collection of data before
sorting it may help to speed up certain sorting algorithms such as quicksort.
About the name: Various types of data structures called double ended heaps
(also known as double ended queues) have already been discovered, but it is
not clear that a biheap can be made into a double ended heap since it
is not clear at the moment whether or not O(log n) pop and push operations
 exist for biheaps. In this sense, the biheap data structure defined here may
not be a specific type of double ended heap; it is nevertheless called a biheap
due to it, by definition, simultaneously being both a min heap and a max heap.

Suppose that we have a positive even number total_num_nodes of nodes
 vec[0 : total_num_nodes) (i.e. nodes vec[0], ..., vec[total_num_nodes - 1])
 and call each vec[i] a "node".
1) Let F_n (L_n) denote the First (resp. Last) n nodes of vec[0 : total_num_nodes)
2) Define a binary tree N_n on F_n rooted at vec[0] in the usual way and let
   PN_n denote the edges of this tree.
   - The letter N in N_n was chosen since it will be a miN heap at vec[0].
   - The letters P and N in PN_n were chosen since PN_n will form the edges of
     what will be called the "Pure miN heap".
3) Let FlipCo(i) := total_num_nodes - 1 - i.
4) Define a binary tree X_n on L_n (not necessarily sharing any edges with binary
   tree H_n defined in 2)), rooted at vec[total_num_nodes - 1] in the usual way,
   which we now describe:
   - declare that for any vec[i], vec[j] in L_n, vec[j] is a parent of vec[i]
    in X_n if and only if vec[FlipCo(j)] is a parent of vec[FlipCo(i)] in N_n
    (i.e. in the tree rooted at 0, which was defined in 2)).
   Let PX_n denote denote the edges of this tree on L_n
   - The letter X in X_n was chosen since it will be a maX heap at
     vec[total_num_nodes - 1].
   - The letters P and X in PX_n were chosen since PX_n will form the edges of
     what will be called the "Pure maX heap".
5) We define the "pure min heap" (resp. the "pure max heap"), denoted by PN
   (resp. PX) to be the graph formed by the edges PN_Mid (resp. PX_Mid) where
   Mid = floor(total_num_nodes / 2) if total_num_nodes is odd and
   Mid = total_num_nodes / 2 - 1    if total_num_nodes is even.
   Note that:
   - If total_num_nodes is odd then PN and PX intersect at exactly one node:
     vec[Mid]. While if total_num_nodes is even the PN and PX are disjoint.
   - Every node is contained in some edge of the either PN and/or PX and that
     only node vec[Mid] is potentially contained in both.
6) Definition of Pure Heaps: Notice that the graphs PN and PX are isomorphic
   "mirror images" of each other, which is why we will sometimes call one
   "the mirror" of the other.
   In particular, we will call PX "the mirror pure heap" or "the pure
   max heap" while PN will simply be called "the pure heap" or "the pure min
   heap."
   - The function PrintBiHeap() in biheap_ostream.h will display a biheap (to
     std::cout by default). Looking at any of its output, the top half
     (including the middle node if total_num_nodes is odd) is the "pure (min)
     heap" while the bottom half (which again includes the middle node if
     total_num_nodes is odd) is "the mirror pure heap" or "the pure max heap".
     See also PrintPureMinHeap() and PrintPureMaxHeap() in the biheap_ostream
     class in biheap_ostream.h.
7) Let U_n denote the union of edges PN_n and PX_n.
8) Let HN_n (resp. HX_n) denote the restriction of U_n to F_n (resp. L_n).
   - The letters H and N in PN_n (resp. H and X in PX_n) were chosen since HN_n
     will form the edges of what will be called the "miN Heap" (resp. the
     "maX Heap").
9) Let K be the unique largest integer such that (F_K, HN_K) is a binary tree,
   where we will call this tree "the heap" or "the min heap".
   K can equivalently be defined as the unique largest integer such that
   (L_K, HX_K) is a binary tree,  where we will call this tree "the max heap"
   or "the mirror heap".
   Note that:
   - This unique integer K will often be denoted by: num_nodes_in_heap.
   - The pure min heap (resp. the pure max heap) is always a subgraph of the
     min heap (resp. the max heap). i.e. K >= Mid always holds.
   - Furthermore, if K > Mid if and only if total_num_nodes > 3. That is,
     the pure min heap (resp. the pure max heap) is a strict subgraph of the
     min heap (resp. the max heap) if and only if total_num_nodes > 3.
   - The edges in U_K cover all nodes in vec[0 : total_num_nodes) and that the
     edges of U_K form a connected graph.
   - The number of nodes in the min heap equals the number of nodes in the max
     heap. This common value is the output of the non-trivial function
     GetNumNodesInHeapContainedInBiheap() found in biheap_common.h. This value
     is fundamental to biheaps and will often be denoted by num_nodes_in_heap.
     The formulas defining this value can be found defined in the function
     GetNumNodesInHeapContainedInBiheap() that is in the header biheap_common.h.
10) The graph that the edges U_K form will be called "the biheap graph (of size
    total_num_nodes)"
11) Note that the biheap graph is nothing more than to union of the min heap
    graph and the max heap graph.
    - This is the key observation that determines the BiHeapify() algorithm by
      reducing most of it down to a combination of sifts moving elements
      up/down min/max heaps.

NOTE: There will be nodes that simultaneously belong to BOTH heaps (F_K, HN_K)
 and (L_K, HX_K), which makes the usual definitions of "parent of", "left child
 of", "shift down", etc. ambiguous. In such cases we will qualify the term with
 "(pure) min heap" or "(pure) max heap". Also, the following definitions will
 help remedy this ambiguity.

Definition of "heap coordinate":
 Consider the biheap graph on
  vec[0 : total_num_nodes - 1] of size total_num_nodes, where as before the
  nodes of this graph are vec[0], vec[1], ..., vec[total_num_nodes - 1].
  For all 0 <= i < total_num_nodes, we will call i the "heap coordinate" of
  the node vec[i].
Definition of "mirror coordinate":
 If a node has heap coordinate i then we will call FlipCo(i) the "mirror
  coordinate" of this node, where recall that
   FlipCo(i) := total_num_nodes - 1 - i.
If a variable represents a node's heap (resp. mirror) coordinate then it will
 often be suffixed by _hc (resp. _mc).

The following example of a biheap of size 10 should clarify these definitions.
 Here O's represent nodes in the biheap graph of size 10 and to the right of
  this biheap graph are the nodes labeled in heap coordinates and then to the
  right of that are the nodes labeled in mirror coordinates.
Biheap Graph  Heap Coordinates   Mirror Coordinates
     O               0                   9
  O     O         1     2             8     7
O   O           3   4               6   5

      O   O           5   6               4   3
  O     O         7     8             2     1
     O               9                   0
The following shows the heap coordinates with the mirror coordinates in
 parentheses:
               0(9)
     1(8)                2(7)
3(6)      4(5)

                    5(4)      6(3)
     7(2)                8(1)
               9(0)
The upper half is the pure min heap while the lower half is the pure max heap.
The heap size for this biheap graph of size 10 is 7 so that
 the nodes with heap   coordinates 0, 1, ..., 6 form the min heap while
 the nodes with mirror coordinates 0, 1, ..., 6 form the min heap.
- To draw the min heap (resp. the max heap), cover up all nodes other than nodes
  whose heap (resp. mirror) coordinates are not 0, 1, ..., 6 and then draw lines
  between the remaining 7 nodes as is ordinarily done. Explicitly:
- To draw the min heap, draw the following lines
    0(9) - 1(8), 0(9) - 2(7),
    1(8) - 3(6), 1(8) - 4(5) (we have now drawn the pure min heap)
    2(7) - 5(4), 2(7) - 6(3) (we have now drawn the min heap).
- We know that we have finished drawing the min heap since if, continuing the
  pattern, we were to draw the line 3(6) - 7(2) then upon ignoring
  (a) all nodes other than 0(9), 1(8), ... 7(2), (i.e. ignore 8(1) and 9(0)) and
  (b) all edges adjacent to these ignored nodes
  then we would be left with a graph that has a loop (i.e. a graph that is not
  a binary tree and so cannot, in particular, be the graph of a heap.)
- To draw the max heap, draw the following lines
    9(0) - 8(1), 9(0) - 7(2),
    8(1) - 6(3), 8(1) - 5(4) (we have now drawn the pure max heap)
    7(2) - 4(5), 7(2) - 3(6) (we have now drawn the max heap).
- We know that we have finished drawing the max heap for reasons copmletely
  analogous to how we knew that we finished drawing the min heap (see above).
- Notice the symmetry involved in drawing in the min heap and the max heap.
  This always occurs and it is due to both the definition of the min and max
  heaps as well as to the definitions of heap and mirror coordinates.
Note that the node 0(9) will contain this biheap's minimum value while the node
 9(0) will contain its maximum value.

DEFINITION OF A BIHEAP:
 Given a biheap graph on nodes vec[0, total_num_nodes) (i.e. nodes vec[0], ...,
 vec[total_num_nodes - 1]), then we will call this graph a "biheap" and say that
 it "has the biheap property" if
 (1) on its heap (a subgraph of the biheap graph) these nodes form a min heap
     with a minimum value at vec[0], and
 (2) on its mirror heap (also a subgraph of the biheap graph) these nodes form
     a max heap with a maximum value at vec[total_num_nodes - 1].

Although the above definition of a biheap is relatively complicated, the author
 considers it to be a natural data structure to define. It is also a definition
 that appears to have gone unnoticed until now. There are still many questions
 to be asked about biheaps, including:
 1) Do there exist O(log n) push and pop operations for biheaps?
 2) Is the working hypothesis that BiHeapify() always produces a biheap true?
*/

#ifndef BIHEAPIFY_H_
#define BIHEAPIFY_H_

#include <iostream>

#include "biheap_common.h"

#include "biheapify_even.h"
#include "biheapify_odd.h"

template<class RAI>
void BiHeapify(RAI first, size_type total_num_nodes) {
  if (total_num_nodes % 2 == 0)
    BiHeapifyEven(first, total_num_nodes);
  else
    BiHeapifyOdd(first, total_num_nodes);
}

//This will BiHeapify all nodes in [0, total_num_nodes).
template<class RAI>
inline void BiHeapifySafe(RAI first, size_type total_num_nodes) {
  //If it's small enough that it's easiest to just sort everything.
  if (total_num_nodes < 2)
    return ;

  if (total_num_nodes % 2 == 0)
    BiHeapifySafeEven(first, total_num_nodes);
  else
    BiHeapifySafeOdd(first, total_num_nodes);
  return ;
}

#endif /* BIHEAPIFY_H_ */
