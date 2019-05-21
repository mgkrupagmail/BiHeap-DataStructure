/*
 * biheapify.h
 *
 *  Created on: Jun 12, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */
/* BIHEAP DEFINITION (QUICK):
We're given N elements indexed by 0, ..., N - 1.
These elements form a BiHeap if:
(1) the elements [0, ..., heap_size) form a min heap with the element
    at 0 being the min, and
(2) the elements [N - heap_size, ..., N)
   form a max heap with the element at N - 1 being the max,
 where heap_size is obtained from the function
GetNumNodesInHeapContainedInBiheap(), which is defined in biheap_common.h.
The quantity heap_size is a fundamentally important quantity associated
 with the biheap on N nodes. The formulas that define this
 quantity, although complicated, stem from the intuitive and natural
 graph-theoretic definition of a biheap given in the detailed definition below.
Given an iterator V, the function IsBiheap(), defined in biheap_common.h,
 may be used to check that the V N elements iterated by V
 (i.e. *V, *(V + 1), ..., *(V + N - 1))
 form a biheap.
It would likely be instructive to read the definition of
 IsBiheap(RAI V, SizeType N) found in biheap_common.h
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

There are still many questions to be asked and answered about biheaps, including:
 1) Do there exist O(log n) push and pop operations for biheaps?
 2) How do BiHeaps relate to the median(s) of a set.

About the name: Various types of data structures called double ended heaps
(also known as double ended queues) have already been discovered, but it is
not clear that a biheap can be made into a double ended heap since it
is not clear at the moment whether or not O(log n) pop and push operations
 exist for biheaps. In this sense, the biheap data structure defined here may
not be a specific type of double ended heap; it is nevertheless called a biheap
due to it, by definition, simultaneously being both a min heap and a max heap.

Suppose that we have a positive even number N of nodes
 vec[0 : N) (i.e. nodes vec[0], ..., vec[N - 1])
 and call each vec[i] a "node".
1) Let F_n (L_n) denote the First (resp. Last) n nodes of vec[0 : N)
2) Define a binary tree N_n on F_n rooted at vec[0] in the usual way and let
   PN_n denote the edges of this tree.
   - The letter N in N_n was chosen since it will be a miN heap at vec[0].
   - The letters P and N in PN_n were chosen since PN_n will form the edges of
     what will be called the "Pure miN heap".
3) Let FlipCo(i) := N - 1 - i.
4) Define a binary tree X_n on L_n (not necessarily sharing any edges with binary
   tree H_n defined in 2)), rooted at vec[N - 1] in the usual way,
   which we now describe:
   - declare that for any vec[i], vec[j] in L_n, vec[j] is a parent of vec[i]
    in X_n if and only if vec[FlipCo(j)] is a parent of vec[FlipCo(i)] in N_n
    (i.e. in the tree rooted at 0, which was defined in 2)).
   Let PX_n denote denote the edges of this tree on L_n
   - The letter X in X_n was chosen since it will be a maX heap at
     vec[N - 1].
   - The letters P and X in PX_n were chosen since PX_n will form the edges of
     what will be called the "Pure maX heap".
5) We define the "pure min heap" (resp. the "pure max heap"), denoted by PN
   (resp. PX) to be the graph formed by the edges PN_Mid (resp. PX_Mid) where
   Mid = floor(N / 2) if N is odd and
   Mid = N / 2 - 1    if N is even.
   Note that:
   - If N is odd then PN and PX intersect at exactly one node:
     vec[Mid]. While if N is even the PN and PX are disjoint.
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
     (including the middle node if N is odd) is the "pure (min)
     heap" while the bottom half (which again includes the middle node if
     N is odd) is "the mirror pure heap" or "the pure max heap".
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
   - This unique integer K will often be denoted by: heap_size.
   - The pure min heap (resp. the pure max heap) is always a subgraph of the
     min heap (resp. the max heap). i.e. K >= Mid always holds.
   - Furthermore, if K > Mid if and only if N > 3. That is,
     the pure min heap (resp. the pure max heap) is a strict subgraph of the
     min heap (resp. the max heap) if and only if N > 3.
   - The edges in U_K cover all nodes in vec[0 : N) and that the
     edges of U_K form a connected graph.
   - The number of nodes in the min heap equals the number of nodes in the max
     heap. This common value is the output of the non-trivial function
     GetNumNodesInHeapContainedInBiheap() found in biheap_common.h. This value
     is fundamental to biheaps and will often be denoted by heap_size.
     The formulas defining this value can be found defined in the function
     GetNumNodesInHeapContainedInBiheap() that is in the header biheap_common.h.
10) The graph that the edges U_K form will be called "the biheap graph (of size
    N)"
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
  vec[0 : N - 1] of size N, where as before the
  nodes of this graph are vec[0], vec[1], ..., vec[N - 1].
  For all 0 <= i < N, we will call i the "heap coordinate" of
  the node vec[i].
Definition of "mirror coordinate":
 If a node has heap coordinate i then we will call FlipCo(i) the "mirror
  coordinate" of this node, where recall that
   FlipCo(i) := N - 1 - i.
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
 Given a biheap graph on nodes vec[0, N) (i.e. nodes vec[0], ...,
 vec[N - 1]), then we will call this graph a "biheap" and say that
 it "has the biheap property" if
 (1) on its heap (a subgraph of the biheap graph) these nodes form a min heap
     with a minimum value at vec[0], and
 (2) on its mirror heap (also a subgraph of the biheap graph) these nodes form
     a max heap with a maximum value at vec[N - 1].

Although the above definition of a biheap is relatively complicated, the author
 considers it to be a natural data structure to define. It is also a definition
 that appears to have gone unnoticed until now. There are still many questions
 to be asked about biheaps, including:
 1) Do there exist O(log n) push and pop operations for biheaps?
*/

#ifndef BIHEAPIFY_H_
#define BIHEAPIFY_H_

#include <algorithm>

//#define FLIP(a) ((N - 1) - (a))

/*
 * ================== START: Definition of BiHeapify ====================
 */

template<typename SizeType = std::size_t>
inline SizeType LeftChild(SizeType node) { return (2 * node) + 1; }

template<typename SizeType = std::size_t>
inline SizeType RightChild(SizeType node) {return 2 * (node + 1); }

//Note that even if SizeType is unsigned, Parent(0) == 0.
template<typename SizeType = std::size_t>
inline SizeType Parent(SizeType node) {
  return static_cast<SizeType>((static_cast<long long>(node) - 1) / 2);
}

//Unlike Parent(), this assumes that node > 0 (hence "NotRoot").
//Note that if SizeType is unsigned, then ParentNotRoot(0) != 0,
// which is the reason for assuming that node > 0.
//If node > 0 then the static_cast<long long> found in Parent()
// can be avoided.
template<typename SizeType = std::size_t>
inline SizeType ParentNotRoot(SizeType node) {
  return (node - 1) / 2;
}

template<typename SizeType = std::size_t>
inline SizeType HeapSize(SizeType N) {
  return N - static_cast<SizeType>(N / 3);
}

template<typename SizeType = std::size_t>
inline SizeType IndexOfLastMinHeapNodeGivenHeapSize(SizeType heap_size) {
  return heap_size - 1;
}

//Assumes that pos_mc is a node in the max heap.
//The node with max heap coordinate last_node_in_biheap_mc is a
// node in the BiHeap constructed so far such that if v is any
// node whose max heap coordinate is < last_node_in_biheap_mc,
// then v does NOT belong to the BiHeap constructed so far.
template<class RAI, typename SizeType = std::size_t>
inline void SiftUpMaxHeapMC(RAI V, SizeType N, SizeType N_minus1, SizeType pos_mc,
                            SizeType last_node_in_biheap_mc) {
  SizeType parent_mc;
  if (pos_mc == 0 ||
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) < last_node_in_biheap_mc)
    return ;
  auto pos_it    = V + (N_minus1 - pos_mc);
  auto pos_value = *pos_it;
  do {
    auto parent_it = V + (N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent_mc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
     (parent_mc = ParentNotRoot<SizeType>(pos_mc)) >= last_node_in_biheap_mc);
  return ;
}

//Assumes that pos_hc is a node in the min heap.
template<class RAI, typename SizeType = std::size_t>
inline void SiftUpMinHeapHC(RAI V, SizeType pos_hc,
                            SizeType first_node_in_biheap_hc) {
  SizeType parent_hc;
  if (pos_hc == 0 ||
      (parent_hc = ParentNotRoot<SizeType>(pos_hc)) < first_node_in_biheap_hc)
    return ;
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto parent_it = V + parent_hc;
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent_hc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 &&
        (parent_hc = ParentNotRoot<SizeType>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

//Assumes that the node pos_hc belongs to the min heap and that
// pos_hc <= last_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifySiftFromMinToMax(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_hc,
                                      SizeType last_node_in_biheap_hc) {
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto left_child_hc  = LeftChild<SizeType>(pos_hc);
    auto right_child_hc = left_child_hc + 1;
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;

    //assert((left_child_hc  < heap_size) && (left_child_hc  <= last_node_in_biheap_hc));
    bool is_right_child_valid = right_child_hc <= last_node_in_biheap_hc &&
                                right_child_hc < heap_size;
    decltype(pos_value) smaller_value = *left_it, right_value;
    RAI smaller_it;
    if (is_right_child_valid && (right_value = *right_it) < smaller_value) {
      smaller_value = right_value;
      smaller_it    = right_it;
      pos_hc        = right_child_hc;
    } else {
      smaller_it    = left_it;
      pos_hc        = left_child_hc;
    }
    if (smaller_value < pos_value) {
      std::iter_swap(pos_it, smaller_it);
      pos_it        = smaller_it;
    } else
      return ;
  } while (pos_hc < first_in_node) ;
  //Start sifting up the max heap.
  //At this point pos is an In node.
  auto pos_mc                 = N_minus1 - pos_hc;
  auto last_node_in_biheap_mc = N_minus1 - last_node_in_biheap_hc;
  SizeType parent_mc          = ParentNotRoot<SizeType>(pos_mc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_mc == 0 comparison since N > 8 implies this pos_mc > 0.
  if (pos_mc == 0 || parent_mc < last_node_in_biheap_mc)
    return ;
  do {
    auto parent_it = V + (N_minus1 - parent_mc);
    if (*parent_it < pos_value) {
      std::iter_swap(pos_it, parent_it);
      pos_mc = parent_mc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_mc > 0 &&
      (parent_mc = ParentNotRoot<SizeType>(pos_mc)) >= last_node_in_biheap_mc);
  return ;
}

//Assumes that the node pos_mc belongs to the max heap and that
// FLIP(pos_mc) >= first_node_in_biheap_hc.
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapifySiftFromMaxToMin(RAI V, SizeType N, SizeType N_minus1,
                                      SizeType heap_size,
                                      SizeType first_in_node,
                                      SizeType pos_mc,
                                      SizeType first_node_in_biheap_hc) {
  auto pos_hc    = N_minus1 - pos_mc;
  auto pos_it    = V + pos_hc;
  auto pos_value = *pos_it;
  do {
    auto left_child_mc  = LeftChild<SizeType>(pos_mc);
    auto right_child_mc = left_child_mc + 1; //= RightChild(pos_mc);
    auto left_child_hc  = N_minus1 - left_child_mc;//= pos_hc - pos_mc - 1
    auto right_child_hc = left_child_hc - 1; //= FLIP(right_child_mc)
    auto left_it        = V + left_child_hc;
    auto right_it       = V + right_child_hc;

    //assert((left_child_mc < heap_size) && (left_child_hc >= first_node_in_biheap_hc) && (right_child_hc >= first_node_in_biheap_hc));
    bool is_right_child_valid = right_child_mc < heap_size;
    decltype(pos_value) larger_value = *left_it, right_value;
    RAI larger_it;
    if (is_right_child_valid && larger_value < (right_value = *right_it)) {
      larger_value = right_value;
      larger_it    = right_it;
      pos_hc       = right_child_hc;
      pos_mc       = right_child_mc;
    } else {
      larger_it    = left_it;
      pos_hc       = left_child_hc;
      pos_mc       = left_child_mc;
    }
    if (pos_value < larger_value){
      std::iter_swap(pos_it, larger_it);
      pos_it       = larger_it;
    }
    else
      return ;
  } while (pos_mc < first_in_node) ;
  //Start sifting up the min heap.
  //At this point pos is an In node.
  SizeType parent_hc = ParentNotRoot<SizeType>(pos_hc);
  //Note: If you initially sort all BiHeaps of size N < 9 and then return,
  // then you can skip this pos_hc == 0 comparison since N > 8 implies this pos_hc > 0.
  if (pos_hc == 0 || parent_hc < first_node_in_biheap_hc)
    return ;
  do {
    auto parent_it = V + parent_hc;
    if (pos_value < *parent_it) {
      std::iter_swap(pos_it, parent_it);
      pos_hc = parent_hc;
      pos_it = parent_it;
    } else {
      return ;
    }
  } while (pos_hc > 0 &&
     (parent_hc = ParentNotRoot<SizeType>(pos_hc)) >= first_node_in_biheap_hc);
  return ;
}

/* This will BiHeapify all nodes in [0, N).
 */
/*
 * Remark:
 *  This algorithm has complexity O(N). To see why, recall the
 *   argument showing that the heapify operation has O(n) complexity (e.g. as
 *   found on pp. 115 - 116 of "The Algorithm Design Manual" 2nd edition); this
 *   argument generalizes to prove that this algorithm also runs in O(n) times.
 *  The key observation is that the biheap is constructed by adding one node
 *   at a time with this node alternating between a node in a min heap and a
 *   node in the max heap. The complexity of this biheapification is easily
 *   seen to be twice the complexity of the above mentioned heapification
 *   operation plus a constant.
 */
template<class RAI, typename SizeType = std::size_t>
inline void BiHeapify(RAI V, SizeType N) {
  if(N < 2)
    return ;
  SizeType heap_size     = HeapSize(N);
  SizeType first_in_node = N - heap_size;
  //Ignore all extended in arrows, unless N % 3 == 2, in which
  // case ignore all but the middle two extended in arrows.
  SizeType last_node_in_biheap_hc  = IndexOfLastMinHeapNodeGivenHeapSize(heap_size)
                                     - (N % 3 == 2);
  SizeType N_minus1                = N - 1;
  SizeType first_node_in_biheap_hc = N_minus1 - last_node_in_biheap_hc;
  //NOTE: It's not necessary to pass N_minus1, heap_size, or first_in_node to
  // these functions (since they can be computed from N) but we do so since
  // these variables won't actually be placed on the stack each time due to these
  // functions being both inlined and templates. Any half decent optimizer would
  // optimize away allocating stack space and copy these values.
  while (first_node_in_biheap_hc > 0) {
    BiHeapifySiftFromMinToMax<RAI, SizeType>(V, N, N_minus1, heap_size, first_in_node,
       --first_node_in_biheap_hc, last_node_in_biheap_hc);
    BiHeapifySiftFromMaxToMin<RAI, SizeType>(V, N, N_minus1, heap_size, first_in_node,
        N_minus1 - (++last_node_in_biheap_hc), first_node_in_biheap_hc);
  }
  return ;
}

/*
 * ================== END: Definition of BiHeapify ====================
 */

/*
 * ================== START: Check if nodes form a BiHeap ====================
 */

/* Checks whether or not the V total_um_nodes given given the iterator
 *  V define a BiHeap.
 */
template<class RAI, typename SizeType = std::size_t>
bool IsBiHeap(RAI V, SizeType N) {
  if (N <= 1)
    return true;
  SizeType heap_size = HeapSize<SizeType>(N);
  //Check that the nodes V, ..., V + heap_size - 1 form
  // a min heap with min at V. This is half of the biheap condition.
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + i);
      //Check that the parent and left child satisfy the min heap condition.
      if (parent_value > *(V + (right_child - 1)))
        return false;

      //Check that the parent and right child satisfy the min heap condition.
      if (parent_value > *(V + right_child))
        return false;
    }
    //If the min heap's last element is an only child then check that it and
    // its parent satisfy the min heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < heap_size
          && *(V + i) > *(V + left_child))
        return false;
    }
  }
  //Check that the nodes FLIP(V), ...,
  // FLIP(V + heap_size - 1) form a max heap with max
  // at V + N - 1.
  SizeType N_minus1 = N - 1;
  {
    SizeType i = 0;
    for (SizeType right_child; (right_child = RightChild<SizeType>(i))
                                                    < heap_size; i++) {
      auto parent_value = *(V + (N_minus1 - i));
      SizeType mirror_left_child_hc = N_minus1 - (right_child - 1);
      //Check that the parent and left child satisfy the max heap condition.
      if (parent_value < *(V + mirror_left_child_hc))
        return false;

      //right_child_hc = FLIP(right_child_mc) = FLIP(left_child + 1)
      // = N - 1 - (left_child + 1) = mirror_left_child_hc - 1.
      //Check that the parent and right child satisfy the max heap condition.
      if (parent_value < *(V + (mirror_left_child_hc - 1)))
        return false;
    }
    //If the max heap's last element is an only child then check that it and
    // its parent satisfy the max heap condition (i.e. the biheap condition).
    {
      SizeType left_child;
      if ((left_child = LeftChild<SizeType>(i)) < heap_size
          && *(V + (N_minus1 - i)) < *(V + (N_minus1 - left_child)))
        return false;
    }
  }
  return true;
}

/*
 * ================== END: Check if nodes form a BiHeap ====================
 */

//#undef FLIP

#endif /* BIHEAPIFY_H_ */
