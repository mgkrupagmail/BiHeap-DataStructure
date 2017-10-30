/*
 * biheap_common.h
 *
 *  Created on: Jun 23, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 * This file contains function definitions that are used by other biheapify
 *  files, such as biheapify_even.h, biheapify_odd.h, and biheapify_simple.h.
 */

#ifndef BIHEAP_COMMON_H_
#define BIHEAP_COMMON_H_

#include <algorithm>

#define FLIP(a) (total_num_nodes - 1 - (a))

template<typename size_type = std::size_t>
inline size_type IndexOfLastMinHeapNode(size_type total_num_nodes) {
  return HeapSize(total_num_nodes) - 1;
}

template<typename size_type = std::size_t>
inline size_type PureHeapSize(size_type total_num_nodes) {
  return (total_num_nodes + 1) / 2;
}

template<typename size_type = std::size_t>
inline size_type IndexOfLastNodeInPureHeapGivenPureHeapSize(size_type size_of_pure_heap) {
  return size_of_pure_heap - 1;
}

template<typename size_type = std::size_t>
inline size_type IndexOfLastNodeInPureHeap(size_type total_num_nodes) {
  return IndexOfLastNodeInPureHeapGivenPureHeapSize(PureHeapSize(total_num_nodes));
}

template<typename size_type = std::size_t>
inline size_type NumberOfExtendedInNodes(size_type total_num_nodes, size_type size_of_heap) {
  return total_num_nodes - 2 * (total_num_nodes - size_of_heap);
        //= 2 * size_of_heap - total_num_nodes;// but this may overflow
}

template<typename size_type = std::size_t>
inline size_type NumberOfExtendedInNodes(size_type total_num_nodes) {
  return NumberOfExtendedInNodes(total_num_nodes, HeapSize(total_num_nodes));
}

template<typename size_type = std::size_t>
inline size_type IndexOfFirstExtendedInNode(size_type total_num_nodes, size_type size_of_heap) {
  return total_num_nodes - size_of_heap;
        //= 2 * size_of_heap - total_num_nodes;// but this may overflow
}

template<typename size_type = std::size_t>
inline size_type IndexOfFirstExtendedInNode(size_type total_num_nodes) {
  return IndexOfFirstExtendedInNode(total_num_nodes, HeapSize(total_num_nodes));
}

template<typename size_type = std::size_t>
inline size_type IndexOfFirstLeafGivenPureHeapSize(size_type size_of_pure_heap) {
  return size_of_pure_heap / 2;
}

template<typename size_type = std::size_t>
inline size_type IndexOfFirstLeaf(size_type total_num_nodes) {
  return IndexOfFirstLeaveGivenPureHeapSize(PureHeapSize(total_num_nodes));
}

template<typename size_type = std::size_t>
inline size_type NumberOfLeavesGivenPureHeapSize(size_type size_of_pure_heap) {
  return 2 * ((size_of_pure_heap + 1) / 2);
}

template<typename size_type = std::size_t>
inline size_type NumberOfLeaves(size_type total_num_nodes) {
  return NumberOfLeavesGivenPureHeapSize(PureHeapSize(total_num_nodes));
}

template<typename size_type = std::size_t>
inline size_type GetNumberOfPureMinHeapLeavesWith2MinHeapChildren(size_type total_num_nodes) {
  size_type pure_heap_size = GetSizeOfPureHeap(total_num_nodes);
  return ((pure_heap_size / 2) + 2 * (pure_heap_size % 2) - 1 - (total_num_nodes % 2) - ((total_num_nodes + 2) % 3)) / 3;
  //Which equals (pure_heap_size - (pure_heap_size / 2) - (total_num_nodes % 2) - 1 + (pure_heap_size % 2) - ((total_num_nodes + 2) % 3)) / 3;
}

#undef FLIP
#endif /* BIHEAP_COMMON_H_ */
