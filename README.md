# BiHeap
A BiHeap is a new independently developed data structure a part of whose data is lies a min heap, a part that lies in a max heap, and a part that lies simultaneously in both. A BiHeap is a doubled ended heap, having both a min and a max that are readily accessible. An exposition on BiHeaps, including their definition, can be found in "BiHeaps and Pivot Selection.pdf". We use BiHeaps (and a derived data structure call "almost BiHeaps") to define a double ended priority queue we call a BiQueue that has amortized O(log N) insertions and pops (of the min and max elements) and can be constructed in O(N) time. 

For the implementation of a double ended priority queue with both a min and a max and amortized O(log N) inserts and pops see the class BiQueue in the file biqueue.h

In biheap_select.h we define BiHeapSelect(), which does the same thing as the QuickSelect algorithm. However, unlike the Quickselect algorithm, BiHeapSelect() is O(N) (this is something that I'm still in the middle of proving mathematically). Empirical testing by the functions in biheap_select_measure_complexity.h show that BiHeapSelect() emplaces the element in the desired position using at most 15 * N swap operations.

This project, among other things, implements the BiHeapify() function which given a random access iterator and a length, forms a BiHeap. 
This new data structure:
<BR>
1) is not excessively computationally costly to create, and<BR>
2) is not an excessively complicated operation, being not much more difficult to implement than a min heap and a max heap; indeed the BiHeapify() operation is primarily a combination of the four standard sifting operations, which are sifting elements up/down in a min/max heap.

These properties make the BiHeap definition and the BiHeapify operation easier to implement than double-ended heaps, its nearest counterpart, while also having good performance. The author considers the definition of a biheap to be a very natural data structure to define. It is also a definition that appears to have gone unnoticed until now.

The definition of a BiHeap is highly detailed but we can give a quick summary of it. We're given N elements indexed by 0, ..., N - 1. These elements form a BiHeap if:

 1) the elements [0, ..., HeapSize(N)) form a min heap with the element at 0 being the min, and<BR>
 2) the elements [N - HeapSize(N), ..., N) form a max heap with the element at N - 1 being the max,
 
where HeapSize(N) is defined in "BiHeaps and Pivot Selection.pdf".

The quantity HeapSize(N) is a fundamentally important quantity associated with the biheap on total_num_nodes nodes. The formulas that define this quantity stems from the intuitive and natural graph-theoretic definition of a BiHeap given in detail in "BiHeaps and Pivot Selection.pdf".

Besides making the minimum and maximum elements readily available in O(N) time, the author has noticed that, as expected, the BiHeapify() functions tends to take randomly distributed elements (following a uniform distribution) and reorders them so that many of the larger elements are near the max and many of the smaller elements are near the min and the value in middle of the data to be close to the median. In fact, we define a function called BiHeapifyInwards() that produces a better-than-random pivot to help to speed up certain sorting algorithms such as quicksort.

About the name: Various types of data structures called double ended heaps (also known as double ended queues) have already been discovered, but it is not clear that a BiHeap can be made into a double ended heap (although this is possible with almost BiHeaps) since it is not clear at the moment whether or not O(log N) pop and push operations exist for biheaps. In this sense, the BiHeap data structure defined here may not be a specific type of double ended heap; it is nevertheless called a BiHeap due to it, by definition, simultaneously being both a min heap and a max heap.<BR>

There are still many questions to be asked and answered about biheaps, including:
 1) Do there exist O(log N) push and pop operations for BiHeap?
 2) How do BiHeaps relate to the median of a set?

Copyright Matthew Gregory Krupa
