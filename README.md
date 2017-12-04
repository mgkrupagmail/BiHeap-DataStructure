# BiHeap
A BiHeap is a new independently developed data structure a part of whose data is lies a min heap, a part that lies in a max heap, and a part that lies simultaneously in both. An exposition on BiHeaps, including their definition, can be found in "BiHeaps and Pivot Selection.pdf".The minimum and maximum elements in a BiHeap are readily accessible. We use BiHeaps (and a derived data structure call "almost BiHeaps") to define a double ended queue, called a BiQueue, that has amortized O(log N) insertions and pops (of the min and max elements).  In biheapify_inwards_and_outwards.h, we give two O(N) functions that based on empircal testing (see biheapify_inwards_and_outwards_testing.h) always simultaneously both emplace the median in the middle of a vector and partition the elements so that everything smaller (resp. larger) than the median lies to the left (resp. right) of the median; it has not been shown theoretically that this always happens but on the billions of random data sets that the author has tested so far, this has always been the case. What the author has shown theoretically (see "BiHeaps and Pivot Selection.pdf") is that that the BiHeapifyInwardsNicerMath() function finds a pivot value with at least 2^(log_3(N)) elements above it and at least this amount below it (although the author recommends using BiHeapifyInwards() instead of BiHeapifyInwardsNicerMath()).

This project, among other things, implements the BiHeapify() function which given a random access iterator and a length, forms a BiHeap. 
This new data structure:
<BR>
1) is not excessively computationally costly to create, and<BR>
2) is not an excessively complicated operation, being not much more difficult to implement than a min heap and a max heap; indeed the BiHeapify() operation is primarily a combination of the four standard sifting operations, which are sifting elements up/down in a min/max heap.

These properties make the BiHeap definition and the BiHeapify operation easier to implement than double-ended heaps, its nearest counterpart, while also having good performance. The author considers the definition of a biheap to be a very natural data structure to define. It is also a definition that appears to have gone unnoticed until now.

The definition of a biheap is highly detailed but we can give a quick summary of it. We're given N elements indexed by 0, ..., N - 1. These elements form a BiHeap if:

 1) the elements [0, ..., num_nodes_in_heap) form a min heap with the element at 0 being the min, and<BR>
 2) the elements [total_num_nodes - num_nodes_in_heap, ..., N) form a max heap with the element at N - 1 being the max,
 
where num_nodes_in_heap is defined in "BiHeaps and Pivot Selection.pdf".

The quantity num_nodes_in_heap is a fundamentally important quantity associated with the biheap on total_num_nodes nodes. The formulas that define this quantity stem from the intuitive and natural graph-theoretic definition of a biheap given in detail in the comments at the top of biheapify.h.

Besides making the minimum and maximum elements readily available in O(N) time, the author has noticed that, as expected, the BiHeapify() functions tends to take randomly distributed elements (following a uniform distribution) and reorders them so that many of the larger elements are near the max and many of the smaller elements are near the min and the value in middle of the data to be close to the median. In fact, we define a function called BiHeapifyInwards() that produces a better-than-random pivot to help to speed up certain sorting algorithms such as quicksort.

About the name: Various types of data structures called double ended heaps (also known as double ended queues) have already been discovered, but it is not clear that a biheap can be made into a double ended heap since it is not clear at the moment whether or not O(log N) pop and push operations exist for biheaps. In this sense, the biheap data structure defined here may not be a specific type of double ended heap; it is nevertheless called a biheap due to it, by definition, simultaneously being both a min heap and a max heap.<BR>

There are still many questions to be asked and answered about biheaps, including:
 1) Do there exist O(log N) push and pop operations for biheaps?
 2) How do BiHeaps relate to the median of a set?

Copyright Matthew Gregory Krupa
