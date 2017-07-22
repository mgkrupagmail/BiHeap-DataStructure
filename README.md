# BiHeap
A biheap is a new independently developed data structure, with a part of its data in a min heap, a part in a max heap, and a part that is simultaneously in both. The minimum and maximum elements in a BiHeap are readily accessible. This project implements the BiHeapify() function which given a random access iterator and a length, forms a biheap. See the comment at the top of biheapify.h for the definition of a BiHeap.
This new data structure:
<BR>
1) is not excessively computationally costly to create, especially when there are an even number of elements being considered, and<BR>
2) is not an excessively complicated operation, being not much more difficult to implement than a min heap and a max heap; indeed the BiHeapify() operation is primarily a combination of the four standard sifting operations, which are sifting elements up/down in a min/max heap operations.

These properties make the biheap definition and the biheapify operation easier to implement than double-ended heaps, its nearest counterpart, while also having good performance. The author considers the definition of a biheap to be a very natural data structure to define. It is also a definition that appears to have gone unnoticed until now.

Besides making the minimum and maximum elements readily available in O(n) time (based on experimentation), the author has noticed that, as expected, the BiHeapify() functions tends to take randomly distributed elements (following a uniform distribution) and reorders them so that many of the larger elements are near the max and many of the smaller elements are near the min and the value in middle of the data to be close to the median. Therefore, it is reasonable to expect that applying BiHeapify to a collection of data before sorting it may help to speed up certain sorting algorithms such as quicksort.

There are still many questions to be asked and answered about biheaps, including:
 1) Do there exist O(log n) push and pop operations for biheaps?<BR>
 2) Is the working hypothesis (obtained from emperical testing) that BiHeapify() always produces a biheap true?<BR>
 3) Is BiHeapifyOdd() an O(n) operation? If not then is it an O(n log n) operation?<BR>
 4) Does there exist any O(n) function that biheapifies odd-sized data.

About the name: Various types of data structures called double ended heaps (also known as double ended queues) have already been discovered, but it is not clear that a biheap can be made into a double ended heap since it is not clear at the moment whether or not O(log n) pop and push operations exist for biheaps. In this sense, the biheap data structure defined here may not be a specific type of double ended heap; it is nevertheless called a biheap due to it, by definition, simultaneously being both a min heap and a max heap.<BR>

The definition of a biheap is highly detailed but we can give a quick summary of it. We're given total_num_nodes elements indexed by 0, ..., total_num_nodes - 1. These elements form a BiHeap if:

 1) the elements [0, ..., num_nodes_in_heap) form a min heap with the element at 0 being the min, and<BR>
 2) the elements [total_num_nodes - num_nodes_in_heap, ..., total_num_nodes) form a max heap with the element at total_num_nodes - 1 being the max,
 
where num_nodes_in_heap is obtained from the function GetNumNodesInHeapContainedInBiheap(), which is defined in biheap_common.h. 

The quantity num_nodes_in_heap is a fundamentally important quantity associated with the biheap on total_num_nodes nodes. The formulas that define this quantity, although complicated, stem from the intuitive and natural graph-theoretic definition of a biheap given in detail in the comments at the top of biheapify.h.
