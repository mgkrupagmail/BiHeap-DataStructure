'''
Created on Jul 16, 2017

@author: Matthew Gregory Krupa

Excluding IsBiHeap(), the functions given below  
are the minimum needed to biheapify a collection
of comparable objects.
For the commented versions of these functions 
containing explanations of the source code, 
see the C++ files biheapify_common.h, biheapify_even.h
biheapify_odd.h, and biheapify.h.
'''

def GetLeftChild(node):
    return (2 * node) + 1

def GetRightChild(node):
    return 2 * (node + 1)

def GetParentNotRoot(node):
    return (node - 1) // 2

def GetNumNodesInHeap(total_num_nodes):
    if total_num_nodes <= 2:
        return total_num_nodes
    if total_num_nodes % 2  == 0:
        i = total_num_nodes // 2
        return 4*((i-(i%3))//3) + (i%3) + 1 - (((i+2)%3)//2)
    n = total_num_nodes % 6
    total = 4*(total_num_nodes // 6) + (n + 1)//2
    if n == 5:
        return total + 1
    return total

#FlipCo is short for Flip Coordinate to (resp. from) an hc coordinate from (resp. to) an mc coordinate
def FlipCo(co, total_num_nodes):
    return total_num_nodes - 1 - co

def SiftUpMinHeapHC(li, pos_hc, smallest_node_in_biheap_hc):
    if pos_hc <= 0:
        return 
    parent = GetParentNotRoot(pos_hc)
    if parent < smallest_node_in_biheap_hc:
        return
    while True:
        if li[pos_hc] < li[parent]:
            li[pos_hc], li[parent] = li[parent], li[pos_hc]
            pos_hc = parent
        else:
            return
        if pos_hc <= 0:
            return
        parent = GetParentNotRoot(pos_hc)
        if parent < smallest_node_in_biheap_hc:
            return
    return

def SiftUpMaxHeapMC(li, total_num_nodes, pos_mc, smallest_node_in_biheap_mc):
    if pos_mc <= 0:
        return 
    parent_mc = GetParentNotRoot(pos_mc)
    if parent_mc < smallest_node_in_biheap_mc:
        return
    pos_hc = FlipCo(pos_mc, total_num_nodes)
    parent_hc = FlipCo(parent_mc, total_num_nodes)
    while True:
        if li[pos_hc] > li[parent_hc]:
            li[pos_hc], li[parent_hc] = li[parent_hc], li[pos_hc]
            pos_mc = parent_mc
            pos_hc = parent_hc
        else:
            return
        if pos_mc <= 0:
            return
        parent_mc = GetParentNotRoot(pos_mc)
        if parent_mc < smallest_node_in_biheap_mc:
            return
        parent_hc = FlipCo(parent_mc, total_num_nodes)
    return

def SiftUpMaxHeapHC(li, total_num_nodes, pos_hc, smallest_node_in_biheap_mc):
    SiftUpMaxHeapMC(li, total_num_nodes, FlipCo(pos_hc, total_num_nodes), smallest_node_in_biheap_mc)
    return 

def SiftFromMinToMaxEven(li, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap, pos_hc, largest_node_in_biheap_hc):
    while pos_hc < first_node_in_mirror_heap:
        left_child   = GetLeftChild(pos_hc)
        right_child  = left_child + 1
        smaller      = None
        is_right_child_valid = (right_child <= largest_node_in_biheap_hc) and (right_child < num_nodes_in_heap)
        if is_right_child_valid and (li[right_child] < li[left_child]):
            smaller = right_child
        else:
            smaller = left_child
        if li[pos_hc] > li[smaller]:
            li[pos_hc], li[smaller] = li[smaller], li[pos_hc]
            pos_hc = smaller
        else:
            return
    SiftUpMaxHeapHC(li, total_num_nodes, pos_hc, FlipCo(largest_node_in_biheap_hc, total_num_nodes))
    return

def SiftFromMaxToMinEven(li, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap, pos_mc, smallest_node_in_biheap_hc):
    pos_hc = FlipCo(pos_mc, total_num_nodes)
    while pos_mc < first_node_in_mirror_heap:
        left_child_mc  = GetLeftChild(pos_mc)
        right_child_mc = left_child_mc + 1
        left_child_hc  = FlipCo(left_child_mc, total_num_nodes)
        right_child_hc = left_child_hc - 1
        larger_hc      = None
        is_right_child_valid = (right_child_hc >= smallest_node_in_biheap_hc) and (right_child_mc < num_nodes_in_heap)
        if is_right_child_valid and (li[right_child_hc] > li[left_child_hc]):
            larger_hc = right_child_hc
            pos_mc = right_child_mc
        else:
            larger_hc = left_child_hc
            pos_mc = left_child_mc
        if li[pos_hc] < li[larger_hc]:
            li[pos_hc], li[larger_hc] = li[larger_hc], li[pos_hc]
            pos_hc = larger_hc
        else:
            return
    SiftUpMinHeapHC(li, pos_hc, smallest_node_in_biheap_hc)
    return
        
def SiftFromMinToMaxOdd(li, total_num_nodes, num_nodes_in_heap, pos_hc, largest_node_in_biheap_hc):
    while pos_hc <= total_num_nodes // 2:
        left_child   = GetLeftChild(pos_hc)
        right_child  = left_child + 1
        smaller      = None
        is_left_child_valid = (left_child <= largest_node_in_biheap_hc) and (left_child < num_nodes_in_heap)
        if is_left_child_valid == False:
            break
        is_right_child_valid = (right_child <= largest_node_in_biheap_hc) and (right_child < num_nodes_in_heap)
        if is_right_child_valid and (li[right_child] < li[left_child]):
            smaller = right_child
        else:
            smaller = left_child
        if li[pos_hc] > li[smaller]:
            li[pos_hc], li[smaller] = li[smaller], li[pos_hc]
            pos_hc = smaller
        else:
            return
    SiftUpMaxHeapHC(li, total_num_nodes, pos_hc, FlipCo(largest_node_in_biheap_hc, total_num_nodes))
    return

def SiftFromMaxToMinOdd(li, total_num_nodes, num_nodes_in_heap, pos_mc, smallest_node_in_biheap_hc):
    pos_hc = FlipCo(pos_mc, total_num_nodes)
    while pos_mc <= total_num_nodes // 2:
        left_child_mc  = GetLeftChild(pos_mc)
        right_child_mc = left_child_mc + 1
        left_child_hc  = FlipCo(left_child_mc, total_num_nodes)
        right_child_hc = left_child_hc - 1
        larger_hc      = None
        is_left_child_valid = (left_child_hc >= smallest_node_in_biheap_hc) and (left_child_mc < num_nodes_in_heap)
        if is_left_child_valid == False:
            break
        is_right_child_valid = (right_child_hc >= smallest_node_in_biheap_hc) and (right_child_mc < num_nodes_in_heap)
        if is_right_child_valid and (li[right_child_hc] > li[left_child_hc]):
            larger_hc = right_child_hc
            pos_mc = right_child_mc
        else:
            larger_hc = left_child_hc
            pos_mc = left_child_mc
        if li[pos_hc] < li[larger_hc]:
            li[pos_hc], li[larger_hc] = li[larger_hc], li[pos_hc]
            pos_hc = larger_hc
        else:
            return
    SiftUpMinHeapHC(li, pos_hc, smallest_node_in_biheap_hc)
    return

def BiHeapifyEven(li, total_num_nodes):
    if total_num_nodes < 2:
        return
    num_nodes_in_heap = GetNumNodesInHeap(total_num_nodes)
    first_node_in_mirror_heap = total_num_nodes - num_nodes_in_heap
    smallest_node_in_biheap_hc = total_num_nodes // 2
    largest_node_in_biheap_hc = smallest_node_in_biheap_hc - 1
    
    while smallest_node_in_biheap_hc > 0:
        smallest_node_in_biheap_hc = smallest_node_in_biheap_hc - 1
        SiftFromMinToMaxEven(li, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap, smallest_node_in_biheap_hc, largest_node_in_biheap_hc)

        largest_node_in_biheap_hc = largest_node_in_biheap_hc + 1
        #Note that FlipCo(largest_node_in_biheap_hc, total_num_nodes) = smallest_node_in_biheap_hc
        SiftFromMaxToMinEven(li, total_num_nodes, num_nodes_in_heap, first_node_in_mirror_heap, smallest_node_in_biheap_hc, smallest_node_in_biheap_hc)
    return

def BiHeapifyOdd(li, total_num_nodes):
    if total_num_nodes < 2:
        return
    num_nodes_in_heap = GetNumNodesInHeap(total_num_nodes)
    smallest_node_in_biheap_hc = (total_num_nodes // 2) + 1
    largest_node_in_biheap_hc = smallest_node_in_biheap_hc
    
    while smallest_node_in_biheap_hc > 0:
        smallest_node_in_biheap_hc = smallest_node_in_biheap_hc - 1
        SiftFromMinToMaxOdd(li, total_num_nodes, num_nodes_in_heap, smallest_node_in_biheap_hc, largest_node_in_biheap_hc)
        if largest_node_in_biheap_hc < total_num_nodes - 1:
            largest_node_in_biheap_hc = largest_node_in_biheap_hc + 1
            SiftFromMaxToMinOdd(li, total_num_nodes, num_nodes_in_heap, FlipCo(largest_node_in_biheap_hc, total_num_nodes), smallest_node_in_biheap_hc)
    return

def BiHeapify(li, total_num_nodes):
    if total_num_nodes % 2 == 0:
        BiHeapifyEven(li, total_num_nodes)
    else:
        BiHeapifyOdd(li, total_num_nodes)
        
def IsBiheap(li, total_num_nodes):
    if total_num_nodes <= 3:
        if total_num_nodes <= 1:
            return True
        elif total_num_nodes == 2:
            return li[0] <= li[1]
        else:
            return li[0] <= li[1] and li[1] <= li[2]
    num_nodes_in_heap = GetNumNodesInHeap(total_num_nodes)
    #Check that the first num_nodes_in_heap nodes form a min heap.
    i = 0
    right_child = GetRightChild(i)
    while right_child < num_nodes_in_heap:
        if li[i] > li[right_child - 1]:
            return False
        if li[i] > li[right_child]:
            return False
        i = i + 1
        right_child = GetRightChild(i)
    #if the min heap has a single left child then check that it satisfies the min heap condition.
    left_child = GetLeftChild(i)
    if (left_child < num_nodes_in_heap and li[i] > li[left_child]):
        return False
    #Check that the last num_nodes_in_heap nodes form a max heap.
    i_mc = 0
    right_child_mc = GetRightChild(i_mc)
    while right_child_mc < num_nodes_in_heap:
        i_hc = FlipCo(i_mc, total_num_nodes)
        right_child_hc = FlipCo(right_child_mc, total_num_nodes)
        left_child_hc = FlipCo(right_child_mc - 1, total_num_nodes)
        if li[i_hc] < li[left_child_hc]:
            return False
        if li[i_hc] < li[right_child_hc]:
            return False
        i_mc = i_mc + 1
        right_child_mc = GetRightChild(i_mc)
    #if the max heap has a single left child then check that it satisfies the max heap condition.
    left_child_mc = GetLeftChild(i_mc)
    left_child_hc = FlipCo(left_child_mc, total_num_nodes)
    i_hc = FlipCo(i_mc, total_num_nodes)
    if (left_child_mc < num_nodes_in_heap and li[i_hc] < li[left_child_hc]):
        return False
    return True


#Example using BiHeapify():    
import random

def BiheapifyTestCorrectness():
    list_size_start = 1
    list_size_end = 2**15
    list_size_increment = 1
    list_size = list_size_start
    num_random_lists_per_vec_size = 2**8
    while list_size <= list_size_end:
        for _ in range(num_random_lists_per_vec_size):
            li = [random.random() for _ in range(list_size)]
            BiHeapify(li, list_size)
            if IsBiheap(li, list_size) == False:
                print('List of size ', list_size, ' is not a biheap:')
                print(li)
        list_size += list_size_increment
        print('All', num_random_lists_per_vec_size, 'lists of size list_size =', 
              list_size, 'have been successfully biheapified.')
    return 

BiheapifyTestCorrectness()

            
