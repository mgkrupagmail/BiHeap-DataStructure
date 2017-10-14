//============================================================================
// Name        : main.cpp
// Author      : Matthew Gregory Krupa
// Version     :
// Copyright   : Copyright
// Description : For the definition of a BiHeap, see the comments at the top
//  of biheapify.h.
// This file contains calls to functions that time and test the various
//  BiHeapify functions; these are:
// (1) TimeBiHeapifies(), which times the various BiHeapify functions, and
// (2) MeasureBiHeapifySuccessRate(), which is used to measure how many times
//      a BiHeapify pass (i.e. one of BiHeapifyEven(),
//      BiHeapifyOdd(), or BiHeapifySimpleSinglePass()) must be
//      called before a BiHeap is formed.
//    - Notice that for BiHeapify() this always require only one call to
//      biheapify.
//============================================================================

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "biheapify.h"
#include "biheap_sift.h"

#include "biheap_tikz_graph.h"

#include "biheapify_single_pass_success_rate.h"
#include "biheap_sift_test_correctness.h"

void PrintTikzGraphs() {
  std::vector<int> vec = {0,11,1,12,13,2,3,14,15,16,17,4,5,6,7,18,19,8,9,20,10,21};
  std::cout << GetTikzGraph<int, std::size_t>(22, -6, 6, -3, 3, vec);

  vec = {0,23,1,24,25,2,3,26,27,28,29,4,5,6,7,30,31,32,33,34,35,36,37,8,9,10,11,12,13,14,15,38,39,40,41,16,17,18,19,42,43,20,21,44,22,45};
  std::cout << GetTikzGraph<int, std::size_t>(46, -6, 6, -3, 3, vec);

  std::cout << GetTikzGraph<int, std::size_t>(11, -8, 8, -4, 4);
  std::cout << GetTikzGraph<int, std::size_t>(46, -8, 8, -4, 4);
  std::cout << GetTikzGraph<int, std::size_t>(70, -10, 10, -6, 6);
  return ;
}

int main() {
  long start_total_num_nodes = 1;
  long increment_size        = 1; //Set this to 2 to test vectors that are only
                                 //odd or only even in size.
  //The two functions below will be applied to biheaps of sizes:
  // start_total_num_nodes + (muliple) * increment_size
  // that are <= end_total_num_nodes
  long end_total_num_nodes = static_cast<int>(1u << 10);

  //For each size, the two functions below will biheapify num_vecs_to_try
  // vectors of that size.
  long num_vecs_to_try = static_cast<int>(1u << 18);

  //For MeasureBiHeapifySuccessRate() only. Print information only after you've
  // gone through print_multiple new vector sizes.
  long print_multiple = 1;;
  bool reset_after_print = true; //Don't show a cumulative success and failure
                                 // counts.
  bool verbose = false;

  //This function will go through each of the sizes:
  // start_total_num_nodes, start_total_num_nodes + increment_size,
  // start_total_num_nodes + 2 * increment_size, ...
  // while start_total_num_nodes + # * increment_size
  // remains <= end_total_num_nodes.
  //For each size, it will biheapify num_vecs_to_try vectors of that size.
  //However, it will only display information about the success and failure
  // counts when # = print_multiple, 2 * print_multiple, ....
  //Note that when a vector's size is even then the biheapify success rate
  // on the first try is 100%.
  MeasureBiHeapifySuccessRate<int>(start_total_num_nodes, end_total_num_nodes,
                                   num_vecs_to_try, increment_size,
                                   print_multiple, reset_after_print, verbose);
  BiHeapSiftTestCorrectness<int>(start_total_num_nodes, end_total_num_nodes,
                                 num_vecs_to_try, increment_size);
  PrintTikzGraphs();
  return 0;
}
