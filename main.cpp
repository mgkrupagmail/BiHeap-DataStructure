//============================================================================
// Name        : main.cpp
// Author      : Matthew Gregory Krupa
// Version     :
// Copyright   : MIT License
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

#include "biheapify_single_pass_success_rate.h"
#include "biheapify_time.h"
#include "biheap_sift_test_correctness.h"

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
  long num_vecs_to_try = static_cast<int>(1u << 3);

  //For TimeBiHeapifies() only. For each std::vector that whose biheapification
  // is to be timed, repeat this process num_repititions_per_vec times.
  long num_repititions_per_vec = static_cast<int>(1u << 11);

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
  long long divisor = num_repititions_per_vec; //Set this to 0 to get the the
                             //average time to biheapify each individual vector.
  TimeBiHeapifies<int>(start_total_num_nodes, end_total_num_nodes,
                       num_vecs_to_try, num_repititions_per_vec, increment_size,
                       divisor);
  BiHeapSiftTestCorrectness<int>(start_total_num_nodes, end_total_num_nodes,
                                 num_vecs_to_try, increment_size);
  return 0;
}
