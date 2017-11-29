//============================================================================
// Name        : main.cpp
// Author      : Matthew Gregory Krupa
// Copyright   : Matthew Gregory Krupa
// Description : For the definition of a BiHeap, see the comments at the top
//  of biheapify.h.
// This file contains calls to functions that time and test the various
//  BiHeapify functions; these are:
// (1) MeasureBiHeapifyInwardsPivotProperties(), which measures how good
//      of a pivot value the BiHeapifyInwards() algorithm produces.
//     The most important quantities that are outputted by this function
//      are the smallest value of min(# elements <= pivot, # elements >= pivot),
//      that was encountered and the average of value of
//      min(# elements <= pivot, # elements >= pivot).
// (2) MeasureBiHeapifySuccessRate(), which verifies that our implementation
//      of the BiHeapify algorithm works correctly.
// (3) BiHeapSiftTestCorrectness(), which verifies that our implementation
//      of the BiHeapSift() algorithm works correctly.
// (4) TimeBiHeapifies(), which times the various BiHeapify functions, and
//============================================================================

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "biheapify.h"
#include "biheap_sift.h"

#include "biheapify_single_pass_success_rate.h"
#include "biheap_sift_test_correctness.h"
#include "biheapify_inwards_pivot_testing.h"
#include "biheapify_inwards_and_outwards_testing.h"

#include "biheap_tikz_graph.h"

#include "biheap_ostream.h"
#include "biheapify_time.h"

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
  
  //PrintTikzGraphsExampleCalls();

  MeasureBiHeapifyInwardsAndOutwardsProperties();
  
  //The two key quantities to look at are:
  //smallest min(# elements <= pivot, # elements >= pivot), and
  //average of min(# elements <= pivot, # elements >= pivot)
  MeasureBiHeapifyInwardsPivotProperties();

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
  //For TimeBiHeapifies() only. For each std::vector that whose biheapification
  // is to be timed, repeat this process num_repititions_per_vec times.
  long num_repititions_per_vec = static_cast<int>(1u << 4);
  long long divisor = num_repititions_per_vec; //Set this to 0 to get the the
                             //average time to biheapify each individual vector.
  TimeBiHeapifies<int>(start_total_num_nodes, end_total_num_nodes,
                       num_vecs_to_try, num_repititions_per_vec, increment_size,
                       divisor);
  return 0;
}
