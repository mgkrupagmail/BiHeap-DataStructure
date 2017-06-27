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
//      a BiHeapify pass (i.e. one of BiHeapifyEvenSinglePass(),
//      BiHeapifyOddSinglePass(), or BiHeapifySimpleSinglePass()) must be
//      called before a BiHeap is formed.
//    - Notice that for BiHeapifyEvenSinglePass() this appears to always
//       require only one pass while for BiHeapifyOddSinglePass() this
//       sometimes requires multiple passes.
//============================================================================

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "biheapify.h"
#include "biheapify_simple.h"
#include "biheapify_single_pass_success_rate.h"
#include "biheapify_time.h"

int main() {
  int start_total_num_nodes = 2;
  int increment_size        = 2; //Set this to 2 to test vectors that are only
                                 //odd or only even in size.
  //Will be applied to BiHeaps of sizes:
  // start_total_num_nodes + (muliple) * increment_size
  int end_total_num_nodes = static_cast<int>(1u << 15);
  int num_vecs_to_try = static_cast<int>(1u << 9);

  //For TimeBiHeapifies() only. For each std::vector that whose biheapification
  // is to be timed, repeat this process num_repititions_per_vec times.
  int num_repititions_per_vec = static_cast<int>(1u << 11);

  //For MeasureBiHeapifySuccessRate() only. Print information only after you've
  // gone through print_multiple new vector sizes.
  int print_multiple = 32;

  TimeBiHeapifies<int>(start_total_num_nodes, end_total_num_nodes,
                       num_vecs_to_try, num_repititions_per_vec, increment_size);
  MeasureBiHeapifySuccessRate<int>(start_total_num_nodes, end_total_num_nodes,
                                   num_vecs_to_try, increment_size,
                                   print_multiple);
  return 0;
}
