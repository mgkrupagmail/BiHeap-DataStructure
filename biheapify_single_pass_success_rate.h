/*
 * biheapify_single_pass_success_rate.h
 *
 *  Created on: Jun 27, 2017
 *      Author: Matthew Gregory Krupa
 *
 *  MeasureBiHeapifySuccessRate(), which is used to measure how many times
//   a BiHeapify pass (i.e. one of BiHeapifyEvenSinglePass(),
//   BiHeapifyOddSinglePass(), or BiHeapifySimpleSinglePass()) must be
//   called before a BiHeap is formed.
// - Notice that for BiHeapifyEvenSinglePass() this appears to always
//    require only one pass while for BiHeapifyOddSinglePass() this sometimes
//    requires multiple passes.
 */

#ifndef BIHEAPIFY_SINGLE_PASS_SUCCESS_RATE_H_
#define BIHEAPIFY_SINGLE_PASS_SUCCESS_RATE_H_

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "biheap_common.h"
#include "biheapify.h"
#include "biheapify_simple.h"
#include "random_helpers.h"

template<class T>
std::string GetDesciption(std::vector<T> fail_counter,
                          std::vector<T> call_counter) {
  int max_num_calls = 1; //The number of non-zero entries in call_counter;
  for ( ; max_num_calls < static_cast<int>(call_counter.size())
          && call_counter[max_num_calls] != 0; max_num_calls++) ;

  fail_counter.resize(max_num_calls);
  call_counter.resize(max_num_calls);
  std::vector<T> success_counter(max_num_calls, static_cast<T>(0));
  for (int i = 0; i < max_num_calls; i++)
    success_counter[i] = call_counter[i] - fail_counter[i];

  //success_by_num_call[i] = success_counter[1] + ... + success_counter[i]
  std::vector<T> sum_success_by_call_num(max_num_calls, static_cast<T>(0));
  std::vector<T> sum_calls_by_call_num(max_num_calls, static_cast<T>(0));
  sum_success_by_call_num[1] = success_counter[1];
  sum_calls_by_call_num[1]   = call_counter[1];
  for (int i = 1; i < max_num_calls; i++) {
    sum_success_by_call_num[i] = sum_success_by_call_num[i - 1] + success_counter[i];
    sum_calls_by_call_num[i]   = sum_calls_by_call_num[i - 1] + call_counter[i];
  }
  std::stringstream strstrm;
  auto max_width_success_counter = 1u;
  auto max_width_call_counter = 1u;
  auto max_width_sum_success_by_call_num = 1u;
  auto max_width_sum_calls_by_call_num = 1u;

  for (int i = 1; i < max_num_calls; i++) {
    if (max_width_success_counter < std::to_string(success_counter[i]).length())
      max_width_success_counter = std::to_string(success_counter[i]).length();
    if (max_width_call_counter < std::to_string(call_counter[i]).length())
      max_width_call_counter = std::to_string(call_counter[i]).length();
    if (max_width_success_counter < std::to_string(sum_success_by_call_num[i]).length())
      max_width_sum_success_by_call_num = std::to_string(sum_success_by_call_num[i]).length();
    if (max_width_call_counter < std::to_string(sum_calls_by_call_num[i]).length())
      max_width_sum_calls_by_call_num = std::to_string(sum_calls_by_call_num[i]).length();

    strstrm << std::string("call ")        << i
            << std::string(" \tsuccess[")  << i
            << std::string("]/num_calls[") << i << std::string("] = ");
    strstrm << std::setw(max_width_success_counter) << success_counter[i]
            << std::string(" / ")
            << std::setw(max_width_call_counter)    << call_counter[i];
    if (call_counter[i] != 0) {
      strstrm << std::string(" = \t")
              << std::to_string(static_cast<long double>(success_counter[i])
                                / static_cast<long double>(call_counter[i]));
    }
    strstrm << std::string("  \tsum_success_leq[") << std::to_string(i)
            << std::string("]/sum_num_calls_leq[") << std::to_string(i)
            << std::string("] = ");
    strstrm << std::setw(max_width_sum_success_by_call_num)
            << std::to_string(sum_success_by_call_num[i])
            << std::string(" / ")
            << std::setw(max_width_sum_calls_by_call_num)
            << std::to_string(sum_calls_by_call_num[i]);
    if (sum_calls_by_call_num[i] != 0) {
      strstrm << std::string(" = \t")
         << std::to_string(static_cast<long double>(sum_success_by_call_num[i])
                         / static_cast<long double>(sum_calls_by_call_num[i]));
    }
    strstrm << std::string("\n");
  }

  auto total_success = sum_success_by_call_num[sum_success_by_call_num.size() - 1];
  auto total_calls   = sum_calls_by_call_num[sum_calls_by_call_num.size() - 1];

  strstrm << std::string("total_success/total_calls \t= ")
          << std::to_string(total_success)
          << std::string(" / ") << std::to_string(total_calls);
  if (total_calls != 0) {
      strstrm << std::string(" = \t")
              << std::to_string(static_cast<long double>(total_success)
                                / static_cast<long double>(total_calls));
  }
  return strstrm.str();
}

/* By default this function measures the success rate of BiHeapifySinglePass().
 * To measure the success rate of BiHeapifyEvenSinglePass() then comment
 *  out the two instances of BiHeapifySinglePass() and un-comment out the two
 *  instances of BiHeapifyEvenSinglePass().
 *  Ditto for BiHeapifySimpleSinglePass() and BiHeapifyOddSinglePass().
 */
template<class T> void MeasureBiHeapifySuccessRate(int start_total_num_nodes,
                                                   int end_total_num_nodes,
                                                   int num_vecs_to_try = 1,
                                                   int increment_size = 1,
                                                   int print_multiple = 32) {
  std::vector<long long> fail_counter(100, 0l);
  std::vector<long long> try_counter(100, 0l);
  long long total_tries = 0;
  int total_num_nodes = start_total_num_nodes;
  while (total_num_nodes <= end_total_num_nodes) {
    auto initial_size = total_num_nodes;
    for (int i = 0; i < print_multiple; i++, total_num_nodes += increment_size) {
      if (total_num_nodes > end_total_num_nodes)
        break ;
      for (int vec_counter = 0; vec_counter < num_vecs_to_try; vec_counter++) {
        int try_num = 1;
        std::vector<T> vec(total_num_nodes);
        randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                  std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

        //BiHeapifySimpleSinglePass(vec.begin(), vec.size());
        //BiHeapifyEvenSinglePass(vec.begin(), vec.size());
        //BiHeapifyOddSinglePass(vec.begin(), vec.size());
        BiHeapifySinglePass(vec.begin(), vec.size());

        total_tries++;
        try_counter[try_num]++;
        while (IsBiheap(vec.begin(), vec.size(), 0, 0, false) == false) {
          fail_counter[try_num]++;
          if (try_num >= 10) {
            std::cout << "Tried and failed to biheapify vector " << try_num
                      << " times. Quitting this vector." << std::endl;
            break;
          }
          try_num++;

          //BiHeapifySimpleSinglePass(vec.begin(), vec.size());
          //BiHeapifyEvenSinglePass(vec.begin(), vec.size());
          //BiHeapifyOddSinglePass(vec.begin(), vec.size());
          BiHeapifySinglePass(vec.begin(), vec.size());

          if (static_cast<unsigned int>(try_num) >= try_counter.size()) {
            fail_counter.resize(2 * try_num);
            try_counter.resize(2 * try_num);
          }
          try_counter[try_num]++;
          total_tries++;
        }
        if (IsBiheap(vec.begin(), vec.size(), 0, 0, false) == false) {
          std::cout << "Failed to BiHeapify() vector of size "
                    << vec.size() << std::endl;
          if (total_num_nodes < (1 << 10))
            PrintBiHeap(vec.begin(), i);
          return ;
        }
      }
    }
    std::cout << "BiHeap start size = " << initial_size << std::endl;
    std::cout << GetDesciption(fail_counter, try_counter) << std::endl;
    std::cout << std::endl;
  }
  return ;
}



#endif /* BIHEAPIFY_SINGLE_PASS_SUCCESS_RATE_H_ */
