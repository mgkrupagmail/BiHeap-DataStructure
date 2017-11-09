/*
 * biheap_sift_test_correctness.h
 *
 *  Created on: Jul 25, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAP_SIFT_TEST_CORRECTNESS_H_
#define BIHEAP_SIFT_TEST_CORRECTNESS_H_

#include <random>
#include "random_helpers.h"
#include "biheap_sift.h"

template<class T, typename size_type = std::size_t>
bool BiHeapSiftTestCorrectness(long start_total_num_nodes,
                               long end_total_num_nodes,
                               long num_vecs_to_try = 1,
                               long increment_size = 1) {
  int total_num_nodes = start_total_num_nodes;
  std::random_device rnd_device;
  std::mt19937 generator(rnd_device());
  while (total_num_nodes <= end_total_num_nodes) {
    std::cout << "Sifting in BiHeaps of size " << total_num_nodes << " \t";
    std::cout.flush();
    std::vector<T> vec(total_num_nodes + 1);
    std::uniform_int_distribution<T> dist(0, total_num_nodes - 1);
    for (auto vec_counter = 0l; vec_counter < num_vecs_to_try; vec_counter++) {
      randomhelpers::FillWithRandomNumbers(vec.begin(), vec.end(),
                0, static_cast<T>(4*total_num_nodes));//std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      size_type pos_hc;
      pos_hc = dist(generator);
      //Form a BiHeap.
      BiHeapify<typename std::vector<T>::iterator, size_type>(vec.begin(), total_num_nodes);
      //Change one randomly selected element to some random value (which was
      // stored outside of the biheap in vec[total_num_nodes]).
      vec[pos_hc] = vec[total_num_nodes];
      //Sift the changed element.
      BiHeapSift<typename std::vector<T>::iterator, size_type>(vec.begin(), total_num_nodes, pos_hc);
      //Check that we once again have a biheap.
      bool is_biheap = IsBiHeap<typename std::vector<T>::iterator, size_type>(vec.begin(), total_num_nodes);
      if (!is_biheap) {
        std::cout << "Sift() failed to form a biheap. "
                  << "total_num_nodes = " << total_num_nodes << std::endl;
        return false;
      }
    }
    std::cout << "Successfully Sifted in BiHeaps of size " << total_num_nodes << std::endl;
    total_num_nodes += increment_size;
  }
  return true;
}

#endif /* BIHEAP_SIFT_TEST_CORRECTNESS_H_ */
