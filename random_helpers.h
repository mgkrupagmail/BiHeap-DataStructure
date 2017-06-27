/*
 * random_helpers.h
 *
 *  Created on: May 5, 2017
 *      Author: Matthew Gregory Krupa
 */

#ifndef RANDOM_HELPERS_H_
#define RANDOM_HELPERS_H_

#include <iostream>
#include <random>
#include <vector>

namespace randomhelpers {
#define FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE false
#define FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM   std::cout
#define FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE "Finished filling vector with random goodness.\n"

template<typename T>
static inline void FillVectorWithRandomNumbersPrintStartMessage(const std::vector<T> &vec,
                const T a, const T b, const bool verbose, std::ostream &ostrm) {
  if (verbose) {
    ostrm << "Started filling vector of size " << vec.size()
          << " with random numbers from the range [" << a << ", " << b << ")." << '\n';
    ostrm.flush();
  }
}

template<typename T>
void FillVectorWithRandomNumbers(std::vector<T> &vec, T a, T b,
                                 bool verbose = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE,
                                 std::ostream &ostrm = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM) {
  std::random_device rnd_device;
  std::mt19937 generator(rnd_device());
  std::uniform_int_distribution<T> dist(a, b);
  FillVectorWithRandomNumbersPrintStartMessage(vec, a, b, verbose, ostrm);
  for (unsigned int i = 0; i < vec.size(); i++)
    vec[i] = dist(generator);
  if (verbose) {
    ostrm << FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE;
    ostrm.flush();
  }
}

void FillVectorWithRandomNumbers(std::vector<float> &vec, float a, float b,
                                 bool verbose = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE,
                                 std::ostream &ostrm = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<float> dis(a, b);
  FillVectorWithRandomNumbersPrintStartMessage(vec, a, b, verbose, ostrm);
  for (unsigned int i = 0; i < vec.size(); i++)
    vec[i] = dis(generator);
  if (verbose) {
    ostrm << FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE;
    ostrm.flush();
  }
}

void FillVectorWithRandomNumbers(std::vector<double> &vec, double a, double b,
                                 bool verbose = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE,
                                 std::ostream &ostrm = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<double> dis(a, b);
  FillVectorWithRandomNumbersPrintStartMessage(vec, a, b, verbose, ostrm);
  for (unsigned int i = 0; i < vec.size(); i++)
    vec[i] = dis(generator);
  if (verbose) {
    ostrm << FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE;
    ostrm.flush();
  }
}

void FillVectorWithRandomNumbers(std::vector<long double> &vec, long double a, long double b,
                                 bool verbose = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE,
                                 std::ostream &ostrm = FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM) {
  std::random_device rnd_device;
  std::mt19937_64 generator(rnd_device());
  std::uniform_real_distribution<long double>  dis(a, b);
  FillVectorWithRandomNumbersPrintStartMessage(vec, a, b, verbose, ostrm);
  for (unsigned int i = 0; i < vec.size(); i++)
    vec[i] = dis(generator);
  if (verbose) {
    ostrm << FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE;
    ostrm.flush();
  }

#undef FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_VERBOSE
#undef FILL_VECTOR_WITH_RANDOM_NUMBERS_DEFAULT_OSTRM
#undef FILL_VECTOR_WITH_RANDOM_NUMBERS_FINISHED_MESSAGE
}

} //END NAMESPACE randomhelpers

#endif /* RANDOM_HELPERS_H_ */
