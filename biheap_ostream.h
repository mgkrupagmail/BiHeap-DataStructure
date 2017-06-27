/*
 * biheap_ostream.h
 *
 *  Created on: Jun 12, 2017
 *      Author: Matthew Gregory Krupa
 *
 * See biheapify.h for the definition of a biheap.
 * This file's class, biheap_ostream, outputs to any std::ostream a string that
 *  displays the first biheap_size_ nodes of a biheap graph of size biheap_size_.
 */

#ifndef BIHEAP_OSTREAM_H_
#define BIHEAP_OSTREAM_H_

#include <iomanip>
#include <iostream>
#include <iterator>

namespace biheap_ostream_ns {

#define BIHEAP_OSTREAM_DEFAULT std::cout

template<class iterator>
class biheap_ostream {

  typedef std::size_t size_type;

  iterator first_, one_past_last_; //The iterators defining the biheap.
  size_type biheap_size_;
  size_type num_nodes_in_heap_;
  size_type num_nodes_in_pure_heap_;
  size_type first_node_in_mirror_heap_;
  size_type num_nodes_in_last_row_of_pure_heap_;
  size_type max_level_of_even_pure_heap_;   //= number of rows in the pure heap.
  size_type max_width_of_elements_;
  size_type min_num_spaces_between_elements_ = 1;
  bool is_even_pure_heap_full_;

  std::ostream &ostrm_;
public:
  biheap_ostream(iterator first, iterator one_past_last, std::ostream &ostrm = BIHEAP_OSTREAM_DEFAULT) :
      ostrm_(ostrm), first_(first), one_past_last_(one_past_last) {
    ComputeMemberVariables();
  }
  ~biheap_ostream() {}

  static bool IsFullHeap(size_type num_nodes_in_heap) {
    return MaxLevelOfHeapOfGivenSize(num_nodes_in_heap) != MaxLevelOfHeapOfGivenSize(num_nodes_in_heap + 1);
  }
  static size_type MaxLevelOfHeapOfGivenSize(size_type num_nodes_in_heap) {
    if (num_nodes_in_heap <= 2)
      return num_nodes_in_heap;
    size_type max_level = 0;
    unsigned long long num_nodes = static_cast<unsigned long long>(num_nodes_in_heap);
    while (num_nodes != 0ull) {
      num_nodes = num_nodes >> 1;
      max_level++;
    }
    return max_level;
  }
  static size_type MaxLevelOfLargestFullSubheapOfHeapOfGivenSize(size_type num_nodes_in_heap) {
    auto max_level = MaxLevelOfHeapOfGivenSize(num_nodes_in_heap);
    if (IsFullHeap(num_nodes_in_heap))
      return max_level;
    else
      return max_level - 1;
  }
  static size_type NumNodesInFullHeapOfGivenNumOfLevels(size_type num_levels_of_full_heap) {
    if (num_levels_of_full_heap <= 2)
      return num_levels_of_full_heap;
    return static_cast<size_type>(1ull << (num_levels_of_full_heap - 1));
  }
  static size_type NumNodesInLargestFullHeapContainingHeapOfSize(size_type num_nodes_in_heap) {
    if(num_nodes_in_heap <= 1)
      return num_nodes_in_heap;
    auto num_nodes = num_nodes_in_heap;
    auto leading_one_counter = 0u;
    while (num_nodes != 1) {
      num_nodes = num_nodes >> 1;
      leading_one_counter++;
    }
    auto num_elements_in_full_heap_contained_in_given_heap = 1u;
    if (0x1u << (leading_one_counter) != 0)
      num_elements_in_full_heap_contained_in_given_heap = (0x1u << (leading_one_counter)) - 1;
    else
      num_elements_in_full_heap_contained_in_given_heap = -1;
    return num_elements_in_full_heap_contained_in_given_heap;
  }
  static size_type NumNodesInLastRowOfSmallestFullHeapContainingHeapOfSize(size_type num_nodes_in_heap) {
    if(num_nodes_in_heap <= 1)
      return num_nodes_in_heap;
    auto num_nodes = num_nodes_in_heap;
    auto leading_one_counter = 0u;
    while (num_nodes != 0) {
      num_nodes = num_nodes >> 1;
      leading_one_counter++;
    }
    return 0x1u << (leading_one_counter - 1);
  }
  static size_type NumNodesInLastRowOfLargestFullHeapContainedInHeapOfSize(size_type num_nodes_in_heap) {
    return NumNodesInLastRowHeapOfSize(NumNodesInFullHeapOfGivenNumOfLevels(MaxLevelOfLargestFullSubheapOfHeapOfGivenSize(num_nodes_in_heap)));
  }
  static size_type NumNodesInLastRowHeapOfSize(size_type num_nodes_in_heap) {
    if(num_nodes_in_heap <= 1)
      return num_nodes_in_heap;
    auto num_nodes = num_nodes_in_heap;
    auto leading_one_counter = 0u;
    while (num_nodes != 1) {
      num_nodes = num_nodes >> 1;
      leading_one_counter++;
    }
    auto num_elements_in_full_heap_contained_in_given_heap = 1u;
    if (0x1u << (leading_one_counter) != 0)
      num_elements_in_full_heap_contained_in_given_heap = (0x1u << (leading_one_counter)) - 1;
    else
      num_elements_in_full_heap_contained_in_given_heap = -1;
    return num_nodes_in_heap - num_elements_in_full_heap_contained_in_given_heap;
  }
  //level - The number of rows in the heap. The root is at level 0, the root's
  // two children at at level 1, etc.
  static size_type NumNodesInGivenLevelOfFullHeap(size_type level) {
    unsigned long long num_nodes_ull = 1ull << level;
    return static_cast<size_type>(num_nodes_ull);
  }

  void ComputeMemberVariables() {
    biheap_size_ = std::distance(first_, one_past_last_);
    auto i = biheap_size_ / 2;
    if (biheap_size_ <= 2)
      num_nodes_in_heap_ = biheap_size_;
    else if (biheap_size_ % 2 == 0)
      num_nodes_in_heap_ = 4*((i-(i%3))/3) + (i%3) + 1 - (((i+2)%3)/2);
    else {
      i--;
      num_nodes_in_heap_ = 4*((i-(i%3))/3) + (i%3);
    }
    num_nodes_in_pure_heap_ = (biheap_size_ / 2) + (biheap_size_ % 2);
    num_nodes_in_last_row_of_pure_heap_ = NumNodesInLastRowHeapOfSize(biheap_size_ / 2);
    first_node_in_mirror_heap_  = biheap_size_ - num_nodes_in_heap_ - (biheap_size_ % 2);

    max_level_of_even_pure_heap_ = MaxLevelOfHeapOfGivenSize(biheap_size_ / 2); //The "even" means that we use biheap_size_ / 2 regardless of whether biheap_size_ is even or odd/
    is_even_pure_heap_full_ = IsFullHeap(biheap_size_ / 2);
    max_width_of_elements_ = FindMaxWidthOfElements(first_, one_past_last_);
    if (min_num_spaces_between_elements_ < max_width_of_elements_)
      min_num_spaces_between_elements_ = max_width_of_elements_;
  }

  std::vector<char> ConstructHeapMiddleString(size_type biheap_size,
                                              size_type max_width_of_elements,
                                              size_type vec_size,
                                              iterator ele_it) const {
    auto char_vec = std::vector<char>(vec_size, ' ');
    auto row_str_it = char_vec.begin();
    auto num_chars_in_row = vec_size;
    row_str_it += (num_chars_in_row - max_width_of_elements) / 2;
    auto ele_str = std::to_string(*ele_it);
    std::copy(ele_str.begin(), ele_str.end(), row_str_it);
    return char_vec;
  }

  std::vector<std::vector<char>> ConstructHeap(size_type biheap_size,
                                               size_type num_spaces_in_between_elements_in_last_level,
                                               size_type max_width_of_elements) const {
    iterator ele_it = first_;
    auto strs = ConstructUpperHeap(biheap_size/2,
                                   num_spaces_in_between_elements_in_last_level,
                                   max_width_of_elements, ele_it);
    auto original_strs_size = strs.size();
    auto new_strs_size = 2 * strs.size() + 1;
    strs.resize(new_strs_size);


    if (biheap_size % 2 == 1) {
      /*auto row_num = original_strs_size;
      strs[row_num] = std::vector<char>(strs[0].size(), ' ');
      auto row_str_it = strs[row_num].begin();
      auto num_chars_in_row = strs[0].size();
      row_str_it += (num_chars_in_row - max_width_of_elements) / 2;
      auto ele_str = std::to_string(*ele_it);
      std::copy(ele_str.begin(), ele_str.end(), row_str_it);
      */
      strs[original_strs_size] = ConstructHeapMiddleString(biheap_size,
                                  max_width_of_elements, strs[0].size(), ele_it);
      ele_it++;
    }

    original_strs_size++;

    auto strs_lower = ConstructLowerHeap(biheap_size/2,
                                         num_spaces_in_between_elements_in_last_level,
                                         max_width_of_elements, ele_it);
    for (size_type i = 0; i < strs_lower.size(); i++) {
      strs[original_strs_size + i] = strs_lower[i];
    }
    return strs;
  }

  //Note: ele_it should usually be initialized to equal first_. After ConstructUpperHeap()
  // completes, ele_it will point to the first node in the pure max heap.
  std::vector<std::vector<char>> ConstructLowerHeap(size_type num_nodes_in_heap,
                                                    size_type num_spaces_in_between_elements_in_last_level,
                                                    size_type max_width_of_elements,
                                                    iterator &ele_it) const {
    if (num_spaces_in_between_elements_in_last_level <= 0) {
      num_spaces_in_between_elements_in_last_level = max_width_of_elements;
    }
    auto num_chars_in_row = GetRowSize(num_nodes_in_heap,
                              num_spaces_in_between_elements_in_last_level,
                              max_width_of_elements);
    std::vector<std::vector<char>> strs(max_level_of_even_pure_heap_, std::vector<char>(num_chars_in_row, ' '));
    size_type node_num = 0;
    auto element_str_width = max_width_of_elements + num_spaces_in_between_elements_in_last_level;
    for (size_type row_num = 0; row_num < static_cast<size_type>(strs.size()); row_num++) {
      auto row_str_it = strs[row_num].begin();
      auto i = GetNumLeadingSpaces(max_level_of_even_pure_heap_ - 1 - row_num,
                 num_chars_in_row,
                 num_spaces_in_between_elements_in_last_level,
                 max_width_of_elements, num_nodes_in_heap);
      row_str_it += i;
      auto num_spaces_in_between = GetNumSpacesInbetweenElementsInGivenRow(
                                    max_level_of_even_pure_heap_ - 1 - row_num,
                                    num_chars_in_row,
                                    num_spaces_in_between_elements_in_last_level,
                                    max_width_of_elements,
                                    num_nodes_in_heap);
      auto num_nodes_this_level = NumNodesInGivenLevelOfFullHeap(max_level_of_even_pure_heap_ - 1 - row_num);
      if (row_num == 0) {
        num_nodes_this_level = NumNodesInLastRowHeapOfSize(num_nodes_in_heap);
        auto num_nodes_in_last_level_of_full_heap = NumNodesInLastRowOfSmallestFullHeapContainingHeapOfSize(num_nodes_in_heap);
        row_str_it += (num_nodes_in_last_level_of_full_heap - num_nodes_this_level) *
            (element_str_width + num_spaces_in_between);
      }
      for (size_type i = 0; i < num_nodes_this_level; i++) {
        auto ele_str = std::to_string(*ele_it);
        std::copy(ele_str.begin(), ele_str.end(), row_str_it);
        ele_it++;
        row_str_it += element_str_width;
        if (i != num_nodes_this_level - 1)
          row_str_it += num_spaces_in_between;
      }
    }
    return strs;
  }

  //Note: ele_it should usually be initialized to equal first_. After ConstructUpperHeap()
  // completes, ele_it will point to the first node in the pure max heap.
  std::vector<std::vector<char>> ConstructUpperHeap(size_type num_nodes_in_heap,
                                                    size_type num_spaces_in_between_elements_in_last_level,
                                                    size_type max_width_of_elements,
                                                    iterator &ele_it) const {
    if (num_spaces_in_between_elements_in_last_level <= 0) {
      num_spaces_in_between_elements_in_last_level = max_width_of_elements;
    }
    auto num_chars_in_row = GetRowSize(num_nodes_in_heap,
                              num_spaces_in_between_elements_in_last_level,
                              max_width_of_elements);
    std::vector<std::vector<char>> strs(max_level_of_even_pure_heap_, std::vector<char>(num_chars_in_row, ' '));
    size_type node_num = 0;
    auto element_str_width = max_width_of_elements + num_spaces_in_between_elements_in_last_level;
    for (size_type row_num = 0; row_num < static_cast<size_type>(strs.size()); row_num++) {
      auto row_str_it = strs[row_num].begin();
      auto i = GetNumLeadingSpaces(row_num, num_chars_in_row,
                 num_spaces_in_between_elements_in_last_level,
                 max_width_of_elements, num_nodes_in_heap);
      row_str_it += i;
      auto num_spaces_in_between = GetNumSpacesInbetweenElementsInGivenRow(
                                    row_num, num_chars_in_row,
                                    num_spaces_in_between_elements_in_last_level,
                                    max_width_of_elements,
                                    num_nodes_in_heap);
      auto num_nodes_this_level = NumNodesInGivenLevelOfFullHeap(row_num);
      if (row_num == strs.size() - 1)
        num_nodes_this_level = NumNodesInLastRowHeapOfSize(num_nodes_in_heap);
      for (size_type i = 0; i < num_nodes_this_level; i++) {
        auto ele_str = std::to_string(*ele_it);
        std::copy(ele_str.begin(), ele_str.end(), row_str_it);
        ele_it++;
        row_str_it += element_str_width;
        if (i != num_nodes_this_level - 1)
          row_str_it += num_spaces_in_between;
      }
    }
    return strs;
  }

  size_type ElementWidth(iterator it) { return static_cast<size_type>(std::to_string(*it).length()); }

  size_type FindMaxWidthOfElements(iterator first, iterator one_past_last) {
    size_type max = 0;
    for (auto i = first; i != one_past_last; i++) {
      size_type width = ElementWidth(i);
      if (width > max)
        max = width;
    }
    return max;
  }

  static size_type GetNumLeadingSpaces(size_type row_num,
                                                  size_type num_chars_in_row,
                                                  size_type num_spaces_in_between_elements_in_last_level,
                                                  size_type max_width_of_elements,
                                                  size_type num_nodes_in_heap) {
    auto num_spaces_in_between_elements = GetNumSpacesInbetweenElementsInGivenRow(row_num,
                                            num_chars_in_row, num_spaces_in_between_elements_in_last_level,
                                            max_width_of_elements, num_nodes_in_heap);
    return num_spaces_in_between_elements / 2;
}

  static size_type GetNumSpacesInbetweenElementsInGivenRow(size_type row_num,
                                                  size_type num_chars_in_row,
                                                  size_type num_spaces_in_between_elements_in_last_level,
                                                  size_type max_width_of_elements,
                                                  size_type num_nodes_in_heap) {
    auto max_level = MaxLevelOfHeapOfGivenSize(num_nodes_in_heap);
    if (num_nodes_in_heap <= 1 || max_level == 0 || row_num == max_level || num_chars_in_row <= 0)
      return 1;

    auto num_nodes_in_level = NumNodesInGivenLevelOfFullHeap(row_num);
    if (num_nodes_in_level == 0)
      return 0;
    auto node_length = max_width_of_elements + num_spaces_in_between_elements_in_last_level;
    auto num_chars_used_by_nodes_in_level = num_nodes_in_level * node_length;
    size_type num_spaces = num_chars_in_row - num_chars_used_by_nodes_in_level;
    num_spaces /= num_nodes_in_level;
    return num_spaces;
  }

  static size_type GetRowSize(size_type num_nodes_in_heap,
                                         size_type num_spaces_in_between_elements,
                                         size_type max_width_of_elements) {
    size_type row_size = 0;
    auto num_nodes_in_last_row = NumNodesInLastRowOfSmallestFullHeapContainingHeapOfSize(num_nodes_in_heap);

    row_size += num_nodes_in_last_row
        * (max_width_of_elements + num_spaces_in_between_elements);
    return row_size;
  }

  void Print(bool print_extra_info = true, bool print_extra_new_line_at_end = true) {
    if (print_extra_info)
      ostrm_ << "num nodes in biheap = " << biheap_size_
             << " \tnum nodes in heap = " << num_nodes_in_heap_
             << " \tfirst node in mirror heap = " << first_node_in_mirror_heap_ << std::endl;
    if (biheap_size_ <= 3) {
      if (biheap_size_ >= 1) {
        ostrm_ << *first_ << '\n';
        if (biheap_size_ == 2)
          ostrm_ << '\n' << *(first_ + 1) << '\n';
        else
          ostrm_ << *(first_ + 1) << '\n' << *(first_ + 2) << '\n';
      }
      if (print_extra_new_line_at_end)
        ostrm_ << '\n';
      ostrm_.flush();
      return ;
    }

    auto max_width_of_elements = FindMaxWidthOfElements(first_, one_past_last_);
    auto strs = ConstructHeap(biheap_size_, 0, max_width_of_elements);
    PrintCharVectors(strs);
    if (print_extra_new_line_at_end)
      ostrm_ << std::endl;
    return ;
  }

  void PrintCharVectors(std::vector<std::vector<char>> strs) {
    for (size_type i = 0; i < static_cast<size_type>(strs.size()); i++) {
      auto &str = strs[i];
      for (size_type char_num = 0; char_num < static_cast<size_type>(str.size()); char_num++)
        ostrm_ << str[char_num];
      ostrm_ << '\n';
    }
    ostrm_.flush();
    return ;
  }

  void PrintPureMaxHeap() {
    if (biheap_size_ <= 2) {
      if (biheap_size_ == 1)
        ostrm_ << *first_ << '\n';
      return ;
    } else if (biheap_size_ <= 4) {
      ostrm_ << *first_ << '\n'
             << *(first_ + 1) << '\n';
      return ;
    }
    auto max_width_of_elements = FindMaxWidthOfElements(first_, one_past_last_);
    auto strs = ConstructUpperHeap(num_nodes_in_pure_heap_, 0, max_width_of_elements);
    PrintCharVectors(strs);
    if (biheap_size_ % 2 == 1) {
      strs[0] = ConstructHeapMiddleString(biheap_size_, max_width_of_elements,
                                          strs[0].size(),
                                          first_ + biheap_size_/2);
      strs.resize(1);
      PrintCharVectors(strs);
    }
    return ;
  }

  void PrintPureMinHeap() {
    if (biheap_size_ == 0)
      return ;
    if (biheap_size_ == 1) {
      ostrm_ << *first_ << '\n';
      return ;
    } else if (biheap_size_ <= 3) {
      ostrm_ << *(first_ + 1) << '\n';
      if (biheap_size_ == 3)
        ostrm_ << *(first_ + 2) << '\n';
      return ;
    } else if (biheap_size_ == 4) {
      ostrm_ << *(first_ + 2) << '\n'
             << *(first_ + 3) << '\n';
      return ;
    }
    auto max_width_of_elements = FindMaxWidthOfElements(first_, one_past_last_);
    auto strs = ConstructLowerHeap(num_nodes_in_pure_heap_, 0, max_width_of_elements);
    if (biheap_size_ % 2 == 1) {
      std::vector<std::vector<char>> str_middle(1);
      str_middle[0] = ConstructHeapMiddleString(biheap_size_, max_width_of_elements,
                                                strs[0].size(),
                                                first_ + biheap_size_/2);
      PrintCharVectors(str_middle);
    }
    PrintCharVectors(strs);
    return ;
  }
};

} //biheap_ostream_ns

namespace {
typedef std::size_t size_type;
}

template<class iterator>
void PrintBiHeap(iterator first, size_type total_num_nodes, std::ostream &ostrm = BIHEAP_OSTREAM_DEFAULT) {
  biheap_ostream_ns::biheap_ostream<iterator> bho(first, first + total_num_nodes, ostrm);
  bho.Print();
  return ;
}

#undef BIHEAP_OSTREAM_DEFAULT

#endif /* BIHEAP_OSTREAM_H_ */
