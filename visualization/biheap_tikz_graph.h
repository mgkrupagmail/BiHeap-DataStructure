/*
 * biheap_tikz_graph.h
 *
 * The GetTikzGraph() function creates the LaTeX tikz
 *  code to display a BiHeapGraph. The nodes of the graph
 *  may either have values assigned to them or else their
 *  text will be the coordinates of that node.
 *
 *  Created on: Sep 26, 2017
 *      Author: Matthew Gregory Krupa
 *   Copyright: Matthew Gregory Krupa
 */

#ifndef BIHEAP_TIKZ_GRAPH_H_
#define BIHEAP_TIKZ_GRAPH_H_

#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../biheapify.h"

#ifndef FLIP
#define FLIP(a) ((total_num_nodes - 1) - (a))
#endif

template<typename size_type = std::size_t>
size_type GetNumNodesInLastRow(size_type total_num_nodes) {
  if (total_num_nodes <= 1)
    return total_num_nodes;
  size_type num_nodes_so_far = 1;
  size_type num_nodes_this_row = 1;
  while (num_nodes_so_far < total_num_nodes) {
    num_nodes_this_row *= 2;
    num_nodes_so_far += num_nodes_this_row;
  }
  return num_nodes_this_row;
}

template<typename size_type = std::size_t>
std::string GetNodeNameHc(size_type pos_hc, size_type total_num_nodes, bool with_parentheses = true, bool pad_with_spaces_to_max_length = false) {
  std::stringstream strm;
  if (with_parentheses)
    strm << "(";
  strm << "s" << total_num_nodes << "n" << pos_hc;
  if (with_parentheses)
    strm << ")";
  if (pad_with_spaces_to_max_length) {
    size_type desired_node_str_width = std::to_string(total_num_nodes).length() + std::to_string(total_num_nodes - 1).length() + 2;
    if (with_parentheses)
      desired_node_str_width += 2; //Take into account the opening and closing parentheses.
    while (static_cast<size_type>(strm.str().length()) < desired_node_str_width)
      strm << ' ';
  }
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetNodeNameMc(size_type pos_hc, size_type total_num_nodes, bool with_parentheses = true) {
  return GetNodeNameHc<size_type>(FLIP(pos_hc), total_num_nodes, with_parentheses);
}


template<typename size_type = std::size_t>
std::string GetNodeTextHc(size_type pos_hc, size_type total_num_nodes, bool with_braces = true) {
  std::stringstream strm;
  if (with_braces)
    strm << "{";
  strm << "$\\FIGCO{" << pos_hc << "}{" << (FLIP(pos_hc))<< "}$";
  if (with_braces)
    strm << "}";
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetNodeTextMc(size_type pos_hc, size_type total_num_nodes, bool with_braces = true) {
  return GetNodeNameHc<size_type>(FLIP(pos_hc), total_num_nodes, with_braces);
}

template<typename size_type = std::size_t>
bool IsInLeftBranch(size_type pos_hc) {
  assert(pos_hc != 0);
  while (pos_hc > 2)
    pos_hc = Parent<size_type>(pos_hc);
  if (pos_hc == 2)
    return false;
  else
    return true;
}


template<typename size_type = std::size_t>
bool IsInLeftBranch(size_type pos_hc, size_type subtree_start_hc) {
  assert(pos_hc != 0);
  size_type subtree_left_child = LeftChild<size_type>(subtree_start_hc);
  auto subtree_right_child = subtree_left_child + 1;
  while (pos_hc > subtree_right_child)
    pos_hc = Parent<size_type>(pos_hc);
  if (pos_hc == subtree_right_child)
    return false;
  else if (pos_hc == subtree_left_child)
    return true;
  assert(false);
  return false;
}

//Node 0 has depth 0.
template<typename size_type = std::size_t>
size_type GetDepth(size_type pos) {
  if (pos == 0)
    return 0;
  auto i = 0;
  while (pos != 0) {
    i++;
    pos = Parent<size_type>(pos);
  }
  return i;
}

//total_num_nodes NOT needed for this. Imagine that total_num_nodes == pos_hc.
//y_scale_down_ratio < 1.0 is meant to "spread" adjacent nodes out apart more,
// while y_scale_down_ratio > 1.0 brings then closer together.
//y_scale_down_ratio == 1.0 places nodes directly in the middle
template<typename size_type = std::size_t>
long double GetYCoordinateRec(size_type pos, size_type subtree_start, size_type total_num_nodes,
              long double start_y, long double end_y, long double y_scale_down_ratio = 1.0) {
  auto distance = end_y - start_y;
  auto half_distance = distance / 2;
  auto mid_point = start_y + half_distance;
  if (pos == subtree_start)
    return mid_point;

  auto new_start_y = start_y;
  auto new_end_y = end_y;
  auto new_subtree_start = subtree_start;

  bool is_in_left_branch = IsInLeftBranch<size_type>(pos, subtree_start);
  if (is_in_left_branch) {
    new_end_y = start_y + y_scale_down_ratio * half_distance;
    new_subtree_start = LeftChild<size_type>(subtree_start);
  } else {
    new_start_y = end_y - y_scale_down_ratio * half_distance;
    new_subtree_start = RightChild<size_type>(subtree_start);
  }

  return GetYCoordinateRec<size_type>(pos, new_subtree_start, total_num_nodes,
                                     new_start_y, new_end_y, y_scale_down_ratio);
}

//total_num_nodes NOT needed for this. Imagine that total_num_nodes == pos_hc.
//y_scale_down_ratio < 1.0 is meant to "spread" adjacent nodes out apart more,
// while y_scale_down_ratio > 1.0 brings then closer together.
//y_scale_down_ratio == 1.0 places nodes directly in the middle
template<typename size_type = std::size_t>
long double GetYCoordinate(size_type pos, size_type total_num_nodes, long double start_y,
                       long double end_y, long double y_scale_down_ratio = 1.0) {
  return GetYCoordinateRec<size_type>(pos, 0, total_num_nodes, start_y, end_y, y_scale_down_ratio);
}

//Should have 0 <= x_middle_separation_distance_ratio < 1.0
template<typename size_type = std::size_t>
long double GetXCoordinate(size_type pos, size_type total_num_nodes, long double start_x, long double end_x,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l) {
  auto distance = end_x - start_x;
  auto half_distance = distance / 2;
  auto mid_point = start_x + half_distance;

  if (total_num_nodes % 2 == 1 && pos == total_num_nodes / 2)
    return mid_point;

  auto half_distance_scaled = half_distance * (1.0l - x_middle_separation_distance_ratio);

  size_type size_of_pure_heap_ignoring_middle = (total_num_nodes / 2);// - (total_num_nodes % 2);
  size_type max_depth_of_pure_heap_ignoring_middle = GetDepth<size_type>(size_of_pure_heap_ignoring_middle - 1);
  size_type depth;
  bool is_in_pure_min_heap = pos < total_num_nodes / 2;
  if (is_in_pure_min_heap)
    depth = GetDepth<size_type>(pos);
  else
    depth = GetDepth<size_type>(FLIP(pos));

  long double offset_from_end = static_cast<long double>(depth) * (half_distance_scaled / (max_depth_of_pure_heap_ignoring_middle + x_num_imaginary_middles_nodes_on_one_side));

  if (is_in_pure_min_heap)
    return start_x + offset_from_end;
  else
    return end_x - offset_from_end;
}

template<typename size_type = std::size_t>
void GetCoordinates(size_type total_num_nodes, long double start_x, long double end_x,
    long double start_y, long double end_y,
    std::vector<long double> &x_coordinates,
    std::vector<long double> &y_coordinates,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    long double y_scale_down_ratio = 1.0) {
  assert(start_x <= end_x && start_y <= end_y && total_num_nodes >= 1);
  if (total_num_nodes == 1) {
    x_coordinates[0] = (start_x + end_x) / 2;
    y_coordinates[0] = (start_y + end_y) / 2;
    return ;
  }
  for (size_type pos_hc = 0; pos_hc < total_num_nodes; pos_hc++) {
    auto x_coordinate = GetXCoordinate<size_type>(pos_hc, total_num_nodes, start_x, end_x, x_middle_separation_distance_ratio, x_num_imaginary_middles_nodes_on_one_side);
    long double y_coordinate;
    if (total_num_nodes % 2 == 1 && pos_hc == total_num_nodes / 2) {
      y_coordinate = (start_y + end_y) / 2;
    } else if (pos_hc < total_num_nodes / 2) {
      y_coordinate = GetYCoordinate<size_type>(pos_hc, total_num_nodes, start_y, end_y, y_scale_down_ratio);
    } else {
      y_coordinate = GetYCoordinate<size_type>(FLIP(pos_hc), total_num_nodes, start_y, end_y, y_scale_down_ratio);
      auto offset_from_bottom = y_coordinate - start_y;
      y_coordinate = end_y - offset_from_bottom;
    }
    x_coordinates[pos_hc] = x_coordinate;
    y_coordinates[pos_hc] = y_coordinate;
  }
  return ;
}

template<typename size_type = std::size_t>
long double GetDefaultNodeXCoordinateLength(size_type total_num_nodes, long double start_x, long double end_x,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    long double num_invisible_columns_of_nodes_that_can_fit_between_any_two_actual_columns = 1.0l) {
  auto distance = end_x - start_x;
  auto half_distance = distance / 2;
  auto half_distance_scaled = half_distance * (1.0l - x_middle_separation_distance_ratio);

  size_type size_of_pure_heap_ignoring_middle = (total_num_nodes / 2);//NOT: total_num_nodes/2 - (total_num_nodes % 2);
  if (size_of_pure_heap_ignoring_middle <= 1)
    return half_distance_scaled / 2;
  size_type max_depth_of_pure_heap_ignoring_middle = GetDepth<size_type>(size_of_pure_heap_ignoring_middle - 1);
  //Pretend that in between each two columns there can fit exactly num_invisible_columns_of_nodes_that_can_fit_between_any_two_actual_columns invisible columns of nodes.
  long double total_num_real_and_invisible_columns_of_nodes_per_actual_column = 1.0l + num_invisible_columns_of_nodes_that_can_fit_between_any_two_actual_columns;
  long double total_num_real_and_invisible_columns_of_nodes_ = total_num_real_and_invisible_columns_of_nodes_per_actual_column * static_cast<long double>(max_depth_of_pure_heap_ignoring_middle) + 1.0l;
  return half_distance_scaled / total_num_real_and_invisible_columns_of_nodes_;
}

template<typename size_type = std::size_t>
long double GetDefaultNodeYCoordinateLength(size_type total_num_nodes, long double start_y, long double end_y,
    long double num_invisible_nodes_that_can_fit_between_any_two_actual_nodes_in_last_column = 1.0l) {
  auto distance = end_y - start_y;
  auto half_distance = distance / 2;

  size_type size_of_pure_heap_ignoring_middle = (total_num_nodes / 2);//NOT: total_num_nodes/2 - (total_num_nodes % 2);
  size_type num_nodes_in_last_column_of_one_side = GetNumNodesInLastRow(size_of_pure_heap_ignoring_middle);
  if (size_of_pure_heap_ignoring_middle <= 1)
    return half_distance;
  //Pretend that in between each two node there can fit exactly num_invisible_nodes_that_can_fit_between_any_two_actual_nodes_in_last_column invisible nodes.
  long double total_num_real_and_invisible_nodes_per_actual_node_in_last_column = 1.0l + num_invisible_nodes_that_can_fit_between_any_two_actual_nodes_in_last_column;
  long double total_num_real_and_invisible_nodes_ = num_nodes_in_last_column_of_one_side * static_cast<long double>(num_invisible_nodes_that_can_fit_between_any_two_actual_nodes_in_last_column)  + 1.0l;
  return distance / total_num_real_and_invisible_nodes_;
}

//Assumes that the parent and child are rectangles of size node_x_length by node_y_length
// and that that their x and y coordinates given in x_coordinates and y_coordinates
// are the centers of these rectangles.
//These rectangles are not allowed to be made to intersect and they must also be at least
// min_x_distance_between_nodes and min_y_distance_between_nodes away from each other
// in each of these directions.
//In addition, after execution we will have that the distance between the x-coordinates
// of the parent and child is >= min_x_child_offset_from_parent.
//Only the coordinates of the child node are potentially changed. The coordinates of the parent node are not changed.
template<typename size_type = std::size_t>
void MoveChildCloseToParent(size_type total_num_nodes,
    std::vector<long double> &x_coordinates,
    std::vector<long double> &y_coordinates,
    size_type parent_hc,
    size_type child_hc,
    long double min_x_distance_between_nodes,
    long double min_y_distance_between_nodes,
    long double node_x_length,
    long double node_y_length,
    long double min_x_child_offset_from_parent) {
  size_type num_nodes_on_left = total_num_nodes / 2;
  assert((parent_hc < num_nodes_on_left && child_hc < num_nodes_on_left) || (FLIP(parent_hc) < num_nodes_on_left && FLIP(child_hc) < num_nodes_on_left));
  bool are_nodes_on_left = child_hc < num_nodes_on_left;
  assert((are_nodes_on_left && Parent(child_hc) == parent_hc) || (!are_nodes_on_left && Parent(FLIP(child_hc)) == FLIP(parent_hc)));
  if (are_nodes_on_left && x_coordinates[parent_hc] + node_x_length / 2 + min_x_distance_between_nodes < x_coordinates[child_hc] - node_x_length / 2)
    x_coordinates[child_hc] = x_coordinates[parent_hc] + node_x_length + min_x_distance_between_nodes;
  if (!are_nodes_on_left && x_coordinates[parent_hc] - node_x_length / 2 - min_x_distance_between_nodes > x_coordinates[child_hc] + node_x_length / 2)
    x_coordinates[child_hc]  = x_coordinates[parent_hc] - node_x_length - min_x_distance_between_nodes;
  if (y_coordinates[child_hc] == y_coordinates[parent_hc])
    return ;

  auto top_y_of_child     = y_coordinates[child_hc]  + node_y_length / 2;
  auto bottom_y_of_child  = y_coordinates[child_hc]  - node_y_length / 2;
  auto top_y_of_parent    = y_coordinates[parent_hc] + node_y_length / 2;
  auto bottom_y_of_parent = y_coordinates[parent_hc] - node_y_length / 2;
  //If moving the child furthure would cause the nodes to intersect then return.
  if (y_coordinates[child_hc] < y_coordinates[parent_hc] && top_y_of_child + min_y_distance_between_nodes > bottom_y_of_parent) {
    return ;
  }
  if (y_coordinates[child_hc] > y_coordinates[parent_hc] && bottom_y_of_child - min_y_distance_between_nodes < top_y_of_parent) {
    return ;
  }
  //At this point, the distance between the extreme y-coordinates of the parent and child are >= min_y_distance_between_nodes.
  if (are_nodes_on_left) {
    x_coordinates[child_hc] = x_coordinates[parent_hc] + min_x_child_offset_from_parent;
  } else {
    x_coordinates[child_hc] = x_coordinates[parent_hc] - min_x_child_offset_from_parent;
  }
  return ;
}

template<typename size_type = std::size_t>
void TightenCoordinates(size_type total_num_nodes, long double start_x, long double end_x,
    long double start_y, long double end_y,
    std::vector<long double> &x_coordinates,
    std::vector<long double> &y_coordinates,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    long double y_scale_down_ratio = 1.0,
    long double min_x_distance_between_nodes_as_ratio_of_node_x_length = 1.5,
    long double min_y_distance_between_nodes_as_ratio_of_node_y_length = 0.5,
    long double node_x_length = 0.0l,
    long double node_y_length = 0.0l,
    long double min_x_child_offset_from_parent = -1.0l) {
  assert(start_x <= end_x && start_y <= end_y && total_num_nodes >= 1);
  if (total_num_nodes <= 1)
    return ;
  if (node_x_length <= 0.0l)
    node_x_length = GetDefaultNodeXCoordinateLength(total_num_nodes, start_x, end_x, x_middle_separation_distance_ratio, x_num_imaginary_middles_nodes_on_one_side);
  if (node_y_length <= 0.0l)
    node_y_length = GetDefaultNodeXCoordinateLength(total_num_nodes, start_y, end_y);
  if (min_x_child_offset_from_parent < 0.0l)
    min_x_child_offset_from_parent = node_x_length * (3.0l / 4.0l);

  long double min_x_distance_between_nodes = min_x_distance_between_nodes_as_ratio_of_node_x_length * node_x_length;
  long double min_y_distance_between_nodes = min_y_distance_between_nodes_as_ratio_of_node_y_length * node_y_length;

  size_type last_node_to_check = total_num_nodes / 2 - 1;
  size_type parent_of_last_node_to_check = Parent<size_type>(last_node_to_check);
  for (size_type pos_hc = 0; pos_hc <= parent_of_last_node_to_check; pos_hc++) {
    size_type left_child_hc  = LeftChild<size_type>(pos_hc);
    size_type right_child_hc = RightChild<size_type>(pos_hc);

    size_type child_hc = left_child_hc;
    MoveChildCloseToParent(total_num_nodes, x_coordinates, y_coordinates, pos_hc, child_hc,
        min_x_distance_between_nodes, min_y_distance_between_nodes, node_x_length, node_y_length, min_x_child_offset_from_parent);
    MoveChildCloseToParent(total_num_nodes, x_coordinates, y_coordinates, FLIP(pos_hc), FLIP(child_hc),
        min_x_distance_between_nodes, min_y_distance_between_nodes, node_x_length, node_y_length, min_x_child_offset_from_parent);

    if (right_child_hc > last_node_to_check)
      break ;
    child_hc = right_child_hc;
    MoveChildCloseToParent(total_num_nodes, x_coordinates, y_coordinates, pos_hc, child_hc,
        min_x_distance_between_nodes, min_y_distance_between_nodes, node_x_length, node_y_length, min_x_child_offset_from_parent);
    MoveChildCloseToParent(total_num_nodes, x_coordinates, y_coordinates, FLIP(pos_hc), FLIP(child_hc),
        min_x_distance_between_nodes, min_y_distance_between_nodes, node_x_length, node_y_length, min_x_child_offset_from_parent);
  }
  return ;
}

template<typename size_type = std::size_t>
std::string GetDefineNodeLineHc(size_type pos_hc, size_type total_num_nodes, long double x, long double y, std::string node_text = std::string()) {
  std::stringstream strm;
  int node_width = std::to_string(total_num_nodes).length() + std::to_string(total_num_nodes - 1).length() + 4;
  const int unsigned_precision = 5; //5 decimal digits should be printed. Note that std::fixed makes sure that exactly 5 decimal digits are outputted.
  const int unsigned_coor_width = unsigned_precision + 2; //+2 since it's expected that (int)x and (int)y both have at most 2 digits.
  int coor_width;
  int precision = unsigned_precision;
  strm << "\\node " << std::left << std::setw(node_width) << GetNodeNameHc<size_type>(pos_hc, total_num_nodes, true);
  coor_width = unsigned_coor_width;
  if (x >= 0.0l)
    coor_width++; //To have the decimal points line up in case some x (resp. y) coordinates are positive and others are negative.
  strm << "  at  (" << std::setw(coor_width) << std::fixed << std::setprecision(precision) << x;
  coor_width = unsigned_coor_width;
  if (y >= 0.0l)
    coor_width++;
  strm << ", "      << std::setw(coor_width) << std::fixed << std::setprecision(precision) << y << ")  ";
  strm << "{";
  if (node_text.empty())
    strm << GetNodeTextHc<size_type>(pos_hc, total_num_nodes, false);
  else
    strm << node_text;
  strm << "};\n";

  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetDefineNodeLine(size_type pos_hc, size_type total_num_nodes,
    long double x_coordinate,
    long double y_coordinate,
    std::string node_text = std::string()) {
  return GetDefineNodeLineHc<size_type>(pos_hc, total_num_nodes, x_coordinate, y_coordinate, node_text);
}

template<typename size_type = std::size_t>
std::string GetNodeLines(size_type total_num_nodes,
    const std::vector<long double> &x_coordinates,
    const std::vector<long double> &y_coordinates,
    std::vector<std::string> node_texts = std::vector<std::string>()) {
  std::stringstream strm;
  for (size_type i = 0; i < total_num_nodes; i++) {
    std::string node_text = std::string();
    if (i < static_cast<size_type>(node_texts.size()))
      node_text = node_texts[i];
    strm << GetDefineNodeLine<size_type>(i, total_num_nodes, x_coordinates[i], y_coordinates[i], node_text);
  }
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetPathDefinitionTextHC(size_type total_num_nodes, size_type from, size_type to, bool is_solid, long double start_x, long double end_x,
    long double start_y, long double end_y,
    const std::vector<long double> &x_coordinates, const std::vector<long double> &y_coordinates) {
  std::stringstream strm;
  strm << std::left;

  int node_width = 2 * std::to_string(total_num_nodes).length() + 4;
  const auto total_num_nodes_mod3 = total_num_nodes % 3;
  const bool is_start_of_double_arrow = total_num_nodes_mod3 == 2 && from == (total_num_nodes - 2) / 3;
  assert(!(total_num_nodes_mod3 == 2 && from == (2 * total_num_nodes - 1) / 3));
  if (is_start_of_double_arrow)
    strm << '\n';

  long double amount_to_bend_start_down = 0.0l;
  long double amount_to_bend_start_up   = 0.0l;
  std::stringstream bend_strm;

  if (is_start_of_double_arrow && total_num_nodes % 2 == 1) {
    if (y_coordinates[from] <= y_coordinates[to]) {
      amount_to_bend_start_down = 20.0l;
      bend_strm << "bend right=" << std::left << amount_to_bend_start_down;
    } else {
      amount_to_bend_start_up = 20.0l;
      bend_strm << "bend left =" << std::left << amount_to_bend_start_up;
    }
  } else if (to == LeftChild(from) && y_coordinates[from] < y_coordinates[to]) {
    if (to + 1 < total_num_nodes && y_coordinates[to] >= y_coordinates[to + 1]) {
      amount_to_bend_start_down = 50.0l;
      bend_strm << "bend right=" << std::left << amount_to_bend_start_down;
    }
  } else if (to == RightChild(from) && y_coordinates[from] > y_coordinates[to]) {
    if (to > 0 && y_coordinates[to] <= y_coordinates[to - 1]) {
      amount_to_bend_start_up = 50.0l;
      bend_strm << "bend left=" << std::left << amount_to_bend_start_up;
    }
  }
  std::string bend_string = bend_strm.str();

  strm << std::setw(node_width) << GetNodeNameHc<size_type>(from, total_num_nodes, true) << " edge";

  std::stringstream options_strm;
  int edge_options_width = 28;
  if (is_start_of_double_arrow) {
    options_strm << "[<->,dashed";
    if (!bend_string.empty())
      options_strm << "," << bend_string;
    options_strm << "]";
  } else {
    if (!is_solid) {
      options_strm << "[dashed";
      if (!bend_string.empty())
        options_strm << "," << bend_string;
      options_strm << "]";
    } else {
      if (!bend_string.empty()) {
        options_strm << "[";
        if (!bend_string.empty())
          options_strm << bend_string;
        options_strm << "]";
      }
    }
  }
  strm << std::left << std::setw(edge_options_width) << options_strm.str();

  strm << " node[";
  if (from < to) {
    strm << "right]{}";
  } else {
    strm << "left]{} ";
  }
  strm << " ";
  strm << GetNodeNameHc<size_type>(to, total_num_nodes, true);
  strm << '\n';
  if (is_start_of_double_arrow)
    strm << '\n';
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetPathDefinitionTextMC(size_type total_num_nodes, size_type from_mc, size_type to_mc, bool is_solid, long double start_x, long double end_x,
    long double start_y, long double end_y,
    const std::vector<long double> &x_coordinates, const std::vector<long double> &y_coordinates) {
  std::stringstream strm;
  strm << std::left;
  size_type from_hc = FLIP(from_mc);
  size_type to_hc = FLIP(to_mc);

  int node_width = 2 * std::to_string(total_num_nodes).length() + 4;
  const auto total_num_nodes_mod3 = total_num_nodes % 3;
  const bool is_start_of_double_arrow = total_num_nodes_mod3 == 2 && from_mc == (total_num_nodes - 2) / 3;
  if (is_start_of_double_arrow) //If it is the max heap start of a double arrow then don't draw a second path.
    return std::string("\n");

  long double amount_to_bend_start_down = 0.0l;
  long double amount_to_bend_start_up   = 0.0l;
  std::string bend_string;
  if (to_mc == LeftChild(from_mc) && y_coordinates[from_hc] > y_coordinates[to_hc]) {
    if (to_mc + 1 < total_num_nodes && y_coordinates[to_hc] <= y_coordinates[FLIP(to_mc + 1)]) {
      amount_to_bend_start_down = 50.0l;
      std::stringstream bend_strm;
      bend_strm << "bend right=" << std::left << amount_to_bend_start_down;
      bend_string = bend_strm.str();
    }
  }
  if (to_mc == RightChild(from_mc) && y_coordinates[from_hc] < y_coordinates[to_hc]) {
    if (to_hc > 0 && y_coordinates[to_hc] >= y_coordinates[FLIP(to_mc - 1)]) {
      amount_to_bend_start_up = 50.0l;
      std::stringstream bend_strm;
      bend_strm << "bend left=" << std::left << amount_to_bend_start_up;
      bend_string = bend_strm.str();
    }
  }

  strm << std::setw(node_width) << GetNodeNameMc<size_type>(from_mc, total_num_nodes, true) << " edge";

  std::stringstream options_strm;
  int edge_options_width = 28;
  if (!is_solid) {
    options_strm << "[dashed";
    if (!bend_string.empty())
      options_strm << "," << bend_string;
    options_strm << "]";
  } else {
    if (!bend_string.empty()) {
      options_strm << "[";
      if (!bend_string.empty())
        options_strm << bend_string;
      options_strm << "]";
    }
  }
  strm << std::left << std::setw(edge_options_width) << options_strm.str();

  strm << " node[";
  if (from_mc < to_mc) {
    strm << "left]{} ";
  } else {
    strm << "right]{}";
  }
  strm << " ";
  strm << GetNodeNameMc<size_type>(to_mc, total_num_nodes, true);
  strm << '\n';
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetMinHeapPathDefinitionTextHC(size_type total_num_nodes, size_type from_hc, bool to_left_child, long double start_x, long double end_x,
    long double start_y, long double end_y,
    const std::vector<long double> &x_coordinates, const std::vector<long double> &y_coordinates) {
  auto heap_size = HeapSize<size_type>(total_num_nodes);
  auto pure_heap_size = (total_num_nodes / 2) + (total_num_nodes % 2);
  size_type end_node_hc = LeftChild<size_type>(from_hc);
  if (!to_left_child)
    end_node_hc++;
  if (end_node_hc >= heap_size)
    return std::string();
  bool should_be_solid = end_node_hc < pure_heap_size;

  return GetPathDefinitionTextHC<size_type>(total_num_nodes, from_hc, end_node_hc, should_be_solid, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);
}

template<typename size_type = std::size_t>
std::string GetMinHeapPathDefinitionTextMC(size_type total_num_nodes, size_type from_mc, bool to_left_child, long double start_x, long double end_x,
    long double start_y, long double end_y,
    const std::vector<long double> &x_coordinates, const std::vector<long double> &y_coordinates) {
  auto heap_size = HeapSize<size_type>(total_num_nodes);
  auto pure_heap_size = (total_num_nodes / 2) + (total_num_nodes % 2);
  size_type end_node_mc = LeftChild(from_mc);
  if (!to_left_child)
    end_node_mc++;
  if (end_node_mc >= heap_size)
    return std::string();
  bool should_be_solid = end_node_mc < pure_heap_size;

  return GetPathDefinitionTextMC<size_type>(total_num_nodes, from_mc, end_node_mc, should_be_solid, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);
}

template<typename size_type = std::size_t>
std::string GetPathDefinitions(size_type total_num_nodes, long double start_x, long double end_x,
            long double start_y, long double end_y,
            const std::vector<long double> &x_coordinates, const std::vector<long double> &y_coordinates) {
  std::stringstream strm;
  //Get the min heap paths.
  for (size_type i = 0; i < total_num_nodes; i++) {
    auto string = GetMinHeapPathDefinitionTextHC<size_type>(total_num_nodes, i, true, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);//Left child
    if (string.empty())
      break;
    strm << string;
    string = GetMinHeapPathDefinitionTextHC<size_type>(total_num_nodes, i, false, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);//Right child
    if (string.empty())
      break;
    strm << string;
  }

  //Get the max heap paths.
  for (size_type i = 0; i < total_num_nodes; i++) {
    auto string = GetMinHeapPathDefinitionTextMC<size_type>(total_num_nodes, i, true, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);//Left child
    if (string.empty())
      break;
    strm << string;
    string = GetMinHeapPathDefinitionTextMC<size_type>(total_num_nodes, i, false, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);//Right child
    if (string.empty())
      break;
    strm << string;
  }
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetTikzGraph(size_type total_num_nodes,
    long double start_x, long double end_x, long double start_y, long double end_y,
    std::vector<std::string> node_texts = std::vector<std::string>(),
    bool should_tighten_coordinates = true,
    long double x_middle_separation_distance_ratio = 0.05l,
    long double y_scale_down_ratio = 1.0,
    long double scale = 1.0,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    bool include_commented_out_FIGCO_definition = false) {
  std::stringstream strm;
  std::vector<long double> x_coordinates(total_num_nodes);
  std::vector<long double> y_coordinates(total_num_nodes);
  GetCoordinates<size_type>(total_num_nodes, start_x, end_x, start_y, end_y,
          x_coordinates, y_coordinates,
          x_middle_separation_distance_ratio, x_num_imaginary_middles_nodes_on_one_side, y_scale_down_ratio);
  if (should_tighten_coordinates)
    TightenCoordinates(total_num_nodes, start_x, end_x, start_y, end_y,
        x_coordinates, y_coordinates,
        x_middle_separation_distance_ratio, x_num_imaginary_middles_nodes_on_one_side, y_scale_down_ratio);
  if (include_commented_out_FIGCO_definition) {
    strm << "%\\newcommand{\\FIGCO}[2]{#1(#2)}\n";
  }
  strm << "\\begin{tikzpicture}[scale=" << scale << "]\n";
  strm << GetNodeLines<size_type>(total_num_nodes, x_coordinates, y_coordinates, node_texts);
  strm << "\\path[->,font=\\scriptsize,>=angle 90]\n";
  {
    std::string path_defs = GetPathDefinitions<size_type>(total_num_nodes, start_x, end_x, start_y, end_y, x_coordinates, y_coordinates);
    size_type len = path_defs.length();
    if (len > 1 && path_defs[len - 2] == '\n' && path_defs[len - 1] == '\n')
      path_defs.erase(path_defs.begin() + (len - 1)); //Only have one new line at the end.
    strm << path_defs;
  }
  strm << ";\n";
  strm << "\\end{tikzpicture}\n";
  return strm.str();
}

/*
 * Wrapper that makes the values in the vector "values" become the text in the tikz BiHeap graph.
 * Example call: GetTikzGraph(total_num_nodes, -std::log2(total_num_nodes/2), std::log2(total_num_nodes/2),
                              -total_num_nodes/2, total_num_nodes/2, vec);
 * Example call 2: std::cout << GetTikzGraph<int, std::size_t>(6, -4, 4, -4, 4, {0, 1, 2, 3, 4, 5}, true);
 *
 * Various  input parameters to this function determine the location of
 *  (and the amount of space in between) nearby nodes.
 * These include x_middle_separation_distance_ratio, y_scale_down_ratio, and x_num_imaginary_middles_nodes_on_one_side
 *
 * x_num_imaginary_middles_nodes_on_one_side - If this is set to 1.0l then the function pretends
 *     that there is an column of (imaginary and invisible) nodes in the pure heap, where
 *     these imaginary nodes would be the pure heap's last row.
 */
template<class T, typename size_type = std::size_t>
std::string GetTikzGraph(size_type total_num_nodes,
    long double start_x, long double end_x, long double start_y, long double end_y,
    const std::vector<T> values = std::vector<T>(),
    bool should_tighten_coordinates = true,
    long double x_middle_separation_distance_ratio = 0.05l,
    long double y_scale_down_ratio = 1.0,
    long double scale = 1.0,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    bool include_commented_out_FIGCO_definition = false) {
  std::vector<std::string> node_texts(values.size());
  for(auto i = 0u; i < values.size(); i++) {
    node_texts[i] += std::string("$");
    node_texts[i] += std::to_string(values[i]);
    node_texts[i] += std::string("$");
  }
  return GetTikzGraph<size_type>(total_num_nodes, start_x, end_x, start_y, end_y, node_texts,
                      should_tighten_coordinates,
                      x_middle_separation_distance_ratio, y_scale_down_ratio,
                      scale, x_num_imaginary_middles_nodes_on_one_side,
                      include_commented_out_FIGCO_definition);
}

//Examples of calling GetTikzGraph()
void PrintTikzGraphsExampleCalls() {
  std::cout << GetTikzGraph<int, std::size_t>(45, -8.5, 8.5, -12, 12);

  std::cout << "\n";
  std::vector<int> vec = {0,11,1,12,13,2,3,14,15,16,17,4,5,6,7,18,19,8,9,20,10,21};
  std::cout << GetTikzGraph<int, std::size_t>(22, -6, 6, -3, 3, vec);

  bool tighten_coordinates = true;
  std::cout << "\n";
  std::cout << GetTikzGraph<std::size_t>(21, -8.5, -0.6, -5, 5, {}, tighten_coordinates, 0.05l, 1.0, 1.0, 1.5l);

  std::cout << "\n";
  vec = {0,23,1,24,25,2,3,26,27,28,29,4,5,6,7,30,31,32,33,34,35,36,37,8,9,10,11,12,13,14,15,38,39,40,41,16,17,18,19,42,43,20,21,44,22,45};
  std::cout << GetTikzGraph<int, std::size_t>(46, -6, 6, -3, 3, vec, tighten_coordinates);
  std::cout << std::endl;
  return ;
}

#undef FLIP

#endif /* BIHEAP_TIKZ_GRAPH_H_ */
