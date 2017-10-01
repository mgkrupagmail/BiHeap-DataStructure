/*
 * biheap_tikz_graph.h
 *
 *  Created on: Sep 26, 2017
 *      Author: Matthew Gregory Krupa
 */

#ifndef BIHEAP_TIKZ_GRAPH_H_
#define BIHEAP_TIKZ_GRAPH_H_

#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "biheap_common.h"

#ifndef FLIP
#define FLIP(a) (total_num_nodes - 1 - (a))
#endif

template<typename size_type = std::size_t>
std::string GetNodeNameHc(size_type pos_hc, size_type total_num_nodes, bool with_parentheses = true) {
  std::stringstream strm;
  if (with_parentheses)
    strm << "(";
  strm << "s" << total_num_nodes << "n" << pos_hc;
  if (with_parentheses)
    strm << ")";
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

  bool is_in_left_branch = IsInLeftBranch(pos, subtree_start);
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
  size_type max_depth_of_pure_heap_ignoring_middle = GetDepth(size_of_pure_heap_ignoring_middle - 1);
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
std::string GetDefineNodeLineHc(size_type pos_hc, size_type total_num_nodes, long double x, long double y, std::string node_text = std::string()) {
  std::stringstream strm;
  int width = 2 * std::to_string(total_num_nodes).length() + 4;
  strm << "\\node " << std::setw(width) << GetNodeNameHc(pos_hc, total_num_nodes, true);
  strm << "  at  (" << std::setprecision(5) << x << ", " << std::setprecision(5) << y << ") \t";
  strm << "{";
  if (node_text.empty())
    strm << GetNodeTextHc(pos_hc, total_num_nodes, false);
  else
    strm << node_text;
  strm << "};\n";

  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetDefineNodeLine(size_type pos_hc, size_type total_num_nodes, long double start_x, long double end_x,
    long double start_y, long double end_y,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    long double y_scale_down_ratio = 1.0,
    std::string node_text = std::string()) {
  assert(start_x <= end_x && start_y <= end_y && total_num_nodes >= 1);
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
  return GetDefineNodeLineHc<size_type>(pos_hc, total_num_nodes, x_coordinate, y_coordinate, node_text);
}

template<typename size_type = std::size_t>
std::string GetNodeLines(size_type total_num_nodes, long double start_x, long double end_x, long double start_y, long double end_y,
    long double x_middle_separation_distance_ratio = 0.1l,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l,
    long double y_scale_down_ratio = 1.0,
    std::vector<std::string> node_texts = std::vector<std::string>()) {
  assert(start_x <= end_x && start_y <= end_y && total_num_nodes >= 1);
  if (total_num_nodes == 1) {
    return GetDefineNodeLineHc<size_type>(0, 1, (start_x + end_x) / 2, (start_y + end_y) / 2);
  }
  std::stringstream strm;
  for (size_type i = 0; i < total_num_nodes; i++) {
    std::string node_text = std::string();
    if (i < node_texts.size())
      node_text = node_texts[i];
    strm << GetDefineNodeLine<size_type>(i, total_num_nodes, start_x, end_x, start_y, end_y,
        x_middle_separation_distance_ratio, x_num_imaginary_middles_nodes_on_one_side, y_scale_down_ratio,
        node_text);
  }
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetPathDefinitionTextHC(size_type total_num_nodes, size_type from, size_type to, bool is_solid) {
  std::stringstream strm;
  strm << std::left;

  int node_width = 2 * std::to_string(total_num_nodes).length() + 4;
  strm << std::setw(node_width) << GetNodeNameHc<size_type>(from, total_num_nodes, true) << " edge";
  if (!is_solid) {
    strm << "[dashed] ";
  } else {
    strm << "         ";
  }
  strm << " node[";
  if (from < to) {
    strm << "right]{}";
  } else {
    strm << "left]{} ";
  }
  strm << " ";
  strm << GetNodeNameHc<size_type>(to, total_num_nodes, true);
  strm << '\n';
  return strm.str();
}

template<typename size_type = std::size_t>
std::string GetPathDefinitionTextMC(size_type total_num_nodes, size_type from_mc, size_type to_mc, bool is_solid) {
  std::stringstream strm;
  strm << std::left;

  int node_width = 2 * std::to_string(total_num_nodes).length() + 4;
  strm << std::setw(node_width) << GetNodeNameMc<size_type>(from_mc, total_num_nodes, true) << " edge";
  if (!is_solid) {
    strm << "[dashed] ";
  } else {
    strm << "         ";
  }
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
std::string GetMinHeapPathDefinitionTextHC(size_type total_num_nodes, size_type from_hc, bool to_left_child) {
  auto heap_size = HeapSize(total_num_nodes);
  auto pure_heap_size = (total_num_nodes / 2) + (total_num_nodes % 2);
  //auto last_interior_node = (pure_heap_size - 1) / 2;
  size_type end_node_hc = LeftChild<size_type>(from_hc);
  if (!to_left_child)
    end_node_hc++;
  if (end_node_hc >= heap_size)
    return std::string();
  bool should_be_solid = end_node_hc < pure_heap_size;

  return GetPathDefinitionTextHC<size_type>(total_num_nodes, from_hc, end_node_hc, should_be_solid);
}

template<typename size_type = std::size_t>
std::string GetMinHeapPathDefinitionTextMC(size_type total_num_nodes, size_type from_mc, bool to_left_child) {
  auto heap_size = HeapSize(total_num_nodes);
  auto pure_heap_size = (total_num_nodes / 2) + (total_num_nodes % 2);
  size_type end_node_mc = LeftChild(from_mc);
  if (!to_left_child)
    end_node_mc++;
  if (end_node_mc >= heap_size)
    return std::string();
  bool should_be_solid = end_node_mc < pure_heap_size;

  return GetPathDefinitionTextMC<size_type>(total_num_nodes, from_mc, end_node_mc, should_be_solid);
}

template<typename size_type = std::size_t>
std::string GetPathDefinitions(size_type total_num_nodes) {
  std::stringstream strm;
  //Get the min heap paths.
  for (size_type i = 0; i < total_num_nodes; i++) {
    auto string = GetMinHeapPathDefinitionTextHC<size_type>(total_num_nodes, i, true);//Left child
    if (string.empty())
      break;
    strm << string;
    string = GetMinHeapPathDefinitionTextHC<size_type>(total_num_nodes, i, false);//Right child
    if (string.empty())
      break;
    strm << string;
  }

  //Get the max heap paths.
  for (size_type i = 0; i < total_num_nodes; i++) {
    auto string = GetMinHeapPathDefinitionTextMC<size_type>(total_num_nodes, i, true);//Left child
    if (string.empty())
      break;
    strm << string;
    string = GetMinHeapPathDefinitionTextMC<size_type>(total_num_nodes, i, false);//Right child
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
    bool include_commented_out_FIGCO_definition = false,
    long double x_middle_separation_distance_ratio = 0.05l,
    long double y_scale_down_ratio = 1.0,
    long double scale = 1.0,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l) {
  std::stringstream strm;
  if (include_commented_out_FIGCO_definition) {
    strm << "%\\newcommand{\\FIGCO}[2]{#1(#2)}\n";
  }
  strm << "\\begin{tikzpicture}[scale=" << scale << "]\n";
  strm << GetNodeLines<size_type>(total_num_nodes, start_x, end_x, start_y, end_y, x_middle_separation_distance_ratio,
      x_num_imaginary_middles_nodes_on_one_side, y_scale_down_ratio, node_texts);
  strm << "\\path[->,font=\\scriptsize,>=angle 90]\n";
  strm << GetPathDefinitions<size_type>(total_num_nodes);
  strm << ";\n";
  strm << "\\end{tikzpicture}\n";
  return strm.str();
}

/*
 * Wrapper that makes the values in the vector "values" become the text in the tikz BiHeap graph.
 * Example call: GetTikzGraph(total_num_nodes, -std::log2(total_num_nodes/2), std::log2(total_num_nodes/2),
                              -total_num_nodes/2, total_num_nodes/2, vec);
 * Example call 2: std::cout << GetTikzGraph<int, std::size_t>(6, -4, 4, -4, 4, {0, 1, 2, 3, 4, 5}, true);
 */
template<class T, typename size_type = std::size_t>
std::string GetTikzGraph(size_type total_num_nodes,
    long double start_x, long double end_x, long double start_y, long double end_y,
    const std::vector<T> values = std::vector<T>(),
    bool include_commented_out_FIGCO_definition = false,
    long double x_middle_separation_distance_ratio = 0.05l,
    long double y_scale_down_ratio = 1.0,
    long double scale = 1.0,
    long double x_num_imaginary_middles_nodes_on_one_side = 0.8l) {
  std::vector<std::string> node_texts(values.size());
  for(auto i = 0u; i < values.size(); i++) {
    node_texts[i] += std::string("$");
    node_texts[i] += std::to_string(values[i]);
    node_texts[i] += std::string("$");
  }
  return GetTikzGraph<size_type>(total_num_nodes, start_x, end_x, start_y, end_y, node_texts,
                      include_commented_out_FIGCO_definition,
                      x_middle_separation_distance_ratio, y_scale_down_ratio,
                      scale, x_num_imaginary_middles_nodes_on_one_side);
}


#endif /* BIHEAP_TIKZ_GRAPH_H_ */
