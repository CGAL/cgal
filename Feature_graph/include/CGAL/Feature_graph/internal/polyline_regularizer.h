// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_POLYLINE_REGULARIZER_H
#define CGAL_FEATURE_GRAPH_INTERNAL_POLYLINE_REGULARIZER_H

#include <CGAL/license/Feature_graph.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/number_utils.h>

#include <CGAL/Feature_graph/internal/graph_operation.h>
#include <CGAL/Feature_graph/internal/feature_graph_operation.h>

namespace CGAL
{
namespace Feature_graph
{
namespace internal
{

template <typename Vertex>
using Corner_to_lines = std::map<Vertex, std::vector<std::size_t> >;
template <typename Vertex>
using Line = std::vector<Vertex>;
template <typename Vertex>
using Polyline = std::vector<Line<Vertex> >;

enum class Line_status
{
  not_selected,
  isthmus,
  simple,
  other
};

struct Constant_line_distance_threshold
{
  double distance;

  explicit Constant_line_distance_threshold(double distance)
    : distance(distance)
  {
    CGAL_assertion(distance > 0);
    CGAL_assertion(std::isfinite(distance * distance));
  }

  template <typename Vertex>
  double operator()(const Vertex&, Line_status) const
  {
    return distance;
  }

  template <typename Vertex>
  double maximum_threshold(const Vertex&) const
  {
    return distance;
  }
};

/*!
\brief Selects short or geometrically duplicated polylines.

A line is selected if its length is at most the configured threshold, or if every
vertex of the line is within its threshold of another line. The distance test is
directed: it measures this line against each other line, not vice versa.

The regularizer also calls `collect_potentially_updated_lines()` after a line
changes or is removed. The selector must write every other line index whose
selection result could have changed to the supplied output iterator; it need not
write the changed line's own index. `should_delete_line()` receives the
topological status of the line.
*/
template <typename K, typename ThresholdProvider = Constant_line_distance_threshold>
class Double_line_selector
{
  ThresholdProvider threshold_provider_;

  template <typename Vertex>
  double threshold(const Vertex& vertex, Line_status status) const
  {
    return threshold_provider_(vertex, status);
  }

  template <typename Vertex>
  double update_threshold_squared(const Vertex& vertex) const
  {
    const double value = threshold_provider_.maximum_threshold(vertex);
    return value * value + 1.e-3 * value;
  }

  template <typename Vertex, typename PointMap>
  double squared_distance(const Vertex& a, const Vertex& b, const PointMap& points) const
  {
    return CGAL::to_double(CGAL::squared_distance(get(points, a), get(points, b)));
  }

  template <typename PointMap, typename Vertex>
  double line_length(const Line<Vertex>& line, const PointMap& points, double stop_after) const
  {
    double length = 0.;
    for (std::size_t i = 1; i < line.size(); ++i)
    {
      length += std::sqrt(squared_distance(line[i - 1], line[i], points));
      if (length > stop_after)
        return length;
    }
    return length;
  }

  template <typename PointMap, typename Vertex>
  double point_to_line_squared_distance(const Vertex& point,
                                        const Line<Vertex>& line,
                                        const PointMap& points) const
  {
    const auto& query = get(points, point);
    if (line.size() == 1)
      return CGAL::to_double(CGAL::squared_distance(query, get(points, line.front())));

    double min_distance = (std::numeric_limits<double>::max)();
    for (std::size_t i = 1; i < line.size(); ++i)
    {
      typename K::Segment_3 segment(get(points, line[i - 1]), get(points, line[i]));
      min_distance = (std::min)(min_distance,
        CGAL::to_double(CGAL::squared_distance(query, segment)));
    }
    return min_distance;
  }

  template <typename PointMap, typename Vertex>
  bool is_double_line_against_other_line(const Line<Vertex>& line,
                                         const Line<Vertex>& other_line,
                                         const PointMap& points,
                                         Line_status status) const
  {
    if (line.empty() || other_line.empty())
      return false;
    for (const Vertex& point : line)
    {
      const double distance_threshold = threshold(point, status);
      CGAL_precondition(distance_threshold >= 0);
      if (point_to_line_squared_distance(point, other_line, points) >
          distance_threshold * distance_threshold)
        return false;
    }
    return true;
  }

public:
  explicit Double_line_selector(ThresholdProvider threshold_provider)
    : threshold_provider_(std::move(threshold_provider))
  {}

  template <typename Graph, typename PointMap, typename Vertex>
  bool should_delete_line(const Graph&,
                          std::size_t line_index,
                          const Polyline<Vertex>& lines,
                          Line_status status,
                          const PointMap& point_map) const
  {
    const Line<Vertex>& line = lines[line_index];
    if (line.empty())
      return false;

    double minimum_threshold = (std::numeric_limits<double>::max)();
    for (const Vertex& vertex : line)
    {
      const double distance_threshold = threshold(vertex, status);
      CGAL_precondition(distance_threshold >= 0);
      minimum_threshold = (std::min)(minimum_threshold, distance_threshold);
    }
    if (line_length(line, point_map, minimum_threshold) <= minimum_threshold)
      return true;

    for (std::size_t other_line_index = 0; other_line_index < lines.size(); ++other_line_index)
      if (other_line_index != line_index &&
          is_double_line_against_other_line(
            line, lines[other_line_index], point_map, status))
        return true;
    return false;
  }

  // The relaxed threshold avoids missing candidates near the selection boundary.
  // TODO: Replace this fixed tolerance with a scale-aware comparison.
  // TODO: Use a spatial acceleration structure to avoid scanning all lines.
  template <typename PointMap, typename Vertex, typename OutputIterator>
  void collect_potentially_updated_lines(std::size_t line_index,
                                         const Polyline<Vertex>& lines,
                                         const PointMap& point_map,
                                         OutputIterator output) const
  {
    if (lines[line_index].empty())
      return;

    for (std::size_t other_line_index = 0; other_line_index < lines.size(); ++other_line_index)
    {
      if (other_line_index == line_index || lines[other_line_index].empty())
        continue;
      bool potentially_double = true;
      for (const Vertex& point : lines[other_line_index])
        if (point_to_line_squared_distance(point, lines[line_index], point_map) >
            update_threshold_squared(point))
        {
          potentially_double = false;
          break;
        }
      if (potentially_double)
        *output++ = other_line_index;
    }
  }
};

using Line_adjacency_graph =
  boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;

template <typename Graph>
struct Polyline_vertex_collector
{
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;

  Polyline<Vertex>& lines;

  Polyline_vertex_collector(Polyline<Vertex>& polylines)
    : lines(polylines)
  {}

  void start_new_polyline()
  {
    lines.emplace_back();
  }

  void add_node(Vertex vertex)
  {
    lines.back().push_back(vertex);
  }

  void end_polyline()
  {
    // Isolated vertices are not regularizable lines.
    if (lines.back().size() < 2)
      lines.pop_back();
  }
};

template <typename Graph>
Polyline<typename boost::graph_traits<Graph>::vertex_descriptor>
extract_polylines(const Graph& graph)
{
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  Polyline<Vertex> lines;
  Polyline_vertex_collector<Graph> collector(lines);
  CGAL::split_graph_into_polylines(graph, collector);
  return lines;
}

template <typename Graph>
auto get_neighbors(typename boost::graph_traits<Graph>::vertex_descriptor vertex, const Graph& graph)
{
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  std::set<Vertex> neighbors;
  for (auto edge : boost::make_iterator_range(out_edges(vertex, graph)))
  {
    const Vertex source_vertex = source(edge, graph);
    const Vertex target_vertex = target(edge, graph);
    neighbors.insert(source_vertex == vertex ? target_vertex : source_vertex);
  }
  return neighbors;
}

// returns a map that associates a vector of adjacent lines (index of line in `lines`) for each corner
template <typename Vertex>
Corner_to_lines<Vertex> make_corner_to_lines(
  const Polyline<Vertex>& lines)
{
  Corner_to_lines<Vertex> corner_to_lines;
  for (std::size_t line_index = 0; line_index < lines.size(); ++line_index)
  {
    const Line<Vertex>& line = lines[line_index];
    if (line.empty())
      continue;
    corner_to_lines[line.front()].push_back(line_index);
    if (line.size() > 1 && line.front() != line.back())
      corner_to_lines[line.back()].push_back(line_index);
  }
  return corner_to_lines;
}

template <typename Vertex>
void remove_line_from_corner_lines(std::size_t line_index,
                                   const Vertex& corner,
                                   Corner_to_lines<Vertex>& corner_to_lines)
{
  std::vector<std::size_t>& lines_at_corner = corner_to_lines[corner];
  const auto it = std::find(lines_at_corner.begin(), lines_at_corner.end(), line_index);
  if (it != lines_at_corner.end())
    lines_at_corner.erase(it);
}

template <typename Vertex>
Line_adjacency_graph make_line_adjacency(
  const Corner_to_lines<Vertex>& corner_to_lines,
  std::size_t number_of_lines)
{
  Line_adjacency_graph line_adjacency(number_of_lines);
  for (const auto& entry : corner_to_lines)
  {
    for (std::size_t line_index_a : entry.second)
      for (std::size_t line_index_b : entry.second)
        if (line_index_a < line_index_b)
          add_edge(line_index_a, line_index_b, line_adjacency);
  }
  return line_adjacency;
}

// true iff one corner is incident to this line only (floating loops are considered isthmuses)
template <typename Vertex>
bool is_line_isthmus(std::size_t line_index,
                     const Polyline<Vertex>& lines,
                     const Corner_to_lines<Vertex>& corner_to_lines)
{
  const Line<Vertex>& line = lines[line_index];
  if (line.size() < 2)
    return true;
  return corner_to_lines.at(line.front()).size() == 1 ||
         corner_to_lines.at(line.back()).size() == 1;
}

inline bool is_line_simple(std::size_t line_index,
                           const Line_adjacency_graph& line_adjacency)
{
  using Vertex = typename boost::graph_traits<Line_adjacency_graph>::vertex_descriptor;
  return is_vertex_locally_connected(static_cast<Vertex>(line_index), line_adjacency);
}

template <typename LineSelector, typename Graph, typename PointMap, typename Vertex>
Line_status compute_line_status(
  std::size_t line_index,
  const Polyline<Vertex>& lines,
  const Corner_to_lines<Vertex>& corner_to_lines,
  const Line_adjacency_graph& line_adjacency,
  const LineSelector& line_selector,
  const Graph& graph,
  const PointMap& point_map)
{
  Line_status status;
  if (is_line_isthmus(line_index, lines, corner_to_lines))
    status = Line_status::isthmus;
  else if (is_line_simple(line_index, line_adjacency))
    status = Line_status::simple;
  else
    status = Line_status::other;
  return line_selector.should_delete_line(graph, line_index, lines, status, point_map)
    ? status : Line_status::not_selected;
}

// returns indices of adjacent lines
template <typename Vertex>
std::set<std::size_t> get_line_adjacent_lines(
  std::size_t line_index,
  const Polyline<Vertex>& lines,
  const Corner_to_lines<Vertex>& corner_to_lines)
{
  std::set<std::size_t> adjacent_lines;
  const Line<Vertex>& line = lines[line_index];
  if (line.empty())
    return {};
  const auto front = corner_to_lines.find(line.front());
  if (front != corner_to_lines.end())
    adjacent_lines.insert(front->second.begin(), front->second.end());
  const auto back = corner_to_lines.find(line.back());
  if (back != corner_to_lines.end())
    adjacent_lines.insert(back->second.begin(), back->second.end());
  return adjacent_lines;
}

// move elements of `line_to_push` to `line_container`,
// deletes `line_to_push`.
template <typename Vertex>
void merge_lines_following_order(
  std::size_t line_container_index,
  bool reverse_line_container,
  std::size_t line_to_push_index,
  bool reverse_line_to_push,
  Polyline<Vertex>& lines,
  Corner_to_lines<Vertex>& corner_to_lines,
  Line_adjacency_graph& line_adjacency,
  std::set<std::size_t>& lines_updated_by_merge,
  std::set<std::size_t>& lines_deleted_by_merge,
  std::vector<Line_status>& lines_status)
{
  Line<Vertex>& line_container = lines[line_container_index];
  Line<Vertex>& line_to_push = lines[line_to_push_index];
  if (reverse_line_container)
    std::reverse(line_container.begin(), line_container.end());
  if (reverse_line_to_push)
    std::reverse(line_to_push.begin(), line_to_push.end());

  const Vertex common_corner = line_container.back();
  CGAL_assertion(common_corner == line_to_push.front());

  // move vertices and delete line.
  corner_to_lines[common_corner].clear();
  const Vertex other_corner = line_to_push.back();
  remove_line_from_corner_lines(line_to_push_index, other_corner, corner_to_lines);
  if (other_corner != line_container.front())
    corner_to_lines[other_corner].push_back(line_container_index);
  line_container.insert(line_container.end(), std::next(line_to_push.begin()), line_to_push.end());
  line_to_push.clear();

  // update the dual graph: add all neighbors of the line to push to the line container,
  // and disconnect the line to push
  using Dual_vertex = typename boost::graph_traits<Line_adjacency_graph>::vertex_descriptor;
  const Dual_vertex container = static_cast<Dual_vertex>(line_container_index);
  const Dual_vertex pushed = static_cast<Dual_vertex>(line_to_push_index);
  for (auto edge : boost::make_iterator_range(out_edges(pushed, line_adjacency)))
  {
    const Dual_vertex neighbor = target(edge, line_adjacency);
    if (neighbor != container)
      add_edge(container, neighbor, line_adjacency);
  }
  clear_vertex(pushed, line_adjacency);
  lines_updated_by_merge.insert(line_container_index);
  lines_deleted_by_merge.insert(line_to_push_index);
  lines_status[line_to_push_index] = Line_status::not_selected;
}

// merges `line_a` and `line_b` around the shared corner `corner`
template <typename Vertex>
void merge_lines_on_corner(
  std::size_t line_index_a,
  std::size_t line_index_b,
  const Vertex& corner,
  Polyline<Vertex>& lines,
  Corner_to_lines<Vertex>& corner_to_lines,
  Line_adjacency_graph& line_adjacency,
  std::set<std::size_t>& lines_updated_by_merge,
  std::set<std::size_t>& lines_deleted_by_merge,
  std::vector<Line_status>& lines_status)
{
  if (line_index_a == line_index_b)
    return;
  Line<Vertex>& line_a = lines[line_index_a];
  Line<Vertex>& line_b = lines[line_index_b];
  // do not merge loops (or isolated lines)
  if (line_a.size() < 2 || line_b.size() < 2 ||
      line_a.front() == line_a.back() || line_b.front() == line_b.back())
    return;

  // four way of merging depending on the orientation of A and B.
  if (line_a.front() == corner)
  {
    if (line_b.front() == corner)
    {
      // ( line_A.back() ) --<-- ( line_A.front() [corner_index] line_B.front() ) -->-- ( line_B.back() )
      // take the reverse of line_A and push line_B
      // result: ( line_A.front() ) ----- > ----- ( line_A.back() )
      merge_lines_following_order(line_index_a, true, line_index_b, false, lines,
        corner_to_lines, line_adjacency, lines_updated_by_merge,
        lines_deleted_by_merge, lines_status);
    }
    else if (line_b.back() == corner)
    {
      // ( line_A.back() ) --<-- ( line_A.front() [corner_index] line_B.back() ) --<-- ( line_B.front() )
      // take line_B and push line_A
      // result: ( line_B.back() ) ----- < ----- ( line_B.front() )
      merge_lines_following_order(line_index_b, false, line_index_a, false, lines,
        corner_to_lines, line_adjacency, lines_updated_by_merge,
        lines_deleted_by_merge, lines_status);
    }
  }
  else if (line_a.back() == corner)
  {
    if (line_b.front() == corner)
    {
      // ( line_A.front() ) -->-- ( line_A.back() [corner_index] line_B.front() ) -->-- ( line_B.back() )
      // take line_A and push line_B
      // result: ( line_A.front() ) ----- > ----- ( line_A.back() )
      merge_lines_following_order(line_index_a, false, line_index_b, false, lines,
        corner_to_lines, line_adjacency, lines_updated_by_merge,
        lines_deleted_by_merge, lines_status);
    }
    else if (line_b.back() == corner)
    {
      // ( line_A.front() ) -->-- ( line_A.back() [corner_index] line_B.back() ) --<-- ( line_B.front() )
      // take line_A and push the reverse of line_B
      // result: ( line_A.front() ) ----- > ----- ( line_A.back() )
      merge_lines_following_order(line_index_a, false, line_index_b, true, lines,
        corner_to_lines, line_adjacency, lines_updated_by_merge,
        lines_deleted_by_merge, lines_status);
    }
  }
}

template <typename Graph, typename Vertex>
void erase_vertex(Graph& graph, const Vertex& vertex)
{
  clear_vertex(vertex, graph);
  remove_vertex(vertex, graph);
}

template <typename Graph, typename Vertex>
void remove_graph_path(Graph& graph, const Line<Vertex>& line)
{
  if (line.empty())
    return;

  for (std::size_t i = 1; i < line.size(); ++i)
    remove_edge(line[i - 1], line[i], graph);

  for (std::size_t i = 1; i + 1 < line.size(); ++i)
    erase_vertex(graph, line[i]);

  std::set<Vertex> endpoints;
  endpoints.insert(line.front());
  endpoints.insert(line.back());
  for (const Vertex& endpoint : endpoints)
    if (degree(endpoint, graph) == 0)
      erase_vertex(graph, endpoint);
}

// assumes `middle` is in `line`
template <typename Graph, typename Vertex>
void collapse_graph_path_to_middle(Graph& graph, const Line<Vertex>& line, const Vertex& middle)
{
  if (line.empty())
    return;

  std::set<Vertex> line_vertices(line.begin(), line.end());
  CGAL_assertion(line_vertices.count(middle) > 0);

  for (std::size_t i = 1; i < line.size(); ++i)
    remove_edge(line[i - 1], line[i], graph);

  for (const Vertex& endpoint : { line.front(), line.back() })
    for (const Vertex& neighbor : get_neighbors<Graph>(endpoint, graph))
      if (line_vertices.count(neighbor) == 0 && !edge(middle, neighbor, graph).second)
        add_edge(middle, neighbor, graph);

  for (const Vertex& vertex : line_vertices)
    if (vertex != middle)
      erase_vertex(graph, vertex);
}

// remove from graph, clears `lines[line_index]`,
// also merges adjacent lines if necessary
template <typename Graph, typename Vertex>
void remove_line(
  Graph& graph,
  std::size_t line_index,
  Polyline<Vertex>& lines,
  Corner_to_lines<Vertex>& corner_to_lines,
  Line_adjacency_graph& line_adjacency,
  std::set<std::size_t>& lines_updated_by_merge,
  std::set<std::size_t>& lines_deleted_by_merge,
  std::vector<Line_status>& lines_status)
{
  Line<Vertex>& line = lines[line_index];
  if (line.size() < 2)
  {
    remove_graph_path(graph, line);
    line.clear();
    lines_status[line_index] = Line_status::not_selected;
    return;
  }

  const Vertex front_corner = line.front();
  const Vertex back_corner = line.back();
  remove_line_from_corner_lines(line_index, front_corner, corner_to_lines);
  remove_line_from_corner_lines(line_index, back_corner, corner_to_lines);

  const std::vector<std::size_t> front_neighbors = corner_to_lines[front_corner];
  if (front_neighbors.size() == 2)
    merge_lines_on_corner(front_neighbors.front(), front_neighbors.back(), front_corner,
      lines, corner_to_lines, line_adjacency, lines_updated_by_merge,
      lines_deleted_by_merge, lines_status);

  const std::vector<std::size_t> back_neighbors = corner_to_lines[back_corner];
  if (back_neighbors.size() == 2)
    merge_lines_on_corner(back_neighbors.front(), back_neighbors.back(), back_corner,
      lines, corner_to_lines, line_adjacency, lines_updated_by_merge,
      lines_deleted_by_merge, lines_status);

  remove_graph_path(graph, line);
  line.clear();
  using Dual_vertex = typename boost::graph_traits<Line_adjacency_graph>::vertex_descriptor;
  clear_vertex(static_cast<Dual_vertex>(line_index), line_adjacency);
  lines_status[line_index] = Line_status::not_selected;
}

template <typename Graph, typename Vertex>
void collapse_line_to_middle(
  Graph& graph,
  std::size_t line_index,
  Polyline<Vertex>& lines,
  Corner_to_lines<Vertex>& corner_to_lines,
  Line_adjacency_graph& line_adjacency,
  std::vector<Line_status>& lines_status)
{
  Line<Vertex>& line = lines[line_index];
  if (line.size() < 2)
    return;

  const Vertex middle = line[line.size() / 2];
  const Vertex front_corner = line.front();
  const Vertex back_corner = line.back();
  if (line.size() > 2)
    corner_to_lines[middle].clear();
  remove_line_from_corner_lines(line_index, front_corner, corner_to_lines);
  remove_line_from_corner_lines(line_index, back_corner, corner_to_lines);

  const std::set<std::size_t> adjacent_lines =
    get_line_adjacent_lines(line_index, lines, corner_to_lines);
  std::vector<std::size_t>& middle_adjacent_lines = corner_to_lines[middle];
  middle_adjacent_lines.clear();
  middle_adjacent_lines.insert(
    middle_adjacent_lines.begin(), adjacent_lines.begin(), adjacent_lines.end());
  if (middle != front_corner)
    corner_to_lines[front_corner].clear();
  if (middle != back_corner)
    corner_to_lines[back_corner].clear();

  for (std::size_t adjacent_line_index : middle_adjacent_lines)
  {
    Line<Vertex>& adjacent_line = lines[adjacent_line_index];
    if (adjacent_line.empty())
      continue;
    if (adjacent_line.front() == front_corner || adjacent_line.front() == back_corner)
      adjacent_line.front() = middle;
    if (adjacent_line.back() == front_corner || adjacent_line.back() == back_corner)
      adjacent_line.back() = middle;
  }

  using Dual_vertex = typename boost::graph_traits<Line_adjacency_graph>::vertex_descriptor;
  const Dual_vertex collapsed = static_cast<Dual_vertex>(line_index);
  // assumes that the iterators given by `out_edges` of `collapsed` are not invalidated
  // by `add_edge` on two other vertices.
  for (auto [edge_iterator, edge_end] = out_edges(collapsed, line_adjacency);
      edge_iterator != edge_end; edge_iterator++)
    for (auto next_edge_iterator = std::next(edge_iterator);
        next_edge_iterator != edge_end; next_edge_iterator++)
    {
      Dual_vertex first_neighbor = target(*edge_iterator, line_adjacency);
      Dual_vertex second_neighbor = target(*next_edge_iterator, line_adjacency);
      if (first_neighbor != collapsed && second_neighbor != collapsed)
        add_edge(first_neighbor, second_neighbor, line_adjacency);
    }
  clear_vertex(collapsed, line_adjacency);

  lines_status[line_index] = Line_status::not_selected;
  collapse_graph_path_to_middle(graph, line, middle);
  line.clear();
}

} // namespace internal

/*!
\brief Regularizes the polylines of a thinned feature graph using a line selector.

The method extracts the polylines once, then updates their stored paths and
adjacency as lines are removed, merged, or collapsed. It prioritizes isthmus
lines, then simple lines, then other selected lines.

\tparam LineSelector a model of the line-selection policy. It provides:
\code
template <typename Graph, typename PointMap, typename Vertex>
bool should_delete_line(const Graph&, std::size_t,
                        const internal::Polyline<Vertex>&,
                        internal::Line_status, const PointMap&) const;

template <typename PointMap, typename Vertex, typename OutputIterator>
void collect_potentially_updated_lines(std::size_t,
                                      const internal::Polyline<Vertex>&,
                                      const PointMap&, OutputIterator) const;
\endcode
The second method must output every line index whose selection result may
change after the specified line is changed or removed.

\param graph the feature graph, modified in place
\param point_map a readable property map assigning a point to each graph vertex;
                 it must remain valid while the graph is modified
\param line_selector the policy that selects lines for removal

\pre `Graph` must preserve the descriptors and vertex iterators of remaining
vertices when `remove_vertex()` is called.
*/
template <typename LineSelector, typename PointMap>
void regularize_lines(internal::Construction_feature_graph& graph,
                      const PointMap& point_map,
                      const LineSelector& line_selector)
{
  using Graph = internal::Construction_feature_graph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;

  // extract lines from graph for faster search and deletion
  //    also assigns an index for each line that is kept:
  //      line deletion empties the line or
  //      line merging moves line elements
  internal::Polyline<Vertex> lines = internal::extract_polylines(graph);
  // map that associates a vector of adjacent lines (index of line in `lines`) for each corner.
  // Updated when `lines` updates
  internal::Corner_to_lines<Vertex> corner_to_lines =
    internal::make_corner_to_lines(lines);
  // dual graph of lines to check the simple property on lines. Updated when `lines` updates
  internal::Line_adjacency_graph line_adjacency =
    internal::make_line_adjacency(corner_to_lines, lines.size());
  // cache the status of lines (not selected, isthmus, simple, or other), which changes when an
  // adjacent or close line is removed
  std::vector<internal::Line_status> line_status(lines.size(), internal::Line_status::not_selected);
  for (std::size_t i = 0; i < lines.size(); ++i)
    line_status[i] = internal::compute_line_status(
      i, lines, corner_to_lines, line_adjacency, line_selector, graph, point_map);

  // containers defined outside of the loop to reuse allocated memory.
  std::vector<std::size_t> isthmus_lines;
  std::vector<std::size_t> simple_lines;
  std::vector<std::size_t> other_lines;
  std::set<std::size_t> lines_to_delete;
  std::set<std::size_t> lines_to_update_status;
  std::set<std::size_t> lines_updated_by_merge;
  std::set<std::size_t> lines_deleted_by_merge;
  std::set<std::size_t> adjacent_lines;

  for (;;)
  {
    // fill lines with priority (other lines are not pushed if there is a simple line, and simple
    // lines are not pushed if there are isthmus lines)
    isthmus_lines.clear();
    simple_lines.clear();
    other_lines.clear();
    for (std::size_t line_index = 0; line_index < line_status.size(); ++line_index)
    {
      const internal::Line_status& current_status = line_status[line_index];
      if (current_status == internal::Line_status::isthmus)
        isthmus_lines.push_back(line_index);
      else if (isthmus_lines.empty() && current_status == internal::Line_status::simple)
        simple_lines.push_back(line_index);
      else if (isthmus_lines.empty() && simple_lines.empty() &&
               current_status == internal::Line_status::other)
        other_lines.push_back(line_index);
    }

    // select the lines to remove this iteration by their status.
    //    select in order: isthmus (all lines), simple (longest line), or any other line
    lines_to_delete.clear();
    internal::Line_status deletion_status;
    if (!isthmus_lines.empty())
    {
      deletion_status = internal::Line_status::isthmus;
      lines_to_delete.insert(isthmus_lines.begin(), isthmus_lines.end());
    }
    else if (!simple_lines.empty())
    {
      deletion_status = internal::Line_status::simple;
      std::size_t longest_line = simple_lines.front();
      for (std::size_t line_index : simple_lines)
        if (lines[line_index].size() > lines[longest_line].size())
          longest_line = line_index;
      lines_to_delete.insert(longest_line);
    }
    else if (!other_lines.empty())
    {
      deletion_status = internal::Line_status::other;
      lines_to_delete.insert(other_lines.front());
    }
    else
      break; // the loop ends here when there are no more lines to delete

    // collapse lines
    lines_to_update_status.clear();
    lines_updated_by_merge.clear();
    lines_deleted_by_merge.clear();
    adjacent_lines.clear();
    // Collect adjacent and selector-dependent lines whose statuses may change.
    for (std::size_t line_index : lines_to_delete)
    {
      const std::set<std::size_t> adjacent =
        internal::get_line_adjacent_lines(line_index, lines, corner_to_lines);
      adjacent_lines.insert(adjacent.begin(), adjacent.end());
      lines_to_update_status.insert(adjacent.begin(), adjacent.end());
      line_selector.collect_potentially_updated_lines(
        line_index, lines, point_map,
        std::inserter(lines_to_update_status, lines_to_update_status.end()));
    }
    // apply line collapses depending of the status.
    // Updates `lines`, `corner_to_lines`, `line_adjacency`, and `lines_status`.
    //  also marks the lines updated or deleted by the merge operation
    if (deletion_status == internal::Line_status::isthmus)
    {
      for (std::size_t line_index : lines_to_delete)
        if (lines_updated_by_merge.count(line_index) == 0)
          internal::remove_line(graph, line_index, lines, corner_to_lines,
            line_adjacency, lines_updated_by_merge, lines_deleted_by_merge, line_status);
    }
    else if (deletion_status == internal::Line_status::simple)
      internal::remove_line(graph, *lines_to_delete.begin(), lines, corner_to_lines,
        line_adjacency, lines_updated_by_merge, lines_deleted_by_merge, line_status);
    else
      internal::collapse_line_to_middle(graph, *lines_to_delete.begin(), lines,
        corner_to_lines, line_adjacency, line_status);

    // update line status for:
    //   - lines adjacent to a deleted line (already collected)
    //   - lines that are double because of a deleted line (already collected)
    //   - lines that are double because of an adjacent line (occurs following the merge operation)
    // Collect lines whose selection may have changed because an adjacent line was merged.
    for (std::size_t adjacent_line_index : adjacent_lines)
    {
      if (lines_to_delete.count(adjacent_line_index) != 0)
        continue;
      line_selector.collect_potentially_updated_lines(
        adjacent_line_index, lines, point_map,
        std::inserter(lines_to_update_status, lines_to_update_status.end()));
    }
    // then update status
    for (std::size_t line_index : lines_to_update_status)
    {
      if (lines[line_index].empty() || lines_to_delete.count(line_index) != 0)
        continue;
      line_status[line_index] = internal::compute_line_status(
        line_index, lines, corner_to_lines, line_adjacency, line_selector, graph, point_map);
    }

  }

  // finally, remove isolated vertices
  //  (only this part assumes that `remove_vertex` does not invalidate other vertex iterators)
  for (auto [vertex_iterator, vertex_end] = vertices(graph);
      vertex_iterator != vertex_end;)
  {
    const Vertex vertex = *vertex_iterator++; // increment iterator before erasing
    if (degree(vertex, graph) == 0)
      internal::erase_vertex(graph, vertex);
  }
}

/*!
\brief Regularizes short and geometrically duplicated lines.

This convenience overload invokes `regularize_lines()` with a
\ref CGAL::Feature_graph::internal::Double_line_selector.

\param graph the feature graph, modified in place
\param point_map a readable property map assigning a point to each graph vertex
\param collapse_double_lines_distance the maximum line length and directed
       vertex-to-polyline distance used by the selector; must be positive
*/
template <typename K, typename PointMap>
void regularize_double_lines(internal::Construction_feature_graph& graph,
                             const PointMap& point_map,
                             double collapse_double_lines_distance)
{
  regularize_lines(
    graph, point_map, internal::Double_line_selector<K>(
      internal::Constant_line_distance_threshold(collapse_double_lines_distance)));
}

} // namespace Feature_graph
} // namespace CGAL

#endif // CGAL_FEATURE_GRAPH_INTERNAL_POLYLINE_REGULARIZER_H
