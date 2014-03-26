#ifndef CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


#define CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H

/// @cond CGAL_DOCUMENT_INTERNAL
/**
 * @file Filters.h
 * @brief This file contains 2 filtering methods, which can be used as a template parameter for CGAL::internal::Surface_mesh_segmentation.
 *
 * Also filtering methods can be used by their-own for smoothing associated values with `CGAL Polyhedron` facets.
 */
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <cmath>

#include <boost/optional.hpp>

namespace CGAL
{

namespace internal
{

template<class Polyhedron>
class Neighbor_selector_by_edge;

template<class Polyhedron>
class Neighbor_selector_by_vertex;

/** @brief Applies bilateral filtering on values which are associated with polyhedron facets. */
template <class Polyhedron, class NeighborSelector = Neighbor_selector_by_edge<Polyhedron> >
class Bilateral_filtering
{
public:
  /**
   * Bilateral filtering for values associated with facets.
   * For each facet, takes neighbors in @a window_size, and assigns weighted average of neighbor values as filtered result.
   * For weighting two weights are multiplied:
   *   - spatial: over geodesic distances (number of edges)
   *   - domain : over value distances
   * @param mesh `CGAL Polyhedron` on which @a values are defined
   * @param window_size range of effective neighbors
   * @param[in, out] values `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   */
  template<class ValuePropertyMap>
  void operator()(const Polyhedron& mesh,
                  std::size_t window_size,
                  ValuePropertyMap values,
                  boost::optional<double> spatial_parameter = boost::optional<double>(),
                  boost::optional<double> range_parameter = boost::optional<double>()
                 ) const {
    typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
    typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;

    double spatial_parameter_actual;
    if(!spatial_parameter) {
      spatial_parameter_actual = window_size / 2.0;
    } else                   {
      spatial_parameter_actual = *range_parameter;
    }

    std::vector<double> smoothed_values; // holds smoothed values
    smoothed_values.reserve(mesh.size_of_facets());

    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      std::map<Facet_const_handle, std::size_t> neighbors;
      NeighborSelector()(facet_it, window_size,
                         neighbors); // gather neighbors in the window
      double current_sdf_value = values[facet_it];

      double range_parameter_actual;
      if(!range_parameter) {
        // calculate deviation for range weighting.
        double deviation = 0.0;
        for(typename std::map<Facet_const_handle, std::size_t>::iterator it =
              neighbors.begin(); it != neighbors.end(); ++it) {
          deviation += std::pow(values[it->first] - current_sdf_value, 2);
        }
        deviation = std::sqrt(deviation / neighbors.size());
        if(deviation == 0.0) {
          //this might happen. In case there is no neighbors (i.e. NeighborSelector() returns just the parameter facet)
          //                 . Or all neighbor facets have same sdf value.
          smoothed_values.push_back(current_sdf_value);
          continue;
        }
        range_parameter_actual = 1.5 * deviation;
      } else {
        range_parameter_actual = *range_parameter;
      }

      // smooth
      double total_sdf_value = 0.0, total_weight = 0.0;
      for(typename std::map<Facet_const_handle, std::size_t>::iterator it =
            neighbors.begin(); it != neighbors.end(); ++it) {
        double spatial_weight = gaussian_function(static_cast<double>(it->second),
                                spatial_parameter_actual);
        double range_weight = gaussian_function(values[it->first] - current_sdf_value,
                                                range_parameter_actual);
        // we can use just spatial_weight for Gaussian filtering
        double weight = spatial_weight * range_weight;

        total_sdf_value += values[it->first] * weight;
        total_weight += weight;
      }
      smoothed_values.push_back(total_sdf_value / total_weight);
    }
    // put smoothed values back again to values pmap.
    std::vector<double>::iterator smoothed_value_it = smoothed_values.begin();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++smoothed_value_it) {
      values[facet_it] = *smoothed_value_it;
    }
  }
private:
  /** Gaussian function for weighting. */
  double gaussian_function(double value, double deviation) const {
    return exp(-0.5 * (std::pow(value / deviation, 2)));
  }
};

/** @brief Applies median filtering on values which are associated with polyhedron facets. */
template <class Polyhedron, class NeighborSelector = Neighbor_selector_by_vertex<Polyhedron> >
class Median_filtering
{
public:
  /**
   * Median filtering for values associated with facets.
   * For each facet, takes neighbors in @a window_size, and assigns median of values of neighbors as filtered result.
   *
   * @param mesh `CGAL Polyhedron` on which @a values are defined
   * @param window_size range of effective neighbors
   * @param[in, out] values `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   */
  template<class ValuePropertyMap>
  void operator()(const Polyhedron& mesh,
                  std::size_t window_size,
                  ValuePropertyMap values) const {
    typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
    typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;

    std::vector<double> smoothed_values;
    smoothed_values.reserve(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      std::map<Facet_const_handle, std::size_t> neighbors;
      NeighborSelector()(facet_it, window_size,
                         neighbors); // gather neighbors in the window

      std::vector<double> neighbor_values;
      neighbor_values.reserve(neighbors.size());
      for(typename std::map<Facet_const_handle, std::size_t>::iterator it =
            neighbors.begin(); it != neighbors.end(); ++it) {
        neighbor_values.push_back(values[it->first]);
      }
      // Find median.
      std::size_t half_neighbor_count = neighbor_values.size() / 2;
      std::nth_element(neighbor_values.begin(),
                       neighbor_values.begin() + half_neighbor_count, neighbor_values.end());
      double median_sdf = neighbor_values[half_neighbor_count];
      if(neighbor_values.size() % 2 == 0) {
        median_sdf += *std::max_element(neighbor_values.begin(),
                                        neighbor_values.begin() + half_neighbor_count);
        median_sdf /= 2;
      }
      smoothed_values.push_back(median_sdf);
    }
    // put smoothed values back again to values pmap.
    std::vector<double>::iterator smoothed_value_it = smoothed_values.begin();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it) {
      values[facet_it] = *smoothed_value_it;
    }
  }
};

/** @brief A dummy filter to use as template parameter of `Surface_mesh_segmentation` class that does nothing. */
struct No_filtering {
  /**
   * empty implementation of required operator.
   */
  template<class Polyhedron,class ValuePropertyMap>
  void operator()(const Polyhedron& /* mesh */,
                  std::size_t /* window_size */,
                  ValuePropertyMap /* values */) const {
  }
};

/** @brief A filter that applies the filter passed as template parameter several times. */
template <class Filter, std::size_t nb_iterations = 5>
struct Iterative_filter : public Filter {
  /**
   * empty implementation of required operator.
   */
  template<class Polyhedron,class ValuePropertyMap>
  void operator()(const Polyhedron&  mesh ,
                  std::size_t  window_size ,
                  ValuePropertyMap  values ) const {
    for (std::size_t i=0; i<nb_iterations; ++i)
      Filter::operator()(mesh, window_size, values);
  }
};

/** @brief Gathers neighbors of a facet for a given window range. @see Bilateral_filtering, Median_filtering */
template<class Polyhedron>
class Neighbor_selector_by_edge
{
private:
  typedef typename Polyhedron::Facet::Halfedge_around_facet_const_circulator
  Halfedge_around_facet_const_circulator;
public:
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  /**
   * Breadth-first traversal on facets by treating facets, which share a common edge, are 1-level neighbors.
   * @param facet root facet
   * @param max_level maximum allowed distance (number of levels) between root facet and visited facet
   * @param[out] neighbors visited facets and their distances to root facet
   */
  void operator()(Facet_const_handle facet,
                  std::size_t max_level,
                  std::map<Facet_const_handle, std::size_t>& neighbors) const {
    typedef std::pair<Facet_const_handle, std::size_t> Facet_level_pair;

    std::queue<Facet_level_pair> facet_queue;
    facet_queue.push(Facet_level_pair(facet, 0));
    neighbors.insert(facet_queue.front());

    if(max_level <= 0) {
      return;
    }

    while(!facet_queue.empty()) {
      const Facet_level_pair& pair = facet_queue.front();

      Halfedge_around_facet_const_circulator facet_circulator =
        pair.first->facet_begin();
      do {
        if(!facet_circulator->opposite()->is_border()) {
          Facet_level_pair new_pair(facet_circulator->opposite()->facet(),
                                    pair.second + 1);
          if(neighbors.insert(new_pair).second
              && new_pair.second < max_level) { // first insert new_pair to map
            // if insertion is OK, then check its level
            facet_queue.push(
              new_pair);                                      // if its level is equal to max_level do not put it in
          }                                                                    // queue since we do not want to traverse its neighbors
        }
      } while(++facet_circulator != pair.first->facet_begin());

      facet_queue.pop();
    }
  }
};

/** @brief Gathers neighbors of a facet for a given window range. @see Bilateral_filtering, Median_filtering */
template<class Polyhedron>
class Neighbor_selector_by_vertex
{
private:
  typedef typename Polyhedron::Facet::Halfedge_around_vertex_const_circulator
  Halfedge_around_vertex_const_circulator;
  typedef typename Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Polyhedron::Vertex_const_iterator   Vertex_const_iterator;
public:
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  /**
   * Breadth-first traversal on facets by treating facets, which share a common vertex, are 1-level neighbors.
   * @param facet root facet
   * @param max_level maximum allowed distance (number of levels) between root facet and visited facet
   * @param[out] neighbors visited facets and their distances to root facet
   */
  void operator()(Facet_const_handle facet,
                  std::size_t max_level,
                  std::map<Facet_const_handle, std::size_t>& neighbors) const {
    typedef std::pair<Facet_const_handle, std::size_t> Facet_level_pair;

    std::queue<Facet_level_pair> facet_queue;
    facet_queue.push(Facet_level_pair(facet, 0));
    neighbors.insert(facet_queue.front());

    if(max_level <= 0) {
      return;
    }

    while(!facet_queue.empty()) {
      const Facet_level_pair& pair = facet_queue.front();

      Facet_const_handle facet_front = pair.first;
      Halfedge_const_iterator edge = facet_front->halfedge();
      do { // loop on three vertices of the facet
        Vertex_const_iterator vertex = edge->vertex();
        Halfedge_around_vertex_const_circulator vertex_circulator =
          vertex->vertex_begin();
        do { // for each vertex loop on incoming edges (through those edges loop on neighbor facets which includes the vertex)
          if(!vertex_circulator->is_border()) {
            Facet_level_pair new_pair(vertex_circulator->opposite()->facet(),
                                      pair.second + 1);
            if(neighbors.insert(new_pair).second
                && new_pair.second < max_level) { // first insert new_pair to map
              // if insertion is OK, then check its level
              facet_queue.push(
                new_pair);                                      // if its level is equal to max_level do not put it in
            }                                                                    // queue since we do not want to traverse its childs
          }
        } while(++vertex_circulator != vertex->vertex_begin());
      } while((edge = edge->next()) != facet_front->halfedge());

      facet_queue.pop();
    }
  }
};
}//namespace internal
/// @endcond
}//namespace CGAL
#endif //CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H
