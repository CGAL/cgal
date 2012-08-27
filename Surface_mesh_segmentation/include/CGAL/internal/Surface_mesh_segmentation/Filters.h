#ifndef CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H
#define CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H

/**
 * @file Filters.h
 * @brief This file contains 2 filtering methods, which can be used as a template parameter for CGAL::internal::Surface_mesh_segmentation.
 *
 * Also filtering methods can be used by their-own for smoothing associated values to `CGAL Polyhedron`.
 */
#include <vector>
#include <map>
#include <algorithm>
#include <queue>

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
                  int window_size,
                  ValuePropertyMap values) const {
    typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
    typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;

    std::vector<double> smoothed_values; // holds smoothed values
    smoothed_values.reserve(mesh.size_of_facets());

    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      std::map<Facet_const_handle, int> neighbors;
      NeighborSelector()(facet_it, window_size,
                         neighbors); // gather neighbors in the window

      double total_sdf_value = 0.0, total_weight = 0.0;
      double current_sdf_value = values[facet_it];
      // calculate deviation for range weighting.
      double deviation = 0.0;
      for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
          it != neighbors.end(); ++it) {
        deviation += std::pow(values[it->first] - current_sdf_value, 2);
      }
      deviation = std::sqrt(deviation / neighbors.size());
      if(deviation == 0.0) {
        //this might happen. In case there is no neighbors (i.e. NeighborSelector() returns just the parameter facet)
        //                 . Or all neighbor facets have same sdf value.
        smoothed_values.push_back(current_sdf_value);
        continue;
      }
      for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
          it != neighbors.end(); ++it) {
        double spatial_weight = gaussian_function(it->second,
                                window_size / 2.0); // window_size => 2*sigma
        double domain_weight = gaussian_function(values[it->first] -  current_sdf_value,
                               1.5 * deviation);
        // we can use just spatial_weight for Gauissian filtering
        double weight = spatial_weight * domain_weight;

        total_sdf_value += values[it->first] * weight;
        total_weight += weight;
      }
      smoothed_values.push_back(total_sdf_value / total_weight);
    }
    // put smoothed values back again to values pmap.
    std::size_t index = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++index) {
      values[facet_it] = smoothed_values[index];
    }
  }
private:
  /** Gauissian function for weighting. */
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
                  int window_size,
                  ValuePropertyMap values) const {
    typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
    typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;

    std::vector<double> smoothed_values;
    smoothed_values.reserve(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      //Find neighbors and put their values into a list
      std::map<Facet_const_handle, int> neighbors;
      NeighborSelector()(facet_it, window_size, neighbors);

      std::vector<double> neighbor_values;
      neighbor_values.reserve(neighbors.size());
      for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
          it != neighbors.end(); ++it) {
        neighbor_values.push_back(values[it->first]);
      }
      // Find median.
      int half_neighbor_count = neighbor_values.size() / 2;
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
    std::size_t index = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++index) {
      values[facet_it] = smoothed_values[index];
    }
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
  void operator()(Facet_const_handle& facet, int max_level,
                  std::map<Facet_const_handle, int>& neighbors) const {
    typedef std::pair<Facet_const_handle, int> Facet_level_pair;
    std::queue<Facet_level_pair> facet_queue;
    facet_queue.push(Facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const Facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        Halfedge_around_facet_const_circulator facet_circulator =
          pair.first->facet_begin();
        do {
          if(!facet_circulator->opposite()->is_border()) {
            facet_queue.push(Facet_level_pair(facet_circulator->opposite()->facet(),
                                              pair.second + 1));
          }
        } while(++facet_circulator != pair.first->facet_begin());
      }
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
  void operator()(Facet_const_handle& facet, int max_level,
                  std::map<Facet_const_handle, int>& neighbors) const {
    typedef std::pair<Facet_const_handle, int> Facet_level_pair;
    std::queue<Facet_level_pair> facet_queue;
    facet_queue.push(Facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const Facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        Facet_const_handle facet = pair.first;
        Halfedge_const_iterator edge = facet->halfedge();
        do { // loop on three vertices of the facet
          Vertex_const_iterator vertex = edge->vertex();
          Halfedge_around_vertex_const_circulator vertex_circulator =
            vertex->vertex_begin();
          do { // for each vertex loop on incoming edges (through those edges loop on neighbor facets which includes the vertex)
            if(!vertex_circulator->is_border()) {
              facet_queue.push(Facet_level_pair(vertex_circulator->facet(), pair.second + 1));
            }
          } while(++vertex_circulator != vertex->vertex_begin());
        } while((edge = edge->next()) != facet->halfedge());
      }
      facet_queue.pop();
    }
  }
};
}//namespace internal
}//namespace CGAL
#endif //CGAL_SURFACE_MESH_SEGMENTATION_FILTERS_H