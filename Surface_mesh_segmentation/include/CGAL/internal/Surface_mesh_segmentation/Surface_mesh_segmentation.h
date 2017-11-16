#ifndef CGAL_SURFACE_MESH_SEGMENTATION_H
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
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Ilker O. Yaz


#define CGAL_SURFACE_MESH_SEGMENTATION_H

#include <CGAL/license/Surface_mesh_segmentation.h>


#include <CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h>
#include <CGAL/internal/Surface_mesh_segmentation/Filters.h>
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>
#include <CGAL/internal/Surface_mesh_segmentation/SDF_calculation.h>

#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/property_map.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>

#define CGAL_NORMALIZATION_ALPHA 5.0
#define CGAL_CONVEX_FACTOR 0.08
#define CGAL_SMOOTHING_LAMBDA_MULTIPLIER 100.0

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{

// Post-process functions for sdf values
template<class Polyhedron>
class Postprocess_sdf_values
{
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator face_iterator;

  typedef Bilateral_filtering<Polyhedron>           Default_filter;

public:

  template <class Filter, class SDFPropertyMap>
  std::pair<double, double>
  postprocess(const Polyhedron& mesh, SDFPropertyMap sdf_pmap) {
    check_missing_sdf_values(mesh, sdf_pmap);
    Filter()(mesh, get_window_size(mesh), sdf_pmap);
    return linear_normalize_sdf_values(mesh, sdf_pmap);
  }

  template <class SDFPropertyMap>
  std::pair<double, double>
  postprocess(const Polyhedron& mesh, SDFPropertyMap sdf_pmap) {
    return postprocess<Default_filter>(mesh, sdf_pmap);
  }

  /**
   * Finds facets which have missing (indicated by -1) sdf values.
   * Sdf values on these facets are assigned to average sdf value of its neighbors.
   * If still there is any facet which has no sdf value, assigns minimum sdf value to it.
   * This is meaningful since (being an outlier) zero sdf values might effect normalization & log extremely.
   * @param[in, out] sdf_values `ReadWritePropertyMap` with `Polyhedron::face_descriptor` as key and `double` as value type
   */
  template<class SDFPropertyMap>
  void check_missing_sdf_values(const Polyhedron& mesh,
                                SDFPropertyMap sdf_values) {
    std::vector<face_descriptor> still_missing_facets;
    double min_sdf = (std::numeric_limits<double>::max)();
    // If there is any facet which has no sdf value, assign average sdf value of its neighbors
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      double sdf_value = get(sdf_values, *facet_it);
      CGAL_assertion(sdf_value == -1 || sdf_value >= 0); // validity check
      if(sdf_value != -1.0) {
        min_sdf = (std::min)(sdf_value, min_sdf);
        continue;
      }

      Halfedge_around_face_circulator<Polyhedron> facet_circulator(halfedge(*facet_it,mesh),mesh), done(facet_circulator);

      double total_neighbor_sdf = 0.0;
      std::size_t nb_valid_neighbors = 0;
      do {
        if(!(face(opposite(*facet_circulator,mesh),mesh) == boost::graph_traits<Polyhedron>::null_face())) {
            double neighbor_sdf = get(sdf_values, face(opposite(*facet_circulator,mesh),mesh));
          if(neighbor_sdf != -1) {
            total_neighbor_sdf += neighbor_sdf;
            ++nb_valid_neighbors;
          }
        }
      } while( ++facet_circulator != done);

      if(nb_valid_neighbors == 0) {
        still_missing_facets.push_back(*facet_it);
      } else {
        sdf_value = total_neighbor_sdf / nb_valid_neighbors;
        put(sdf_values, *facet_it, sdf_value);
        // trying to update min_sdf is pointless, since it is interpolated one of the neighbors sdf will be smaller than it
      }
    }
    // If still there is any facet which has no sdf value, assign minimum sdf value.
    // This is meaningful since (being an outlier) 0 sdf values might effect normalization & log extremely.
    for(typename std::vector<face_descriptor>::iterator it =
          still_missing_facets.begin();
        it != still_missing_facets.end(); ++it) {
      put(sdf_values, *it, min_sdf);
    }
  }

  template<class SDFPropertyMap>
  std::pair<double, double> min_max_value(const Polyhedron& mesh,
                                          SDFPropertyMap sdf_values) {
    double min_sdf = (std::numeric_limits<double>::max)();
    double max_sdf = -min_sdf;
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      double sdf_value = get(sdf_values, *facet_it);
      max_sdf = (std::max)(sdf_value, max_sdf);
      min_sdf = (std::min)(sdf_value, min_sdf);
    }
    return std::make_pair(min_sdf, max_sdf);
  }
  /**
   * Normalize sdf values between [0-1].
   * @param sdf_values `ReadWritePropertyMap` with `Polyhedron::face_descriptor` as key and `double` as value type
   * @return minimum and maximum SDF values before normalization
   */
  template<class SDFPropertyMap>
  std::pair<double, double> linear_normalize_sdf_values(const Polyhedron& mesh,
      SDFPropertyMap sdf_values) {
    double min_sdf, max_sdf;
    boost::tie(min_sdf, max_sdf) = min_max_value(mesh, sdf_values);

    if(min_sdf == max_sdf) {
      CGAL_warning(min_sdf == max_sdf && !"Linear normalization is not applicable!");
      return std::make_pair(min_sdf, max_sdf);
    }

    const double max_min_dif = max_sdf - min_sdf;
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      put(sdf_values, *facet_it, (get(sdf_values, *facet_it) - min_sdf) / max_min_dif);
    }
    return std::make_pair(min_sdf, max_sdf);
  }

  /**
   * Simple window-size determination function for smoothing.
   * It is proportional to square root of size of facets in polyhedron.
   * @return size of the window
   *  - 0-2000     -> 1
   *  - 2000-8000  -> 2
   *  - 8000-18000 -> 3
   *  - ...
   */
  std::size_t get_window_size(const Polyhedron& mesh) {
    double facet_sqrt = std::sqrt(num_faces(mesh) / 2000.0);
    return static_cast<std::size_t>(facet_sqrt) + 1;
  }
};

/**
 * @brief Main entry point for mesh segmentation algorithm.
 *
 * It is a connector class which uses:
 *   -  SDF_calculation for calculating sdf values
 *   -  Expectation_maximization for soft clustering
 *   -  An implementation of alpha-expansion graph cut for hard clustering
 *
 * Other than being a connector, it is also responsible for pre-process and postprocess on intermediate data, which are:
 *   - log-normalizing probabilities received from soft clustering
 *   - log-normalizing and calculating dihedral-angle based weights for edges
 *   - smoothing and log-normalizing sdf values received from sdf calculation (Filters.h)
 *   - assigning segment-id for each facet after hard clustering
 *
 * @tparam Polyhedron a CGAL polyhedron
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @tparam GraphCut chosen graph-cut implementation from Alpha_expansion_graph_cut.h
 * @tparam Filter chosen filtering method from Filters.h
 */
template <
class Polyhedron,
      class GeomTraits,
  class VertexPointPmap,
      bool fast_bbox_intersection = true,
#ifndef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
      class GraphCut = Alpha_expansion_graph_cut_boykov_kolmogorov,
#else
      class GraphCut = Alpha_expansion_graph_cut_boost,
#endif
      class Filter = Bilateral_filtering<Polyhedron>
      >
class Surface_mesh_segmentation
{
//type definitions
public:
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor      face_descriptor;

private:
  //typedef typename Polyhedron::Traits Kernel;
  typedef typename GeomTraits::Point_3 Point;

  typedef typename boost::graph_traits<Polyhedron>::edge_iterator     edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator    face_iterator;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;

  typedef SDF_calculation<Polyhedron, VertexPointPmap, GeomTraits, fast_bbox_intersection>
  SDF_calculation_class;

// member variables
private:
  const Polyhedron& mesh;
  GeomTraits traits;
  VertexPointPmap vertex_point_pmap;
// member functions
public:
  /**
   * @pre @a polyhedron.is_pure_triangle()
   * @param mesh `CGAL Polyhedron` on which other functions operate.
   */
  Surface_mesh_segmentation(const Polyhedron& mesh, GeomTraits traits, VertexPointPmap vertex_point_pmap)
    : mesh(mesh), traits(traits), vertex_point_pmap(vertex_point_pmap) {
    //    CGAL_precondition(is_pure_triangle(mesh));
  }

// Use these two functions together
  template <class SDFPropertyMap>
  std::pair<double, double>
  calculate_sdf_values(double cone_angle, std::size_t number_of_rays,
                       SDFPropertyMap sdf_pmap, bool postprocess_req) {
    // calculate sdf values
    SDF_calculation_class sdf_calculator(mesh, vertex_point_pmap,
                                         false, /* build_kd_ree */
                                         true, /* use_diagonal --> set to false to use `AABB_tree::first_intersection()` */
                                         traits);
    sdf_calculator.calculate_sdf_values(faces(mesh).first, faces(mesh).second,
                                        cone_angle, number_of_rays, sdf_pmap);

    Postprocess_sdf_values<Polyhedron> p;
    return postprocess_req ? p.template postprocess<Filter>(mesh,
           sdf_pmap) : // return minimum and maximum sdf values before normalization
           p.min_max_value(mesh, sdf_pmap);
  }

  template <class FacetSegmentMap, class SDFPropertyMap>
  std::size_t partition(std::size_t number_of_centers, double smoothing_lambda,
                        SDFPropertyMap sdf_pmap, FacetSegmentMap segment_pmap,
                        bool clusters_to_segments) {
    smoothing_lambda = (std::max)(0.0, smoothing_lambda); // min zero
    smoothing_lambda *=
      CGAL_SMOOTHING_LAMBDA_MULTIPLIER; // scale it into meaningful range for graph-cut

    // log normalize sdf values
    std::vector<double> sdf_values;
    log_normalize_sdf_values(sdf_pmap, sdf_values);

    // soft clustering using GMM-fitting initialized with k-means
    Expectation_maximization fitter(number_of_centers, sdf_values,
                                    Expectation_maximization::K_MEANS_INITIALIZATION, 1);

    std::vector<std::size_t> labels;
    fitter.fill_with_center_ids(labels);

    std::vector<std::vector<double> > probability_matrix;
    fitter.fill_with_probabilities(probability_matrix);
    log_normalize_probability_matrix(probability_matrix);

    // calculating edge weights
    std::vector<std::pair<std::size_t, std::size_t> > edges;
    std::vector<double> edge_weights;
    calculate_and_log_normalize_dihedral_angles(smoothing_lambda, edges,
        edge_weights);

    // apply graph cut
    GraphCut()(edges, edge_weights, probability_matrix, labels);
    std::vector<std::size_t>::iterator label_it = labels.begin();
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend;
        ++facet_it, ++label_it) {
      put(segment_pmap, *facet_it, *label_it); // fill with cluster-ids
    }
    if(clusters_to_segments) {
      // assign a segment id for each facet
      std::size_t number_of_segments = assign_segments(number_of_centers, sdf_pmap,
                                       segment_pmap);
      return number_of_segments;
    }
    return number_of_centers;
  }

private:

  /**
   * Calculates dihedral angle between facets and normalize them between [0-1] from [0 - 2*pi].
   * Also convex dihedral angles are multiplied by a factor smaller than 1.0 which reduces their effect in graph-cut.
   * @pre parameter @a edge should not be a border.
   * @param edge whose dihedral angle is computed using incident facets
   * @return computed dihedral angle
   */
  double calculate_dihedral_angle_of_edge(halfedge_descriptor edge) const {

    CGAL_precondition( (! (face(edge,mesh)==boost::graph_traits<Polyhedron>::null_face()))
                       && (! (face(opposite(edge,mesh),mesh)==boost::graph_traits<Polyhedron>::null_face())) );
    const Point a = get(vertex_point_pmap,target(edge,mesh));
    const Point b = get(vertex_point_pmap,target(prev(edge,mesh),mesh));
    const Point c = get(vertex_point_pmap,target(next(edge,mesh),mesh));
    const Point d = get(vertex_point_pmap,target(next(opposite(edge,mesh),mesh),mesh));
    // As far as I check: if, say, dihedral angle is 5, this returns 175,
    // if dihedral angle is -5, this returns -175.
    // Another words this function returns angle between planes.
    double n_angle = to_double( ::CGAL::approximate_dihedral_angle(a, b, c, d) );
    n_angle /= 180.0;
    bool concave = n_angle > 0;
    double angle = 1 + ((concave ? -1 : +1) * n_angle);

    if(!concave) {
      angle *= CGAL_CONVEX_FACTOR;
    }
    return angle;
  }

  /**
   * Normalize sdf values using function:
   * normalized_sdf = log( alpha * ( current_sdf - min_sdf ) / ( max_sdf - min_sdf ) + 1 ) / log( alpha + 1 )
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::face_descriptor` as key and `double` as value type
   * @param[out] normalized_sdf_values normalized values stored in facet iteration order
   * Important note: @a sdf_values parameter should contain linearly normalized values between [0-1]
   */
  template<class SDFPropertyMap>
  void log_normalize_sdf_values(SDFPropertyMap sdf_values,
                                std::vector<double>& normalized_sdf_values) {
    normalized_sdf_values.reserve(num_faces(mesh));
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      double log_normalized = log(get(sdf_values, *facet_it) * CGAL_NORMALIZATION_ALPHA +
                                  1) / log(CGAL_NORMALIZATION_ALPHA + 1);
      normalized_sdf_values.push_back(log_normalized);
    }
  }

  /**
   * Receives probability-matrix with probabilities between [0-1], and returns log-normalized probabilities
   * which are suitable to use in graph-cut.
   * @param[in, out] probabilities probability matrix in [center][facet] order
   */
  void log_normalize_probability_matrix(std::vector<std::vector<double> >&
                                        probabilities) const {
    const double epsilon = 5e-6;
    for(std::vector<std::vector<double> >::iterator it_i = probabilities.begin();
        it_i != probabilities.end(); ++it_i) {
      for(std::vector<double>::iterator it = it_i->begin(); it != it_i->end(); ++it) {
        double probability = (std::max)(*it,
                                        epsilon); // give every facet a little probability to be in any cluster
        probability = -log(probability);
        *it = (std::max)(probability, std::numeric_limits<double>::epsilon());
        // zero values are not accepted in max-flow as weights for edges which connects some vertex with Source or Sink (in boost::boykov..)
      }
    }
  }

  /**
   * Calculates dihedral-angle based weight for each edge which is not a border edge.
   * @param smoothing_lambda a factor for each weight (weight *= smoothing_lambda).
   * @param[out] edges list of pair of neighbor facet ids
   * @param[out] edge_weights calculated weight for each edge in @a edges
   */
  void calculate_and_log_normalize_dihedral_angles(double smoothing_lambda,
      std::vector<std::pair<std::size_t, std::size_t> >& a_edges,
      std::vector<double>& edge_weights) const {
    // associate each facet with an id
    // important note: ids should be compatible with iteration order of facets:
    // [0 <- facet_begin(),...., size_of_facets() -1 <- facet_end()]
    // Why ? it is send to graph cut algorithm where other data associated with facets are also sorted according to iteration order.
    std::map<face_descriptor, std::size_t> facet_index_map;
    std::size_t facet_index = 0;
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend;
        ++facet_it, ++facet_index) {
      facet_index_map[*facet_it] = facet_index;
    }

    const double epsilon = 5e-6;
    // edges and their weights. pair<std::size_t, std::size_t> stores facet-id pairs (see above) (may be using boost::tuple can be more suitable)
    edge_iterator edge_it, eend;
    for(boost::tie(edge_it,eend) = edges(mesh);
        edge_it != eend; ++edge_it) {
      halfedge_descriptor hd = halfedge(*edge_it,mesh);
      halfedge_descriptor ohd = opposite(hd,mesh);
      if((face(hd,mesh)==boost::graph_traits<Polyhedron>::null_face())
         || (face(ohd,mesh)==boost::graph_traits<Polyhedron>::null_face()))  {
        continue;  // if edge does not contain two neighbor facets then do not include it in graph-cut
      }
      const std::size_t index_f1 = facet_index_map[face(hd,mesh)];
      const std::size_t index_f2 = facet_index_map[face(ohd,mesh)];
      a_edges.push_back(std::make_pair(index_f1, index_f2));

      double angle = calculate_dihedral_angle_of_edge(hd);

      angle = (std::max)(angle, epsilon);
      angle = -log(angle);
      angle *= smoothing_lambda;

      edge_weights.push_back(angle);
    }
  }

  template<class Pair>
  struct Sort_pairs_with_second {
    bool operator() (const Pair& pair_1, const Pair& pair_2) const {
      return pair_1.second < pair_2.second;
    }
  };

/////////////////////////////////////
//  0 0 1 1 0 0      0 0 2 2 4 4   //
//  0 0 1 1 0 0      0 0 2 2 4 4   //
//  1 0 2 2 1 1  ->  1 0 3 3 5 5   //
//  1 0 2 2 1 1      1 0 3 3 5 5   //
//  cluster-ids  ->  segment-ids   //
  /**
   * Definitions:
   *   -Cluster is a set of facet which can be connected or disconnected.
   *   -Segment is a connected set of facets which are placed under same cluster (after graph-cut).
   * Function takes a map which contains a cluster-id per facet. It then fills the map with segment-ids by giving a unique id to each
   * set of connected facets which are placed under same cluster. Note that returned segment-ids are ordered by average sdf value of segment ascen.
   *
   * @param number_of_clusters cluster-ids in @a segments should be between [0, number_of_clusters -1]
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::face_descriptor` as key and `double` as value type
   * @param[in, out] segments `ReadWritePropertyMap` with `Polyhedron::face_descriptor` as key and `std::size_t` as value type.
   * @return number of segments
   */
  template<class SegmentPropertyMap, class SDFProperyMap>
  std::size_t assign_segments(std::size_t number_of_clusters,
                              SDFProperyMap sdf_values, SegmentPropertyMap segments) {
    // assign a segment-id to each facet
    std::size_t segment_id = number_of_clusters;
    std::vector<std::pair<std::size_t, double> > segments_with_average_sdf_values;
    face_iterator facet_it, fend;
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      if(get(segments, *facet_it) <
          number_of_clusters) { // not visited by depth_first_traversal
        double average_sdf_value = breadth_first_traversal(*facet_it, segment_id,
                                   sdf_values, segments);

        segments_with_average_sdf_values.push_back(std::make_pair(segment_id,
            average_sdf_value));
        ++segment_id;
      }
    }
    // sort segments according to their average sdf value
    sort(segments_with_average_sdf_values.begin(),
         segments_with_average_sdf_values.end(),
         Sort_pairs_with_second<std::pair<std::size_t, double> >());
    // map each segment-id to its new sorted index
    std::vector<std::size_t> segment_id_to_sorted_id_map(
      segments_with_average_sdf_values.size());
    for(std::size_t index = 0; index < segments_with_average_sdf_values.size();
        ++index) {
      std::size_t segment_id = segments_with_average_sdf_values[index].first -
                               number_of_clusters;
      segment_id_to_sorted_id_map[segment_id] = index;
    }
    // make one-pass on facets. First make segment-id zero based by subtracting number_of_clusters
    //                        . Then place its sorted index to pmap
    for(boost::tie(facet_it,fend) = faces(mesh);
        facet_it != fend; ++facet_it) {
      std::size_t segment_id = get(segments, *facet_it) - number_of_clusters;
      put(segments, *facet_it, segment_id_to_sorted_id_map[segment_id]);
    }
    return segment_id - number_of_clusters;
  }

  /**
   * Breadth-first traverse all connected facets which has same cluster with @a facet.
   * Each visited facet assigned to @a segment_id.
   * @param facet root facet
   * @param segment_id segment-id of root facet
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::face_descriptor` as key and `double` as value type
   * @param[in, out] segments `ReadWritePropertyMap` with `Polyhedron::face_descriptor` as key and `std::size_t` as value type.
   * @return average sdf value for segment
   */
  template<class SegmentPropertyMap, class SDFProperyMap>
  double
  breadth_first_traversal(face_descriptor root, std::size_t segment_id,
                          SDFProperyMap sdf_values, SegmentPropertyMap segments) {
    std::queue<face_descriptor> facet_queue;
    facet_queue.push(root);

    std::size_t prev_segment_id = get(segments, root);
    put(segments, root, segment_id);

    double total_sdf_value = get(sdf_values, root);
    std::size_t    visited_facet_count = 1;

    while(!facet_queue.empty()) {
      face_descriptor facet = facet_queue.front();

      Halfedge_around_face_circulator<Polyhedron> facet_circulator(halfedge(facet,mesh),mesh), done(facet_circulator);
      do {
        if(face(opposite(*facet_circulator,mesh),mesh)==boost::graph_traits<Polyhedron>::null_face()) {
          continue;  // no facet to traversal
        }
        face_descriptor neighbor = face(opposite(*facet_circulator,mesh),mesh);
        if(prev_segment_id == get(segments, neighbor)) {
          put(segments, neighbor, segment_id);
          facet_queue.push(neighbor);

          total_sdf_value += get(sdf_values, neighbor);
          ++visited_facet_count;
        }
      } while( ++facet_circulator != done);

      facet_queue.pop();
    }

    return total_sdf_value / visited_facet_count;
  }

};
}//namespace internal
/// @endcond
} //namespace CGAL

#undef CGAL_NORMALIZATION_ALPHA
#undef CGAL_CONVEX_FACTOR
#undef CGAL_SMOOTHING_LAMBDA_MULTIPLIER
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H
