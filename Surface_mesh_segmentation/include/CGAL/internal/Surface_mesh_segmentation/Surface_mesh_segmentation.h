#ifndef CGAL_SURFACE_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_H
/* NEED TO BE DONE
 * About implementation:
 * +) I am not using BGL, as far as I checked there is a progress on BGL redesign
 *    (https://cgal.geometryfactory.com/CGAL/Members/wiki/Features/BGL) which introduces some features
 *    for face-based traversal / manipulation by FaceGraphs
 * +) Deciding on which parameters will be taken from user
 *
 * About paper (and correctness / efficiency etc.):
 * +) Weighting ray distances with inverse of their angles: not sure how to weight exactly
 * +) Anisotropic smoothing: have no idea what it is exactly, should read some material (google search is not enough)
 * +) Deciding how to generate rays in cone: for now using "polar angle" and "accept-reject (square)" and "concentric mapping" techniques
 */


#include <CGAL/internal/Surface_mesh_segmentation/Expectation_maximization.h>
#include <CGAL/internal/Surface_mesh_segmentation/K_means_clustering.h>
#include <CGAL/internal/Surface_mesh_segmentation/Filters.h>
//#include "Alpha_expansion_graph_cut.h"
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

#include <CGAL/internal/Surface_mesh_segmentation/SDF_calculation.h>

#include <CGAL/utility.h>
#include <CGAL/Timer.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include <boost/optional.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>
#include <map>

#define CGAL_NORMALIZATION_ALPHA 5.0
#define CGAL_CONVEX_FACTOR 0.08
#define CGAL_SMOOTHING_LAMBDA_MULTIPLIER 100.0

//IOY: these are going to be removed at the end (no CGAL_ pref)
#define ACTIVATE_SEGMENTATION_DEBUG
#ifdef ACTIVATE_SEGMENTATION_DEBUG
#define SEG_DEBUG(x) x;
#else
#define SEG_DEBUG(x)
#endif
// If defined then profile function becomes active, and called from constructor.
//#define SEGMENTATION_PROFILE

namespace CGAL
{
namespace internal
{
/**
 * @brief Main entry point for mesh segmentation algorithm.
 *
 * It is a connector class which uses:
 *   -  SDF_calculation for calculating sdf values
 *   -  Expectation_maximization for soft clustering
 *   -  An implementation of alpha-expansion graph cut for hard clustering (Alpha_expansion_graph_cut.h)
 *
 * Other than being a connector, it is also responsable for preprocess and postprocess on intermadiate data, which are:
 *   - log-normalizing probabilities received from soft clustering
 *   - log-normalizing and calculating dihedral-angle based weights for edges
 *   - smoothing and log-normalizing sdf values received from sdf calculation (Filters.h)
 *   - assigning segment-id for each facet after hard clustering
 */
template <
class Polyhedron,
      class GraphCut = Alpha_expansion_graph_cut_boost,
      class FacetIndexMap =
      boost::associative_property_map<std::map<typename Polyhedron::Facet_const_handle, int> >,
      class Filter = Bilateral_filtering<Polyhedron>
      >
class Surface_mesh_segmentation
{
private:
  /**
   * An adaptor for Lvalue property-map. It stores a pointer to vector for underlying data-structure,
   * and also stores another property-map which maps each `key` to vector index.
   */
  template<class AnyPolyhedron, class ValueType, class FacetIdPropertyMap>
  struct Polyhedron_property_map_for_facet
      : public boost::put_get_helper<ValueType&,
        Polyhedron_property_map_for_facet<AnyPolyhedron, ValueType, FacetIdPropertyMap>
        > {
  public:
    typedef typename AnyPolyhedron::Facet_const_handle key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::writable_property_map_tag category;

    Polyhedron_property_map_for_facet() : internal_vector(NULL) { }
    Polyhedron_property_map_for_facet(std::vector<ValueType>* internal_vector,
                                      FacetIdPropertyMap facet_id_map)
      : internal_vector(internal_vector), facet_id_map(facet_id_map)  { }

    reference operator[](key_type key) const {
      return (*internal_vector)[facet_id_map[key]];
    }
  private:
    std::vector<ValueType>* internal_vector;
    FacetIdPropertyMap facet_id_map;
  };

//type definitions
public:
  typedef typename Polyhedron::Facet_const_handle      Facet_const_handle;
private:
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Polyhedron::Facet  Vertex;
  typedef typename Kernel::Vector_3   Vector;
  typedef typename Kernel::Point_3    Point;
  typedef typename Kernel::Plane_3   Plane;

  typedef typename Polyhedron::Edge_const_iterator     Edge_const_iterator;
  typedef typename Polyhedron::Halfedge_const_handle   Edge_const_handle;
  typedef typename Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Polyhedron::Facet_const_iterator    Facet_const_iterator;
  typedef typename Polyhedron::Vertex_const_iterator   Vertex_const_iterator;

  typedef SDF_calculation<Polyhedron> SDF_calculation_class;

// member variables
private:
  const Polyhedron& mesh;

// member functions
public:
  /**
   * @pre parameter @a mesh should consist of triangles.
   * @param mesh `CGAL Polyhedron` on which other functions operate.
   */
  Surface_mesh_segmentation(const Polyhedron& mesh)
    : mesh(mesh) {
    CGAL_precondition(mesh.is_pure_triangle());
  }

// Use these two functions together
  template <class SDFPropertyMap>
  std::pair<double, double>
  calculate_sdf_values(SDFPropertyMap sdf_pmap, double cone_angle,
                       int number_of_rays) {
    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())

    SDF_calculation_class().calculate_sdf_values(mesh, sdf_pmap, cone_angle,
        number_of_rays);

    SEG_DEBUG(std::cout << t.time() << std::endl)

    check_zero_sdf_values(sdf_pmap);
    Filter()(mesh, get_window_size(), sdf_pmap);
    std::pair<double, double> min_max_sdf_values = linear_normalize_sdf_values(
          sdf_pmap);

    SEG_DEBUG(std::cout << t.time() << std::endl)
    return min_max_sdf_values;
  }

  template <class FacetSegmentMap, class SDFPropertyMap>
  int partition(SDFPropertyMap sdf_pmap, FacetSegmentMap segment_pmap,
                int number_of_centers,
                double smoothing_lambda) {
    smoothing_lambda = (std::max)(0.0, (std::min)(1.0,
                                  smoothing_lambda)); // clip into [0-1]
    smoothing_lambda *=
      CGAL_SMOOTHING_LAMBDA_MULTIPLIER; // scale it into meaningful range for graph-cut

    // log normalize sdf values
    std::vector<double> sdf_values;
    log_normalize_sdf_values(sdf_pmap, sdf_values);
    // soft clustering using GMM-fitting initialized with k-means
    Expectation_maximization fitter(number_of_centers, sdf_values,
                                    Expectation_maximization::K_MEANS_INITIALIZATION, 1);

    std::vector<int> labels;
    fitter.fill_with_center_ids(labels);

    std::vector<std::vector<double> > probability_matrix;
    fitter.fill_with_probabilities(probability_matrix);
    log_normalize_probability_matrix(probability_matrix);

    // calculating edge weights
    std::vector<std::pair<int, int> > edges;
    std::vector<double> edge_weights;
    calculate_and_log_normalize_dihedral_angles(smoothing_lambda, edges,
        edge_weights);

    // apply graph cut
    GraphCut()(edges, edge_weights, probability_matrix, labels);
    std::vector<int>::iterator label_it = labels.begin();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++label_it) {
      segment_pmap[facet_it] = *label_it;
    }
    // assign a segment id for each facet
    //return number_of_centers;
    int number_of_segments = assign_segments(number_of_centers, sdf_pmap,
                             segment_pmap);
    //std::cout << "ne : " <<number_of_segments << std::endl;
    return number_of_segments;
  }

private:
  /**
   * Going to be removed
   */
  double calculate_dihedral_angle_of_edge(Edge_const_handle& edge) const {
    CGAL_precondition(!edge->is_border_edge());
    Facet_const_handle f1 = edge->facet();
    Facet_const_handle f2 = edge->opposite()->facet();

    const Point& f2_v1 = f2->halfedge()->vertex()->point();
    const Point& f2_v2 = f2->halfedge()->next()->vertex()->point();
    const Point& f2_v3 = f2->halfedge()->prev()->vertex()->point();
    /*
     * As far as I see from results, segment boundaries are occurred in 'concave valleys'.
     * There is no such thing written (clearly) in the paper but should we just penalize 'concave' edges (not convex edges) ?
     * Actually that is what I understood from 'positive dihedral angle'.
     */
    const Point& unshared_point_on_f1 = edge->next()->vertex()->point();
    Plane p2(f2_v1, f2_v2, f2_v3);
    bool concave = p2.has_on_positive_side(unshared_point_on_f1);

    const Point& f1_v1 = f1->halfedge()->vertex()->point();
    const Point& f1_v2 = f1->halfedge()->next()->vertex()->point();
    const Point& f1_v3 = f1->halfedge()->prev()->vertex()->point();
    Vector f1_normal = unit_normal(f1_v1, f1_v2, f1_v3);
    Vector f2_normal = unit_normal(f2_v1, f2_v2, f2_v3);

    double dot = f1_normal * f2_normal;
    if(dot > 1.0)       {
      dot = 1.0;
    } else if(dot < -1.0) {
      dot = -1.0;
    }
    double angle = acos(dot) / CGAL_PI; // [0-1] normalize

    if(!concave) {
      angle *= CGAL_CONVEX_FACTOR;
    }
    return angle;
  }

  /**
   * Calculates dihedral angle between facets and normalize them between [0-1] from [0 - 2*pi].
   * Also convex dihedral angles are multiplied by a factor smaller than 1.0 which reduces their effect in graph-cut.
   * @pre parameter @a edge should not be a border.
   * @param edge whose dihedral angle is computed using incident facets
   * @return computed dihedral angle
   */
  double calculate_dihedral_angle_of_edge_2(Edge_const_handle& edge) const {
    CGAL_precondition(!edge->is_border_edge());
    const Point& a = edge->vertex()->point();
    const Point& b = edge->prev()->vertex()->point();
    const Point& c = edge->next()->vertex()->point();
    const Point& d = edge->opposite()->next()->vertex()->point();
    // As far as I check: if, say, dihedral angle is 5, this returns 175,
    // if dihedral angle is -5, this returns -175.
    // Another words this function returns angle between planes.
    double n_angle = Mesh_3::dihedral_angle(a, b, c, d) / 180.0;
    bool concave = n_angle > 0;
    double angle = 1 + ((concave ? -1 : +1) * n_angle);

    if(!concave) {
      angle *= CGAL_CONVEX_FACTOR;
    }
    return angle;

    //Facet_const_handle f1 = edge->facet();
    //Facet_const_handle f2 = edge->opposite()->facet();
    //
    //const Point& f2_v1 = f2->halfedge()->vertex()->point();
    //const Point& f2_v2 = f2->halfedge()->next()->vertex()->point();
    //const Point& f2_v3 = f2->halfedge()->prev()->vertex()->point();
    ///*
    // * As far as I see from results, segment boundaries are occurred in 'concave valleys'.
    // * There is no such thing written (clearly) in the paper but should we just penalize 'concave' edges (not convex edges) ?
    // * Actually that is what I understood from 'positive dihedral angle'.
    // */
    //const Point& unshared_point_on_f1 = edge->next()->vertex()->point();
    //Plane p2(f2_v1, f2_v2, f2_v3);
    //bool concave = p2.has_on_positive_side(unshared_point_on_f1);
    ////std::cout << n_angle << std::endl;
    ////if(!concave) { return epsilon; } // So no penalties for convex dihedral angle ? Not sure...

    //const Point& f1_v1 = f1->halfedge()->vertex()->point();
    //const Point& f1_v2 = f1->halfedge()->next()->vertex()->point();
    //const Point& f1_v3 = f1->halfedge()->prev()->vertex()->point();
    //Vector f1_normal = unit_normal(f1_v1, f1_v2, f1_v3);
    //Vector f2_normal = unit_normal(f2_v1, f2_v2, f2_v3);
    //
    //double dot = f1_normal * f2_normal;
    //if(dot > 1.0)       { dot = 1.0;  }
    //else if(dot < -1.0) { dot = -1.0; }
    //double angle = acos(dot) / CGAL_PI; // [0-1] normalize
    //std::cout << angle << " " << n_angle << " " << (concave ? "concave": "convex") << std::endl;
    //if(std::abs(angle - folded_angle) > 1e-6)
    //{
    //
    //}
    //if(angle < epsilon) { angle = epsilon; }
    //return angle;
  }

  /**
   * Normalize sdf values using function:
   * normalized_sdf = log( alpha * ( current_sdf - min_sdf ) / ( max_sdf - min_sdf ) + 1 ) / log( alpha + 1 )
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   * @param[out] normalized_sdf_values normalized values stored in facet iteration order
   */
  template<class SDFPropertyMap>
  void log_normalize_sdf_values(SDFPropertyMap sdf_values,
                                std::vector<double>& normalized_sdf_values) {
    normalized_sdf_values.reserve(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      double log_normalized = log(sdf_values[facet_it] * CGAL_NORMALIZATION_ALPHA +
                                  1) / log(CGAL_NORMALIZATION_ALPHA + 1);
      normalized_sdf_values.push_back(log_normalized);
    }
  }

  /**
   * Normalize sdf values between [0-1].
   */
  template<class SDFPropertyMap>
  std::pair<double, double> linear_normalize_sdf_values(SDFPropertyMap
      sdf_values) {
    double max_sdf = -(std::numeric_limits<double>::max)();
    double min_sdf = (std::numeric_limits<double>::max)();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      double sdf_value = sdf_values[facet_it];
      max_sdf = (std::max)(sdf_value, max_sdf);
      min_sdf = (std::min)(sdf_value, min_sdf);
    }
    const double max_min_dif = max_sdf - min_sdf;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      sdf_values[facet_it] = (sdf_values[facet_it] - min_sdf) / max_min_dif;
    }
    return std::pair<double, double>(min_sdf, max_sdf);
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
  int get_window_size() {

    double facet_sqrt = std::sqrt(mesh.size_of_facets() / 2000.0);
    return static_cast<int>(facet_sqrt) + 1;
  }

  /**
   * Finds facets which have zero sdf values.
   * Sdf values on these facets are assigned to average sdf value of its neighbors.
   * If still there is any facet which has no sdf value, assigns minimum sdf value to it.
   * This is meaningful since (being an outlier) zero sdf values might effect normalization & log extremely.
   * @param[in, out] sdf_values `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   */
  template<class SDFPropertyMap>
  void check_zero_sdf_values(SDFPropertyMap sdf_values) {
    std::vector<Facet_const_handle> still_zero_facets;
    double min_sdf = (std::numeric_limits<double>::max)();
    // If there is any facet which has no sdf value, assign average sdf value of its neighbors
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      double sdf_value = sdf_values[facet_it];
      if(sdf_value == 0.0) {
        typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
          facet_it->facet_begin();
        double total_neighbor_sdf = 0.0;
        do {
          if(!facet_circulator->opposite()->is_border()) {
            total_neighbor_sdf += sdf_values[facet_circulator->opposite()->facet()];
          }
        } while( ++facet_circulator !=  facet_it->facet_begin());

        sdf_value = total_neighbor_sdf / 3.0;
        sdf_values[facet_it] = sdf_value;

        if(sdf_value == 0.0) {
          still_zero_facets.push_back(
            facet_it);     // if sdf_value is still zero, put facet to zero-list
        } else                 {
          min_sdf = (std::min)(sdf_value,
                               min_sdf);  // find min_sdf other than zero
        }
      } else {
        min_sdf = (std::min)(sdf_value, min_sdf);  // find min_sdf other than zero
      }
    }
    // If still there is any facet which has no sdf value, assign minimum sdf value.
    // This is meaningful since (being an outlier) 0 sdf values might effect normalization & log extremely.
    for(typename std::vector<Facet_const_handle>::iterator it =
          still_zero_facets.begin();
        it != still_zero_facets.end(); ++it) {
      sdf_values[*it] = min_sdf;
    }
  }

  /**
   * Receives probability-matrix with probabilities betwen [0-1], and returns log-normalized probabilities
   * which are suitable to use in graph-cut.
   * @param[in, out] probabilities probability matrix in [center][facet] order
   */
  void log_normalize_probability_matrix(std::vector<std::vector<double> >&
                                        probabilities) const {
    const double epsilon = 1e-5;
    for(std::vector<std::vector<double> >::iterator it_i = probabilities.begin();
        it_i != probabilities.end(); ++it_i) {
      for(std::vector<double>::iterator it = it_i->begin(); it != it_i->end(); ++it) {
        double probability = (std::max)(*it,
                                        epsilon); // give any facet a little probability to be in any cluster
        probability = -log(probability);
        *it = (std::max)(probability,
                         std::numeric_limits<double>::epsilon()); // zero values are not accepted in max-flow
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
      std::vector<std::pair<int, int> >& edges,
      std::vector<double>& edge_weights) const {
    std::map<Facet_const_handle, int> facet_index_map;
    int facet_index = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++facet_index) {
      facet_index_map[facet_it] = facet_index;
    }

    const double epsilon = 1e-5;
    //edges and their weights. pair<int, int> stores facet-id pairs (see above) (may be using boost::tuple can be more suitable)
    for(Edge_const_iterator edge_it = mesh.edges_begin();
        edge_it != mesh.edges_end(); ++edge_it) {
      if(edge_it->is_border_edge()) {
        continue;  // if edge does not contain two neighbor facets then do not include it in graph-cut
      }
      const int index_f1 = facet_index_map[edge_it->facet()];
      const int index_f2 = facet_index_map[edge_it->opposite()->facet()];
      edges.push_back(std::pair<int, int>(index_f1, index_f2));

      double angle = calculate_dihedral_angle_of_edge_2(edge_it);

      if(angle < epsilon) {
        angle = epsilon;
      }
      angle = -log(angle);
      angle = (std::max)(angle, std::numeric_limits<double>::epsilon());
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
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   * @param[in, out] segments `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type.
   * @return number of segments
   */
  template<class SegmentPropertyMap, class SDFProperyMap>
  int assign_segments(int number_of_clusters, SDFProperyMap sdf_values,
                      SegmentPropertyMap segments) {
    int segment_id = number_of_clusters;
    std::vector<std::pair<int, double> > segments_with_average_sdf_values;

    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      if(segments[facet_it] <
          number_of_clusters) { // not visited by depth_first_traversal
        std::pair<double, int> sdf_facet_count_pair = depth_first_traversal(facet_it,
            segment_id, sdf_values, segments);
        double average_sdf_value_for_segment = sdf_facet_count_pair.first /
                                               sdf_facet_count_pair.second;
        segments_with_average_sdf_values.push_back(std::pair<int, double>(segment_id,
            average_sdf_value_for_segment));
        ++segment_id;
      }
    }
    // sort segments according to their average sdf value
    sort(segments_with_average_sdf_values.begin(),
         segments_with_average_sdf_values.end(),
         Sort_pairs_with_second<std::pair<int, double> >());
    // map each segment-id to its new sorted index
    std::vector<int> segment_id_to_sorted_id_map(
      segments_with_average_sdf_values.size());
    for(std::size_t index = 0; index < segments_with_average_sdf_values.size();
        ++index) {
      int segment_id = segments_with_average_sdf_values[index].first -
                       number_of_clusters;
      segment_id_to_sorted_id_map[segment_id] = index;
    }
    // make one-pass on facets. First make segment-id zero based by subtracting number_of_clusters
    //                        . Then place its sorted index to pmap
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      int segment_id = segments[facet_it] - number_of_clusters;
      segments[facet_it] = segment_id_to_sorted_id_map[segment_id];
    }
    return segment_id - number_of_clusters;
  }

  /**
   * Depth-first traverse all connected facets which has same cluster with @a facet.
   * Each visited facet assigned to @a segment_id.
   * @param facet root facet
   * @param segment_id segment-id of root facet
   * @param sdf_values `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   * @param[in, out] segments `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type.
   * @return pair of first: accumulated sdf values of visited facets, second: number of visited facets
   */
  template<class SegmentPropertyMap, class SDFProperyMap>
  std::pair<double, int>
  depth_first_traversal(Facet_const_handle& facet, int segment_id,
                        SDFProperyMap sdf_values, SegmentPropertyMap segments) {
    int prev_segment_id = segments[facet];
    segments[facet] = segment_id;
    std::pair<double, int> sdf_facet_count_pair(sdf_values[facet], 1);

    typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
      facet->facet_begin();
    do {
      if(facet_circulator->opposite()->is_border()) {
        continue;  // no facet to traversal
      }
      Facet_const_handle neighbor = facet_circulator->opposite()->facet();
      if(prev_segment_id == segments[neighbor]) {
        const std::pair<double, int>& total_pair = depth_first_traversal(neighbor,
            segment_id, sdf_values, segments);
        sdf_facet_count_pair.first += total_pair.first;
        sdf_facet_count_pair.second += total_pair.second;
      }
    } while( ++facet_circulator !=  facet->facet_begin());

    return sdf_facet_count_pair;
  }

};
}//namespace internal
} //namespace CGAL

#undef CGAL_NORMALIZATION_ALPHA
#undef CGAL_CONVEX_FACTOR
#undef CGAL_SMOOTHING_LAMBDA_MULTIPLIER

#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H
