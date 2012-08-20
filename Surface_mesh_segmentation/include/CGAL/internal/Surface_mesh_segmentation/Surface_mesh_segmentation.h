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
//#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>
#include "Alpha_expansion_graph_cut.h"
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
//AF: macros must be prefixed with "CGAL_"
//IOY: done
#define CGAL_DEFAULT_NUMBER_OF_CLUSTERS 5
#define CGAL_DEFAULT_SMOOTHING_LAMBDA 0.23
#define CGAL_DEFAULT_CONE_ANGLE (2.0 / 3.0) * CGAL_PI
#define CGAL_DEFAULT_NUMBER_OF_RAYS 25

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
 * It is a connector class which uses:
 *   -  SDF_calculation for calculating sdf values
 *   -  Expectation_maximization for soft clustering
 *   -  Alpha_expansion_graph_cut for hard clustering
 *
 * Other than being a connector, it is also responsable for preprocess and postprocess on intermadiate data, which are:
 *   - log-normalizing probabilities received from soft clustering
 *   - log-normalizing and calculating dihedral-angle based weights for edges
 *   - smoothing and log-normalizing sdf values received from sdf calculation
 *   - assigning segment-id for each facet after hard clustering
 */
template <
class Polyhedron,
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
  template<class Any_polyhedron, class ValueType, class FacetIdPropertyMap>
  struct Polyhedron_property_map_for_facet
      : public boost::put_get_helper<ValueType&,
        Polyhedron_property_map_for_facet<Any_polyhedron, ValueType, FacetIdPropertyMap>
        > {
  public:
    typedef typename Any_polyhedron::Facet_const_handle key_type;
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
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Polyhedron::Facet  Vertex;
  typedef typename Kernel::Vector_3   Vector;
  typedef typename Kernel::Point_3    Point;

  typedef typename Polyhedron::Edge_const_iterator     Edge_const_iterator;
  typedef typename Polyhedron::Halfedge_const_handle   Edge_const_handle;
  typedef typename Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Polyhedron::Facet_const_iterator    Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle      Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_iterator   Vertex_const_iterator;

  typedef SDF_calculation<Polyhedron> SDF_calculation_class;

private:
  typedef typename Kernel::Plane_3   Plane;

// member variables
public: // going to be private !
  const Polyhedron& mesh;

  FacetIndexMap
  facet_index_map; /**< Property-map which associates each facet with an id */
  std::map<Facet_const_handle, int> facet_index_map_internal;

  std::vector<double> sdf_values;
  std::vector<int>    centers;
  std::vector<int>    segments;

  bool is_pmmap_custom;

// member functions
public:
  Surface_mesh_segmentation(const Polyhedron& mesh)
    : mesh(mesh), facet_index_map(facet_index_map_internal),
      is_pmmap_custom(false) {
    int facet_index = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++facet_index) {
      facet_index_map[facet_it] = facet_index;
    }
  }

  Surface_mesh_segmentation(const Polyhedron& mesh, FacetIndexMap facet_index_map)
    : mesh(mesh), facet_index_map(facet_index_map), is_pmmap_custom(true) {
  }

// Copy-constructor
// Default, we store facet-to-id in property-map adaptor 'FacetIndexMap'.
// It stores a pointer to std::map instance which is a member variable.
// So default copy constructor is not working since copied property-map adaptor still
// points to std::map instance which exists in copied object (rhs).
  Surface_mesh_segmentation(const
                            Surface_mesh_segmentation<Polyhedron, FacetIndexMap>& rhs):
    mesh(rhs.mesh), facet_index_map_internal(rhs.facet_index_map_internal),
    sdf_values(rhs.sdf_values),
    centers(rhs.centers), segments(rhs.segments),
    is_pmmap_custom(rhs.is_pmmap_custom) {
    if(!is_pmmap_custom) {
      facet_index_map = FacetIndexMap(facet_index_map_internal);
    }
  }

// Use another copy of polyhedron but also copy the state of rhs.
// Main difference between copy constructor is rhs.mesh is NOT equal to this->mesh.
// So we cannot copy facet_index_map_internal since facet_iterators of 'this' is different than rhs.
  Surface_mesh_segmentation(const Polyhedron& mesh,
                            const Surface_mesh_segmentation<Polyhedron, FacetIndexMap>& rhs)
    : mesh(mesh), sdf_values(rhs.sdf_values), centers(rhs.centers),
      segments(rhs.segments),
      is_pmmap_custom(rhs.is_pmmap_custom) {
    if(!is_pmmap_custom) {
      facet_index_map = FacetIndexMap(facet_index_map_internal);

      int facet_index = 0;
      for(Facet_const_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end();
          ++facet_it, ++facet_index) {
        facet_index_map[facet_it] = facet_index;
      }
    }
  }
///////////////////////////////////////////////////////////////////
// Use these two functions together
  void calculate_sdf_values(double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                            int number_of_ray = CGAL_DEFAULT_NUMBER_OF_RAYS) {
    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())
    //read_sdf_values("D:/dino.txt");
    //return;

    sdf_values = std::vector<double>(mesh.size_of_facets());
    Polyhedron_property_map_for_facet<Polyhedron, double, FacetIndexMap>
    sdf_pmap(&sdf_values, facet_index_map);

    SDF_calculation_class().calculate_sdf_values(mesh, sdf_pmap, cone_angle,
        number_of_ray);

    SEG_DEBUG(std::cerr << "SDF computation time: " << t.time() << std::endl)
    SEG_DEBUG(t.reset())

    check_zero_sdf_values();
    smooth_sdf_values();
    normalize_sdf_values();

    SEG_DEBUG(std::cerr <<"Normalization and smoothing time: " << t.time() <<
              std::endl)
  }

  void partition(int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
                 double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA) {
    smoothing_lambda = (std::max)(0.0, (std::min)(1.0,
                                  smoothing_lambda)); // clip into [0-1]
    smoothing_lambda *=
      CGAL_SMOOTHING_LAMBDA_MULTIPLIER; // scale it into meaningful range for graph-cut

    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())
    // soft clustering using GMM-fitting initialized with k-means
    Expectation_maximization fitter(number_of_centers, sdf_values,
                                    Expectation_maximization::K_MEANS_INITIALIZATION, 1);

    std::vector<int> labels;
    fitter.fill_with_center_ids(labels);

    std::vector<std::vector<double> > probability_matrix;
    fitter.fill_with_probabilities(probability_matrix);
    log_normalize_probability_matrix(probability_matrix);

    SEG_DEBUG(std::cerr <<"Soft clustering time: " << t.time() << std::endl)
    t.reset();

    // calculating edge weights
    std::vector<std::pair<int, int> > edges;
    std::vector<double> edge_weights;
    calculate_and_log_normalize_dihedral_angles(smoothing_lambda, edges,
        edge_weights);

    // apply graph cut
    Alpha_expansion_graph_cut(edges, edge_weights, probability_matrix, labels);
    centers = labels;
    // assign a segment id for each facet
    assign_segments();
  }
///////////////////////////////////////////////////////////////////
// Use these two functions together
  template <class SDFPropertyMap>
  void calculate_sdf_values(SDFPropertyMap sdf_pmap,
                            double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                            int number_of_rays = CGAL_DEFAULT_NUMBER_OF_RAYS) {
    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())

    sdf_values = std::vector<double>(mesh.size_of_facets());
    Polyhedron_property_map_for_facet<Polyhedron, double, FacetIndexMap>
    sdf_pmap_internal(&sdf_values, facet_index_map);

    SDF_calculation_class().calculate_sdf_values(mesh, sdf_pmap_internal,
        cone_angle, number_of_rays);

    SEG_DEBUG(std::cout << t.time() << std::endl)

    check_zero_sdf_values();
    smooth_sdf_values();
    linear_normalize_sdf_values();

    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      sdf_pmap[facet_it] = get(sdf_values, facet_it);
    }

    SEG_DEBUG(std::cout << t.time() << std::endl)
  }

  template <class FacetSegmentMap, class SDFPropertyMap>
  void partition(SDFPropertyMap sdf_pmap, FacetSegmentMap segment_pmap,
                 int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
                 double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA) {
    smoothing_lambda = (std::max)(0.0, (std::min)(1.0,
                                  smoothing_lambda)); // clip into [0-1]
    smoothing_lambda *=
      CGAL_SMOOTHING_LAMBDA_MULTIPLIER; // scale it into meaningful range for graph-cut

    sdf_values = std::vector<double>(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      get(sdf_values, facet_it) = sdf_pmap[facet_it];
    }

    // log normalize sdf values
    normalize_sdf_values();
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
    Alpha_expansion_graph_cut(edges, edge_weights, probability_matrix, labels);
    centers = labels;
    // assign a segment id for each facet
    assign_segments();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      segment_pmap[facet_it] = get(segments, facet_it);
    }
  }
///////////////////////////////////////////////////////////////////
  double get_sdf_value_of_facet(Facet_const_handle facet) const {
    return sdf_values[facet_index_map[facet]];
  }

  int get_segment_id_of_facet(Facet_const_handle facet) const {
    return segments[facet_index_map[facet]];
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
   */
  void normalize_sdf_values() {
    typedef std::vector<double>::iterator fv_iterator;
    std::pair<fv_iterator, fv_iterator> min_max_pair =
      min_max_element(sdf_values.begin(), sdf_values.end());

    const double max_value = *min_max_pair.second;
    const double min_value = *min_max_pair.first;
    const double max_min_dif = max_value - min_value;
    for(fv_iterator it = sdf_values.begin(); it != sdf_values.end(); ++it) {
      double linear_normalized = (*it - min_value) / max_min_dif;
      double log_normalized = log(linear_normalized * CGAL_NORMALIZATION_ALPHA + 1) /
                              log(CGAL_NORMALIZATION_ALPHA + 1);
      *it = log_normalized;
    }
  }

  /**
   * Normalize sdf values between [0-1].
   */
  void linear_normalize_sdf_values() {
    typedef std::vector<double>::iterator fv_iterator;
    std::pair<fv_iterator, fv_iterator> min_max_pair =
      min_max_element(sdf_values.begin(), sdf_values.end());

    const double max_value = *min_max_pair.second;
    const double min_value = *min_max_pair.first;
    const double max_min_dif = max_value - min_value;
    for(fv_iterator it = sdf_values.begin(); it != sdf_values.end(); ++it) {
      *it = (*it - min_value) / max_min_dif;
    }
  }

  /**
   * Smoothing sdf values.
   * Takes neighbors using `get_window_size()`, and assigns weighted average of neighbors as filtered result.
   * For weighting two weights are multiplied:
   *   - spatial: over geodesic distances (number of edges)
   *   - domain : over sdf value distances
   */
  void smooth_sdf_values() {
    typedef Polyhedron_property_map_for_facet<Polyhedron, double, FacetIndexMap>
    Facet_vector_map;

    Facet_vector_map sdf_input(&sdf_values, facet_index_map);

    std::vector<double> smoothed_sdf(mesh.size_of_facets());
    Facet_vector_map sdf_output(&smoothed_sdf, facet_index_map);

    Filter()(mesh, get_window_size(), sdf_input, sdf_output);
    sdf_values = smoothed_sdf;
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
   */
  void check_zero_sdf_values() {
    // If there is any facet which has no sdf value, assign average sdf value of its neighbors
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      if(get(sdf_values, facet_it) == 0.0) {
        typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
          facet_it->facet_begin();
        double total_neighbor_sdf = 0.0;
        do {
          if(!facet_circulator->opposite()->is_border()) {
            total_neighbor_sdf += get(sdf_values, facet_circulator->opposite()->facet());
          }
        } while( ++facet_circulator !=  facet_it->facet_begin());
        get(sdf_values, facet_it) = total_neighbor_sdf / 3.0;
      }
    }
    // Find minimum sdf value other than 0
    double min_sdf = (std::numeric_limits<double>::max)();
    for(std::vector<double>::iterator it = sdf_values.begin();
        it != sdf_values.end(); ++it) {
      if(*it < min_sdf && *it != 0.0) {
        min_sdf = *it;
      }
    }
    // If still there is any facet which has no sdf value, assign minimum sdf value.
    // This is meaningful since (being an outlier) 0 sdf values might effect normalization & log extremely.
    for(std::vector<double>::iterator it = sdf_values.begin();
        it != sdf_values.end(); ++it) {
      if(*it == 0.0) {
        *it = min_sdf;
      }
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

  /**
   * Assigns a segment-id to each facet.
   * Segment is a connected set of facets which are placed under same cluster (after graph-cut).
   */
  void assign_segments() {
    segments = std::vector<int>(mesh.size_of_facets(), -1);
    int segment_id = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      if(get(segments, facet_it) == -1) { // not visited by depth_first_traversal
        depth_first_traversal(facet_it, segment_id++);
      }
    }
  }

  /**
   * Depth-first traverse all connected facets which has same cluster with @a facet.
   * Each visited facet assigned to @a segment_id.
   * @param facet root facet
   * @param segment_id segment-id of root facet
   */
  void depth_first_traversal(Facet_const_handle& facet, int segment_id) {
    get(segments, facet) = segment_id;
    typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
      facet->facet_begin();
    double total_neighbor_sdf = 0.0;
    do {
      if(facet_circulator->opposite()->is_border()) {
        continue;  // no facet to traversal
      }
      Facet_const_handle neighbor = facet_circulator->opposite()->facet();

      bool not_visited_before = get(segments, neighbor) == -1;
      bool under_same_cluster = get(centers, facet) == get(centers, neighbor);

      if(not_visited_before && under_same_cluster) {
        depth_first_traversal(neighbor, segment_id);
      }
    } while( ++facet_circulator !=  facet->facet_begin());
  }

  /**
   * Every assosiated data with facets are stored using vectors.
   * These data are reached using facet-ids which are received from #facet_index_map.
   * @param data associated data vector with facets
   * @param facet key to reach associated data
   */
  template<class T>
  T& get(std::vector<T>& data, const Facet_const_handle& facet) {
    return data[ facet_index_map[facet] ];
  }



public:
  /**
   * Going to be removed
   */
  void write_sdf_values(const char* file_name) {
    std::ofstream output(file_name);
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      output << get(sdf_values, facet_it) << std::endl;
    }
    output.close();
  }
  /**
   * Going to be removed
   */
  void read_sdf_values(const char* file_name) {
    std::ifstream input(file_name);
    sdf_values = std::vector<double>(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      double sdf_value;
      input >> sdf_value;
      get(sdf_values, facet_it) = sdf_value;
    }
  }
  /**
   * Going to be removed
   */
  void read_center_ids(const char* file_name) {
    //std::ifstream input(file_name);
    //centers.clear();
    //int max_center = 0;
    //for(Facet_const_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it)
    //{
    //    int center_id;
    //    input >> center_id;
    //    centers.insert(std::pair<Facet_const_handle, int>(facet_it, center_id));
    //    if(center_id > max_center) { max_center = center_id; }
    //}
    //number_of_centers = max_center + 1;
  }

  /**
   * Going to be removed
   */
  void read_probabilities(const char* file_name,
                          std::vector<std::vector<double> > & probability_matrix) {
    std::ifstream input(file_name);
    for(std::vector<std::vector<double> >::iterator vec_it =
          probability_matrix.begin(); vec_it != probability_matrix.end(); ++vec_it) {
      for(std::vector<double>::iterator it = vec_it->begin(); it != vec_it->end();
          ++it) {
        input >> (*it);
      }
    }
  }

  /**
   * Going to be removed
   */
  void write_segment_ids(const char* file_name) {
    assign_segments();
    std::ofstream output(file_name);
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      output << get(segments, facet_it) << std::endl;
    }
    output.close();
  }
  /**
   * Going to be removed
   */
  void profile(const char* file_name) {

#ifdef SEGMENTATION_PROFILE
    typedef Listing_intersection_traits_ray_or_segment_triangle
    < typename Tree::AABB_traits, Ray, std::back_insert_iterator< std::list<Object_and_primitive_id> > >
    Traits_with_ray;
    typedef Listing_intersection_traits_ray_or_segment_triangle
    < typename Tree::AABB_traits, Segment, std::back_insert_iterator< std::list<Object_and_primitive_id> > >
    Traits_with_segment;
    std::ofstream output(file_name);

    //for minimum
    for(int i = 0; i < 5; ++i) {
      use_minimum_segment = true;
      multiplier_for_segment = 1.0 + i;

      Timer t;
      t.start();
      calculate_sdf_values();

      double time = t.time();
      long inter_counter = Traits_with_ray::inter_counter +
                           Traits_with_segment::inter_counter;
      long true_inter_counter = Traits_with_ray::true_inter_counter +
                                Traits_with_segment::true_inter_counter;
      long do_inter_counter  = Traits_with_ray::do_inter_counter +
                               Traits_with_segment::do_inter_counter;
      long do_inter_in_segment = Traits_with_segment::do_inter_counter;
      long do_inter_in_ray = Traits_with_ray::do_inter_counter;

      output << "time: " << time << std::endl;
      output << "how many times intersecion(query, primitive) is called: " <<
             inter_counter << std::endl;
      output << "how many times intersecion(query, primitive) returned true: " <<
             true_inter_counter << std::endl;
      output << "how many nodes are visited in segment casting: " <<
             do_inter_in_segment << std::endl;
      output << "how many nodes are visited in ray casting: " << do_inter_in_ray <<
             std::endl;
      output << "how many nodes are visited: " << do_inter_counter << std::endl;
      //output << "how many miss occured: " << miss_counter << std::endl;

      //reset
      //miss_counter = 0;
      Traits_with_ray::inter_counter = 0;
      Traits_with_segment::inter_counter = 0;
      Traits_with_ray::true_inter_counter = 0;
      Traits_with_segment::true_inter_counter = 0;
      Traits_with_ray::do_inter_counter = 0;
      Traits_with_segment::do_inter_counter = 0;
    }
    //for maximum
    for(int i = 0; i < 5; ++i) {
      use_minimum_segment = false;
      multiplier_for_segment = 1.0 + i * 0.2;

      Timer t;
      t.start();
      calculate_sdf_values();

      double time = t.time();
      long inter_counter = Traits_with_ray::inter_counter +
                           Traits_with_segment::inter_counter;
      long true_inter_counter = Traits_with_ray::true_inter_counter +
                                Traits_with_segment::true_inter_counter;
      long do_inter_counter  = Traits_with_ray::do_inter_counter +
                               Traits_with_segment::do_inter_counter;
      long do_inter_in_segment = Traits_with_segment::do_inter_counter;
      long do_inter_in_ray = Traits_with_ray::do_inter_counter;

      output << "time: " << time << std::endl;
      output << "how many times intersecion(query, primitive) is called: " <<
             inter_counter << std::endl;
      output << "how many times intersecion(query, primitive) returned true: " <<
             true_inter_counter << std::endl;
      output << "how many nodes are visited in segment casting: " <<
             do_inter_in_segment << std::endl;
      output << "how many nodes are visited in ray casting: " << do_inter_in_ray <<
             std::endl;
      output << "how many nodes are visited: " << do_inter_counter << std::endl;
      //output << "how many miss occured: " << miss_counter << std::endl;

      //reset
      //miss_counter = 0;
      Traits_with_ray::inter_counter = 0;
      Traits_with_segment::inter_counter = 0;
      Traits_with_ray::true_inter_counter = 0;
      Traits_with_segment::true_inter_counter = 0;
      Traits_with_ray::do_inter_counter = 0;
      Traits_with_segment::do_inter_counter = 0;
    }
    output.close();
#endif
  }
};
}//namespace internal
} //namespace CGAL

#undef CGAL_NORMALIZATION_ALPHA
#undef CGAL_CONVEX_FACTOR
#undef CGAL_DEFAULT_NUMBER_OF_CLUSTERS
#undef CGAL_DEFAULT_SMOOTHING_LAMBDA
#undef CGAL_DEFAULT_CONE_ANGLE
#undef CGAL_DEFAULT_NUMBER_OF_RAYS

#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H
