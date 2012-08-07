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
//AF: macros must be prefixed with "CGAL_"
//IOY: done
#define CGAL_NORMALIZATION_ALPHA 5.0
#define CGAL_CONVEX_FACTOR 0.08
#define CGAL_DEFAULT_NUMBER_OF_CLUSTERS 5
#define CGAL_DEFAULT_SMOOTHING_LAMBDA 23.0
#define CGAL_DEFAULT_CONE_ANGLE (2.0 / 3.0) * CGAL_PI
#define CGAL_DEFAULT_NUMBER_OF_RAYS 25
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

/**
 * It is a connector class which uses soft clustering and graph cut in order to segment meshes.
 * All preprocessing and postprocessing issues are handled here.
 */
template <
class Polyhedron,
      class FacetIndexMap =
      boost::associative_property_map<std::map<typename Polyhedron::Facet_const_handle, int> >
      >
class Surface_mesh_segmentation
{
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

  typedef internal::SDF_calculation<Polyhedron> SDF_calculation;

protected:
  typedef typename Kernel::Plane_3   Plane;

// member variables
public: // going to be protected !
  const Polyhedron& mesh;

  FacetIndexMap facet_index_map;
  std::map<Facet_const_handle, int> facet_index_map_internal;

  std::vector<double> sdf_values;
  std::vector<int>    centers;
  std::vector<int>    segments;
// member functions
  bool is_pmmap_custom;
  std::ofstream log_file;
public:
  Surface_mesh_segmentation(const Polyhedron& mesh)
    : mesh(mesh), facet_index_map(facet_index_map_internal),
      is_pmmap_custom(false) {
    int facet_index = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end();
        ++facet_it, ++facet_index) {
      boost::put(facet_index_map, facet_it, facet_index);
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
    mesh(rhs.mesh), sdf_values(rhs.sdf_values), centers(rhs.centers),
    segments(rhs.segments),
    is_pmmap_custom(rhs.is_pmmap_custom),
    facet_index_map_internal(rhs.facet_index_map_internal) {
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
        boost::put(facet_index_map, facet_it, facet_index);
      }
    }
  }
// Use these two functions together
  void calculate_sdf_values(double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                            int number_of_ray = CGAL_DEFAULT_NUMBER_OF_RAYS) {
    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())
    //read_sdf_values("D:/dino.txt");
    //return;
    typedef std::map<Facet_const_handle, double> internal_map;
    internal_map facet_value_map_internal;
    boost::associative_property_map<internal_map> sdf_pmap(
      facet_value_map_internal);

    SDF_calculation().calculate_sdf_values(cone_angle, number_of_ray, mesh,
                                           sdf_pmap);

    sdf_values = std::vector<double>(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      get(sdf_values, facet_it) = boost::get(sdf_pmap, facet_it);
    }

    SEG_DEBUG(std::cerr << "SDF computation time: " << t.time() << std::endl)
    SEG_DEBUG(t.reset())

    check_zero_sdf_values();
    smooth_sdf_values_with_bilateral();
    normalize_sdf_values();

    SEG_DEBUG(std::cerr <<"Normalization and smoothing time: " << t.time() <<
              std::endl)
  }

  void partition(int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
                 double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA) {
    log_file.open("log_file.txt");
    // soft clustering using GMM-fitting initialized with k-means
    internal::Expectation_maximization fitter(number_of_centers, sdf_values,
        internal::Expectation_maximization::K_MEANS_INITIALIZATION, 1);

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
    internal::Alpha_expansion_graph_cut gc(edges, edge_weights, probability_matrix,
                                           labels);
    centers = labels;
    // assign a segment id for each facet
    assign_segments();
    log_file.close();
  }
///////////////////////////////////////////////////////////////////
// Use these two functions together
  template <class SDFPropertyMap>
  void calculate_sdf_values(SDFPropertyMap sdf_pmap,
                            double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                            int number_of_rays = CGAL_DEFAULT_NUMBER_OF_RAYS) {
    SEG_DEBUG(Timer t)
    SEG_DEBUG(t.start())

    SDF_calculation().calculate_sdf_values(cone_angle, number_of_rays, mesh,
                                           sdf_pmap);

    sdf_values = std::vector<double>(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      get(sdf_values, facet_it) = boost::get(sdf_pmap, facet_it);
    }

    SEG_DEBUG(std::cout << t.time() << std::endl)

    check_zero_sdf_values();
    smooth_sdf_values_with_bilateral();
    linear_normalize_sdf_values();

    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      boost::put(sdf_pmap, facet_it, get(sdf_values, facet_it));
    }

    SEG_DEBUG(std::cout << t.time() << std::endl)
  }

  template <class FacetSegmentMap, class SDFPropertyMap>
  void partition(SDFPropertyMap sdf_pmap, FacetSegmentMap segment_pmap,
                 int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
                 double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA) {
    sdf_values = std::vector<double>(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      get(sdf_values, facet_it) = boost::get(sdf_pmap, facet_it);
    }

    // log normalize sdf values
    normalize_sdf_values();
    // soft clustering using GMM-fitting initialized with k-means
    internal::Expectation_maximization fitter(number_of_centers, sdf_values,
        internal::Expectation_maximization::K_MEANS_INITIALIZATION, 1);

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
    internal::Alpha_expansion_graph_cut gc(edges, edge_weights, probability_matrix,
                                           labels);
    centers = labels;
    // assign a segment id for each facet
    assign_segments();
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      boost::put(segment_pmap, facet_it, get(segments, facet_it));
    }
  }
///////////////////////////////////////////////////////////////////
  double get_sdf_value_of_facet(Facet_const_handle facet) const {
    return sdf_values[boost::get(facet_index_map, facet)];
  }

  int get_segment_id_of_facet(Facet_const_handle facet) const {
    return segments[boost::get(facet_index_map, facet)];
  }

protected:
  double calculate_dihedral_angle_of_edge(Edge_const_handle& edge) const {
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

  double calculate_dihedral_angle_of_edge_2(Edge_const_handle& edge) {
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
    log_file << "v2: " << angle << std::endl;
    log_file << "v1: " << calculate_dihedral_angle_of_edge(edge) << std::endl;
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

  void smooth_sdf_values() {
    std::vector<double> smoothed_sdf_values(mesh.size_of_facets());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
        facet_it->facet_begin();
      double total_neighbor_sdf = 0.0;
      do {
        total_neighbor_sdf += get(sdf_values, facet_circulator->opposite()->facet());
      } while( ++facet_circulator !=  facet_it->facet_begin());

      total_neighbor_sdf /= 3.0;
      get(smoothed_sdf_values, facet_it) = (get(sdf_values,
                                            facet_it) + total_neighbor_sdf) / 2.0;
    }
    sdf_values = smoothed_sdf_values;
  }

  void smooth_sdf_values_with_gaussian() {
    // take neighbors, use weighted average of neighbors as filtered result. (for weights use gaussian kernel with sigma = window_size/2)
    const int window_size = 2;
    const int iteration = 1;

    for(int i = 0; i < iteration; ++i) {
      std::vector<double> smoothed_sdf_values(mesh.size_of_facets());
      for(Facet_const_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end(); ++facet_it) {
        std::map<Facet_const_handle, int> neighbors;
        get_neighbors_by_vertex(facet_it, neighbors, window_size);

        double total_sdf_value = 0.0;
        double total_weight = 0.0;
        for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          double weight = gaussian_function(it->second,
                                            window_size/2.0); // window_size => 2*sigma
          total_sdf_value += get(sdf_values, it->first) * weight;
          total_weight += weight;
        }
        get(smoothed_sdf_values, facet_it) = total_sdf_value / total_weight;
      }
      sdf_values = smoothed_sdf_values;
    }
  }

  void smooth_sdf_values_with_median() {
    // take neighbors, use median sdf_value as filtered one.
    const int window_size = 2;
    const int iteration = 1;

    for(int i = 0; i < iteration; ++i) {
      std::vector<double> smoothed_sdf_values(mesh.size_of_facets());;
      for(Facet_const_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end(); ++facet_it) {
        //Find neighbors and put their sdf values into a list
        std::map<Facet_const_handle, int> neighbors;
        get_neighbors_by_vertex(facet_it, neighbors, window_size);
        std::vector<double> sdf_of_neighbors;
        for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          sdf_of_neighbors.push_back(get(sdf_values, it->first));
        }
        // Find median.
        double median_sdf = 0.0;
        int half_neighbor_count = sdf_of_neighbors.size() / 2;
        std::nth_element(sdf_of_neighbors.begin(),
                         sdf_of_neighbors.begin() + half_neighbor_count, sdf_of_neighbors.end());
        double median_sdf = sdf_of_neighbors[half_neighbor_count];
        if( half_neighbor_count % 2 == 0) {
          median_sdf += *std::max_element(sdf_of_neighbors.begin(),
                                          sdf_of_neighbors.begin() + half_neighbor_count);
          median_sdf /= 2;
        }
        get(smoothed_sdf_values, facet_it) = median_sdf;
      }
      sdf_values = smoothed_sdf_values;
    }
  }

  void smooth_sdf_values_with_bilateral() {
    // take neighbors, use weighted average of neighbors as filtered result.
    // two weights are multiplied:
    // spatial: over geodesic distances
    // domain : over sdf_value distances
    const int window_size = 2;
    const int iteration = 1;

    for(int i = 0; i < iteration; ++i) {
      std::vector<double> smoothed_sdf_values(mesh.size_of_facets());
      for(Facet_const_iterator facet_it = mesh.facets_begin();
          facet_it != mesh.facets_end(); ++facet_it) {
        //Facet_handle facet = facet_it;
        std::map<Facet_const_handle, int> neighbors;
        get_neighbors_by_edge(facet_it, neighbors, window_size);

        double total_sdf_value = 0.0, total_weight = 0.0;
        double current_sdf_value = get(sdf_values, facet_it);
        // calculate deviation for range weighting.
        double deviation = 0.0;
        for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          deviation += std::pow(get(sdf_values, it->first) - current_sdf_value, 2);
        }
        deviation = std::sqrt(deviation / neighbors.size());
        if(deviation == 0.0) {
          deviation = std::numeric_limits<double>::epsilon();  //this might happen
        }
        for(typename std::map<Facet_const_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          double spatial_weight = gaussian_function(it->second,
                                  window_size/2.0); // window_size => 2*sigma
          // deviation was std::sqrt(2.0)*deviation
          double domain_weight = gaussian_function(get(sdf_values,
                                 it->first) -  current_sdf_value, deviation);

          //double domain_weight = exp(-0.5 * (std::pow( (get(sdf_values, it->first) -  current_sdf_value) / 0.1, 2)));
          double weight = spatial_weight * domain_weight;
          total_sdf_value += get(sdf_values, it->first) * weight;
          total_weight += weight;
        }
        get(smoothed_sdf_values, facet_it) = total_sdf_value / total_weight;
      }
      sdf_values = smoothed_sdf_values;
    }
  }

  void get_neighbors_by_edge(Facet_const_handle& facet,
                             std::map<Facet_const_handle, int>& neighbors, int max_level) {
    typedef std::pair<Facet_const_handle, int> facet_level_pair;
    std::queue<facet_level_pair> facet_queue;
    facet_queue.push(facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
          pair.first->facet_begin();
        do {
          facet_queue.push(facet_level_pair(facet_circulator->opposite()->facet(),
                                            pair.second + 1));
        } while(++facet_circulator != pair.first->facet_begin());
      }
      facet_queue.pop();
    }
  }

  void get_neighbors_by_vertex(Facet_const_handle& facet,
                               std::map<Facet_const_handle, int>& neighbors, int max_level) {
    typedef std::pair<Facet_const_handle, int> facet_level_pair;
    std::queue<facet_level_pair> facet_queue;
    facet_queue.push(facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        Facet_const_handle facet = pair.first;
        Halfedge_const_iterator edge = facet->halfedge();
        do { // loop on three vertices of the facet
          Vertex_const_iterator vertex = edge->vertex();
          typename Facet::Halfedge_around_vertex_const_circulator vertex_circulator =
            vertex->vertex_begin();
          do { // for each vertex loop on incoming edges (through those edges loop on neighbor facets which includes the vertex)
            facet_queue.push(facet_level_pair(vertex_circulator->facet(), pair.second + 1));
          } while(++vertex_circulator != vertex->vertex_begin());
        } while((edge = edge->next()) != facet->halfedge());
      }
      facet_queue.pop();
    }
  }

  void check_zero_sdf_values() {
    // If there is any facet which has no sdf value, assign average sdf value of its neighbors
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      if(get(sdf_values, facet_it) == 0.0) {
        typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
          facet_it->facet_begin();
        double total_neighbor_sdf = 0.0;
        do {
          total_neighbor_sdf += get(sdf_values, facet_circulator->opposite()->facet());
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

  void log_normalize_probability_matrix(std::vector<std::vector<double> >&
                                        probabilities) const {
    const double epsilon = 1e-5;
    for(std::vector<std::vector<double> >::iterator it_i = probabilities.begin();
        it_i != probabilities.end(); ++it_i) {
      for(std::vector<double>::iterator it = it_i->begin(); it != it_i->end(); ++it) {
        double probability = (std::max)(*it, epsilon);
        probability = -log(probability);
        *it = (std::max)(probability, std::numeric_limits<double>::epsilon());
      }
    }
  }

  void calculate_and_log_normalize_dihedral_angles(double smoothing_lambda,
      std::vector<std::pair<int, int> >& edges, std::vector<double>& edge_weights) {
    const double epsilon = 1e-5;
    //edges and their weights. pair<int, int> stores facet-id pairs (see above) (may be using boost::tuple can be more suitable)
    for(Edge_const_iterator edge_it = mesh.edges_begin();
        edge_it != mesh.edges_end(); ++edge_it) {
      const int index_f1 = boost::get(facet_index_map, edge_it->facet());
      const int index_f2 = boost::get(facet_index_map, edge_it->opposite()->facet());
      edges.push_back(std::pair<int, int>(index_f1, index_f2));

      double angle = calculate_dihedral_angle_of_edge_2(edge_it);
      log_file << "angle received: " << angle << std::endl;
      if(angle < epsilon) {
        angle = epsilon;
      }
      log_file << "angle after < epsilon-check: " << angle << std::endl;
      angle = -log(angle);
      log_file << "angle after minus-log: " << angle << std::endl;
      angle = (std::max)(angle, std::numeric_limits<double>::epsilon());
      log_file << "angle after > epsilon-check: " << angle << std::endl;
      angle *= smoothing_lambda;
      log_file << "angle after multiplied by lambda: " << angle << std::endl;
      edge_weights.push_back(angle);
    }
  }

  void assign_segments() {
    segments = std::vector<int>(mesh.size_of_facets(), -1);
    int segment_id = 0;
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      if(get(segments, facet_it) == -1) {
        depth_first_traversal(facet_it, segment_id);
        segment_id++;
      }
    }
  }

  void depth_first_traversal(Facet_const_handle& facet, int segment_id) {
    get(segments, facet) = segment_id;
    typename Facet::Halfedge_around_facet_const_circulator facet_circulator =
      facet->facet_begin();
    double total_neighbor_sdf = 0.0;
    do {
      Facet_const_handle neighbor = facet_circulator->opposite()->facet();
      if(get(segments, neighbor) == -1
          && get(centers, facet) == get(centers, neighbor)) {
        depth_first_traversal(neighbor, segment_id);
      }
    } while( ++facet_circulator !=  facet->facet_begin());
  }



  template<class T>
  T& get(std::vector<T>& data, const Facet_const_handle& facet) {
    return data[ boost::get(facet_index_map, facet) ];
  }

  double gaussian_function(double value, double deviation) {
    return exp(-0.5 * (std::pow(value / deviation, 2)));
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
