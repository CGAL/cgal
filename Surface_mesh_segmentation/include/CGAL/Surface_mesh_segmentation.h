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

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>

//AF: macros must be prefixed with "CGAL_"
//IOY: done
#define CGAL_NORMALIZATION_ALPHA 5.0
#define CGAL_CONVEX_FACTOR 0.08

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
template <class Polyhedron, class SDFCalculation = internal::SDF_calculation<Polyhedron> >
class Surface_mesh_segmentation
{
//type definitions
public:
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Polyhedron::Facet  Vertex;
  typedef typename Kernel::Vector_3   Vector;
  typedef typename Kernel::Point_3    Point;
  typedef typename Polyhedron::Facet_iterator  Facet_iterator;
  typedef typename Polyhedron::Facet_handle    Facet_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Edge_iterator   Edge_iterator;
  typedef typename Polyhedron::Vertex_handle   Vertex_handle;

  typedef typename SDFCalculation::Parameters  SDF_Parameters;

protected:
  typedef typename Kernel::Plane_3   Plane;

  typedef std::map<Facet_handle, double>    Face_value_map;
  typedef std::map<Facet_handle, int>       Face_center_map;
  typedef std::map<Facet_handle, int>       Face_segment_map;

  template <typename ValueTypeName>
  struct Compare_second_element {
    bool operator()(const ValueTypeName& v1, const ValueTypeName& v2) const {
      return v1.second < v2.second;
    }
  };

//member variables
public:
  Polyhedron* mesh;

  Face_value_map   sdf_values;
  Face_center_map  centers;
  Face_segment_map segments;

  int    number_of_centers;
  double smoothing_lambda;


  internal::Expectation_maximization fitter;/**< going to be removed */

//member functions
public:
  Surface_mesh_segmentation(Polyhedron* mesh): mesh(mesh) {
#ifdef SEGMENTATION_PROFILE
    profile("profile.txt");
#endif
  }

  void calculate_sdf_values(SDF_Parameters parameters = SDF_Parameters()) {
    SEG_DEBUG(CGAL::Timer t)
    SEG_DEBUG(t.start())

    sdf_values.clear();
    SDFCalculation(parameters).calculate_sdf_values(*mesh, sdf_values);

    SEG_DEBUG(std::cout << t.time() << std::endl)

    check_zero_sdf_values();
    smooth_sdf_values_with_bilateral();
    normalize_sdf_values();

    SEG_DEBUG(std::cout << t.time() << std::endl)
  }

  void partition(int number_of_centers = 5, double smoothing_lambda = 23.0) {
    this->number_of_centers = number_of_centers;
    this->smoothing_lambda = smoothing_lambda;
    centers.clear();

    std::vector<double> sdf_vector;
    sdf_vector.reserve(sdf_values.size());
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it) {
      sdf_vector.push_back(sdf_values[facet_it]);
    }
    // soft clustering using GMM-fitting initialized with k-means
    fitter = internal::Expectation_maximization(number_of_centers, sdf_vector,
             internal::Expectation_maximization::K_MEANS_INITIALIZATION, 1);

    std::vector<int> labels;
    fitter.fill_with_center_ids(labels);

    std::vector<std::vector<double> > probability_matrix;
    fitter.fill_with_probabilities(probability_matrix);
    log_normalize_probability_matrix(probability_matrix);

    // calculating edge weights
    std::vector<std::pair<int, int> > edges;
    std::vector<double> edge_weights;
    calculate_and_log_normalize_dihedral_angles(edges, edge_weights);

    // apply graph cut
    internal::Alpha_expansion_graph_cut gc(edges, edge_weights, probability_matrix,
                                           labels);

    std::vector<int>::iterator center_it = labels.begin();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
    // assign a segment id for each facet
    assign_segments();
  }

public:
  double calculate_dihedral_angle_of_edge(const Halfedge_handle& edge) const {
    Facet_handle f1 = edge->facet();
    Facet_handle f2 = edge->opposite()->facet();

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
    Vector f1_normal = CGAL::unit_normal(f1_v1, f1_v2, f1_v3);
    Vector f2_normal = CGAL::unit_normal(f2_v1, f2_v2, f2_v3);

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

  double calculate_dihedral_angle_of_edge_2(const Halfedge_handle& edge) const {
    const Point& a = edge->vertex()->point();
    const Point& b = edge->prev()->vertex()->point();
    const Point& c = edge->next()->vertex()->point();
    const Point& d = edge->opposite()->next()->vertex()->point();
    // As far as I check: if, say, dihedral angle is 5, this returns 175,
    // if dihedral angle is -5, this returns -175.
    double n_angle = CGAL::Mesh_3::dihedral_angle(a, b, c, d) / 180.0;
    bool concave = n_angle > 0;
    double angle = 1 + ((concave ? -1 : +1) * n_angle);

    if(!concave) {
      angle *= CGAL_CONVEX_FACTOR;
    }
    return angle;

    //Facet_handle f1 = edge->facet();
    //Facet_handle f2 = edge->opposite()->facet();
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
    //Vector f1_normal = CGAL::unit_normal(f1_v1, f1_v2, f1_v3);
    //Vector f2_normal = CGAL::unit_normal(f2_v1, f2_v2, f2_v3);
    //
    //double dot = f1_normal * f2_normal;
    //if(dot > 1.0)       { dot = 1.0;  }
    //else if(dot < -1.0) { dot = -1.0; }
    //double angle = acos(dot) / CGAL_PI; // [0-1] normalize
    //std::cout << angle << " " << n_angle << " " << (concave ? "concave": "convex") << std::endl;
    //if(fabs(angle - folded_angle) > 1e-6)
    //{
    //
    //}
    //if(angle < epsilon) { angle = epsilon; }
    //return angle;
  }

  void normalize_sdf_values() {
    typedef typename Face_value_map::iterator fv_iterator;
    Compare_second_element<typename Face_value_map::value_type> comparator;
    std::pair<fv_iterator, fv_iterator> min_max_pair =
      CGAL::min_max_element(sdf_values.begin(), sdf_values.end(), comparator,
                            comparator);

    double max_value = min_max_pair.second->second,
           min_value = min_max_pair.first->second;
    double max_min_dif = max_value - min_value;
    for(fv_iterator pair_it = sdf_values.begin(); pair_it != sdf_values.end();
        ++pair_it) {
      double linear_normalized = (pair_it->second - min_value) / max_min_dif;
      double log_normalized = log(linear_normalized * CGAL_NORMALIZATION_ALPHA + 1) /
                              log(CGAL_NORMALIZATION_ALPHA + 1);
      pair_it->second = log_normalized;
    }
  }

  void smooth_sdf_values() {
    Face_value_map smoothed_sdf_values;
    for(typename Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      Facet_handle f = pair_it->first;
      typename Facet::Halfedge_around_facet_circulator facet_circulator =
        f->facet_begin();
      double total_neighbor_sdf = 0.0;
      do {
        total_neighbor_sdf += sdf_values[facet_circulator->opposite()->facet()];
      } while( ++facet_circulator !=  f->facet_begin());
      total_neighbor_sdf /= 3.0;
      smoothed_sdf_values[f] = (pair_it->second + total_neighbor_sdf) / 2.0;
    }
    sdf_values = smoothed_sdf_values;
  }

  void smooth_sdf_values_with_gaussian() {
    // take neighbors, use weighted average of neighbors as filtered result. (for weights use gaussian kernel with sigma = window_size/2)
    const int window_size = 2;
    const int iteration = 1;

    for(int i = 0; i < iteration; ++i) {
      Face_value_map smoothed_sdf_values;
      for(typename Face_value_map::iterator pair_it = sdf_values.begin();
          pair_it != sdf_values.end(); ++pair_it) {
        Facet_handle facet = pair_it->first;
        std::map<Facet_handle, int> neighbors;
        get_neighbors_by_vertex(facet, neighbors, window_size);

        double total_sdf_value = 0.0;
        double total_weight = 0.0;
        for(std::map<Facet_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          double weight =  exp(-0.5 * (std::pow(it->second / (window_size/2.0),
                                                2))); // window_size => 2*sigma
          total_sdf_value += sdf_values[it->first] * weight;
          total_weight += weight;
        }
        smoothed_sdf_values[facet] = total_sdf_value / total_weight;
      }
      sdf_values = smoothed_sdf_values;
    }
  }

  void smooth_sdf_values_with_median() {
    // take neighbors, use median sdf_value as filtered one.
    const int window_size = 2;
    const int iteration = 1;

    for(int i = 0; i < iteration; ++i) {
      Face_value_map smoothed_sdf_values;
      for(typename Face_value_map::iterator pair_it = sdf_values.begin();
          pair_it != sdf_values.end(); ++pair_it) {
        //Find neighbors and put their sdf values into a list
        Facet_handle facet = pair_it->first;
        std::map<Facet_handle, int> neighbors;
        get_neighbors_by_vertex(facet, neighbors, window_size);
        std::vector<double> sdf_of_neighbors;
        for(std::map<Facet_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          sdf_of_neighbors.push_back(sdf_values[it->first]);
        }
        // Find median.
        double median_sdf = 0.0;
        int half_neighbor_count = sdf_of_neighbors.size() / 2;
        std::nth_element(sdf_of_neighbors.begin(),
                         sdf_of_neighbors.begin() + half_neighbor_count, sdf_of_neighbors.end());
        if( half_neighbor_count % 2 == 0) {
          double median_1 = sdf_of_neighbors[half_neighbor_count];
          double median_2 = *std::max_element(sdf_of_neighbors.begin(),
                                              sdf_of_neighbors.begin() + half_neighbor_count);
          median_sdf = (median_1 + median_2) / 2;
        } else {
          median_sdf = sdf_of_neighbors[half_neighbor_count];
        }
        smoothed_sdf_values[facet] = median_sdf;
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
      Face_value_map smoothed_sdf_values;
      for(typename Face_value_map::iterator pair_it = sdf_values.begin();
          pair_it != sdf_values.end(); ++pair_it) {
        Facet_handle facet = pair_it->first;
        std::map<Facet_handle, int> neighbors;
        get_neighbors_by_vertex(facet, neighbors, window_size);

        double total_sdf_value = 0.0, total_weight = 0.0;
        double current_sdf_value = sdf_values[facet];
        // calculate deviation for range weighting.
        double deviation = 0.0;
        for(std::map<Facet_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          deviation += std::pow(sdf_values[it->first] - current_sdf_value, 2);
        }
        deviation = std::sqrt(deviation / neighbors.size());
        if(deviation == 0.0) {
          deviation = std::numeric_limits<double>::epsilon();  //this might happen
        }
        for(std::map<Facet_handle, int>::iterator it = neighbors.begin();
            it != neighbors.end(); ++it) {
          double spatial_weight =  exp(-0.5 * (std::pow(it->second / (window_size/2.0),
                                               2))); // window_size => 2*sigma
          double domain_weight  =  exp(-0.5 * (std::pow((sdf_values[it->first] -
                                               current_sdf_value) / (std::sqrt(2.0)*deviation), 2)));
          double weight = spatial_weight * domain_weight;
          total_sdf_value += sdf_values[it->first] * weight;
          total_weight += weight;
        }
        smoothed_sdf_values[facet] = total_sdf_value / total_weight;
      }
      sdf_values = smoothed_sdf_values;
    }
  }

  void get_neighbors_by_edge(const Facet_handle& facet,
                             std::map<Facet_handle, int>& neighbors, int max_level) {
    typedef std::pair<Facet_handle, int> facet_level_pair;
    std::queue<facet_level_pair> facet_queue;
    facet_queue.push(facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        typename Facet::Halfedge_around_facet_circulator facet_circulator =
          pair.first->facet_begin();
        do {
          facet_queue.push(facet_level_pair(facet_circulator->opposite()->facet(),
                                            pair.second + 1));
        } while(++facet_circulator != pair.first->facet_begin());
      }
      facet_queue.pop();
    }
  }

  void get_neighbors_by_vertex(const Facet_handle& facet,
                               std::map<Facet_handle, int>& neighbors, int max_level) {
    typedef std::pair<Facet_handle, int> facet_level_pair;
    std::queue<facet_level_pair> facet_queue;
    facet_queue.push(facet_level_pair(facet, 0));
    while(!facet_queue.empty()) {
      const facet_level_pair& pair = facet_queue.front();
      bool inserted = neighbors.insert(pair).second;
      if(inserted && pair.second < max_level) {
        Facet_handle facet = pair.first;
        Halfedge_handle edge = facet->halfedge();
        do { // loop on three vertices of the facet
          Vertex_handle vertex = edge->vertex();
          typename Facet::Halfedge_around_vertex_circulator vertex_circulator =
            vertex->vertex_begin();
          do { // for each vertex loop on neighbor vertices
            facet_queue.push(facet_level_pair(vertex_circulator->facet(), pair.second + 1));
          } while(++vertex_circulator != vertex->vertex_begin());
        } while((edge = edge->next()) != facet->halfedge());
      }
      facet_queue.pop();
    }
  }
  void check_zero_sdf_values() {
    // If there is any facet which has no sdf value, assign average sdf value of its neighbors
    for(typename Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      if(pair_it->second == 0.0) {
        typename Facet::Halfedge_around_facet_circulator facet_circulator =
          pair_it->first->facet_begin();
        double total_neighbor_sdf = 0.0;
        do {
          total_neighbor_sdf += sdf_values[facet_circulator->opposite()->facet()];
        } while( ++facet_circulator !=  pair_it->first->facet_begin());
        pair_it->second = total_neighbor_sdf / 3.0;
      }
    }
    // Find minimum sdf value other than 0
    double min_sdf = (std::numeric_limits<double>::max)();
    for(typename Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      if(pair_it->second < min_sdf && pair_it->second != 0.0) {
        min_sdf = pair_it->second;
      }
    }
    // If still there is any facet which has no sdf value, assign minimum sdf value.
    // This is meaningful since (being an outlier) 0 sdf values might effect normalization & log extremely.
    for(typename Face_value_map::iterator pair_it = sdf_values.begin();
        pair_it != sdf_values.end(); ++pair_it) {
      if(pair_it->second == 0.0) {
        pair_it->second = min_sdf;
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
        //probability += epsilon;
        //probability = (std::min)(probability, 1.0);
        probability = -log(probability);
        *it = (std::max)(probability, std::numeric_limits<double>::epsilon());
      }
    }
  }

  void calculate_and_log_normalize_dihedral_angles(
    std::vector<std::pair<int, int> >& edges,
    std::vector<double>& edge_weights) const {
    const double epsilon = 1e-5;
    //assign an id for every facet (facet-id)
    std::map<Facet_handle, int> facet_indices;
    int index = 0;
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++index) {
      facet_indices.insert(std::pair<Facet_handle, int>(facet_it, index));
    }
    //edges and their weights. pair<int, int> stores facet-id pairs (see above) (may be using CGAL::Triple can be more suitable)
    for(Edge_iterator edge_it = mesh->edges_begin(); edge_it != mesh->edges_end();
        ++edge_it) {
      int index_f1 = facet_indices[edge_it->facet()];
      int index_f2 = facet_indices[edge_it->opposite()->facet()];
      edges.push_back(std::pair<int, int>(index_f1, index_f2));

      double angle = calculate_dihedral_angle_of_edge(edge_it);
      if(angle < epsilon) {
        angle = epsilon;
      }
      angle = -log(angle);
      angle = (std::max)(angle, std::numeric_limits<double>::epsilon());
      angle *= smoothing_lambda;
      edge_weights.push_back(angle);
    }
  }

  void assign_segments() {
    segments.clear();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      segments.insert(std::pair<Facet_handle, int>(facet_it, -1));
    }
    int segment_id = 0;
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      if(segments[facet_it] == -1) {
        depth_first_traversal(facet_it, segment_id);
        segment_id++;
      }
    }
  }

  void depth_first_traversal(const Facet_handle& facet, int segment_id) {
    segments[facet] = segment_id;
    typename Facet::Halfedge_around_facet_circulator facet_circulator =
      facet->facet_begin();
    double total_neighbor_sdf = 0.0;
    do {
      Facet_handle neighbor = facet_circulator->opposite()->facet();
      if(segments[neighbor] == -1 && centers[facet] == centers[neighbor]) {
        depth_first_traversal(neighbor, segment_id);
      }
    } while( ++facet_circulator !=  facet->facet_begin());
  }

  /**
   * Going to be removed
   */
  void apply_GMM_fitting() {
    centers.clear();
    std::vector<double> sdf_vector;
    sdf_vector.reserve(sdf_values.size());
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it) {
      sdf_vector.push_back(sdf_values[facet_it]);
    }
    SEG_DEBUG(CGAL::Timer t)
    SEG_DEBUG(t.start())
    //internal::Expectation_maximization fitter(number_of_centers, sdf_vector, 10);

    fitter = internal::Expectation_maximization(number_of_centers, sdf_vector);
    //fitter = internal::Expectation_maximization(number_of_centers, sdf_vector, 40);
    SEG_DEBUG(std::cout << "GMM fitting time: " << t.time() << std::endl)
    std::vector<int> center_memberships;
    fitter.fill_with_center_ids(center_memberships);
    std::vector<int>::iterator center_it = center_memberships.begin();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
  }
  /**
   * Going to be removed
   */
  void apply_K_means_clustering() {
    centers.clear();
    std::vector<double> sdf_vector;
    sdf_vector.reserve(sdf_values.size());
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it) {
      sdf_vector.push_back(sdf_values[facet_it]);
    }
    internal::K_means_clustering clusterer(number_of_centers, sdf_vector);
    std::vector<int> center_memberships;
    clusterer.fill_with_center_ids(center_memberships);
    std::vector<int>::iterator center_it = center_memberships.begin();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
    //center_memberships_temp = center_memberships; //remove
  }
  /**
   * Going to be removed
   */
  void apply_GMM_fitting_with_K_means_init() {
    centers.clear();
    std::vector<double> sdf_vector;
    sdf_vector.reserve(sdf_values.size());
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it) {
      sdf_vector.push_back(sdf_values[facet_it]);
    }
    std::vector<int> center_memberships;
    fitter = internal::Expectation_maximization(number_of_centers, sdf_vector,
             internal::Expectation_maximization::K_MEANS_INITIALIZATION, 1);
    center_memberships.clear();
    fitter.fill_with_center_ids(center_memberships);
    std::vector<int>::iterator center_it = center_memberships.begin();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
  }
  /**
   * Going to be removed
   */
  void apply_GMM_fitting_and_K_means() {
    centers.clear();
    std::vector<double> sdf_vector;
    sdf_vector.reserve(sdf_values.size());
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it) {
      sdf_vector.push_back(sdf_values[facet_it]);
    }
    internal::Expectation_maximization gmm_random_init(number_of_centers,
        sdf_vector);

    internal::K_means_clustering k_means(number_of_centers, sdf_vector);
    std::vector<int> center_memberships;
    k_means.fill_with_center_ids(center_memberships);
    internal::Expectation_maximization gmm_k_means_init(number_of_centers,
        sdf_vector,
        internal::Expectation_maximization::K_MEANS_INITIALIZATION, 1);

    if(gmm_k_means_init.final_likelihood > gmm_random_init.final_likelihood) {
      fitter = gmm_k_means_init;
    } else {
      fitter = gmm_random_init;
    }
    center_memberships.clear();
    fitter.fill_with_center_ids(center_memberships);
    std::vector<int>::iterator center_it = center_memberships.begin();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
  }
  /**
   * Going to be removed
   */
  void apply_graph_cut() {

    std::vector<std::pair<int, int> > edges;
    std::vector<double> edge_weights;
    calculate_and_log_normalize_dihedral_angles(edges, edge_weights);


    std::vector<std::vector<double> > probability_matrix;
    fitter.fill_with_probabilities(probability_matrix);

    std::vector<int> labels;
    fitter.fill_with_center_ids(labels);

    log_normalize_probability_matrix(probability_matrix);
    //////////////////////////////////////////////////////////////
    // FOR READING FROM MATLAB, GOING TO BE REMOVED
//   std::ifstream cc("D:/GSoC/Matlab/ccount.txt");
//   cc >> number_of_centers;
//   std::vector<std::vector<double> > probability_matrix(number_of_centers, std::vector<double>(sdf_values.size(), 0.0));
//   read_center_ids("D:/GSoC/Matlab/result.txt");
//   read_probabilities("D:/GSoC/Matlab/probs.txt", probability_matrix);
//   int f = 0;
//   for (Facet_iterator facet_it = mesh->facets_begin(); facet_it != mesh->facets_end(); facet_it++, f++)
    //{
    //    for(int i = 0; i < number_of_centers; ++i)
    //    {
    //        double probability = probability_matrix[i][f];
    //        probability += 1e-8;
//           probability = (std::min)(probability, 1.0);
//           probability = -log(probability);
//           probability_matrix[i][f] = (std::max)(probability, std::numeric_limits<double>::epsilon());
    //    }
    //}
//
//   std::vector<int> labels;
//   for(Facet_iterator facet_it = mesh->facets_begin(); facet_it != mesh->facets_end();
//        ++facet_it)
//   {
//       labels.push_back(centers[facet_it]);
//   }
    //////////////////////////////////////////////////////////////

    //apply graph cut
    internal::Alpha_expansion_graph_cut gc(edges, edge_weights, probability_matrix,
                                           labels);

    std::vector<int>::iterator center_it = labels.begin();
    centers.clear();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end();
        ++facet_it, ++center_it) {
      centers.insert(std::pair<Facet_handle, int>(facet_it, (*center_it)));
    }
  }

  /**
   * Going to be removed
   */
  void select_cluster_number() {
    int min_cluster_count = 3;
    int max_cluster_count = 5;
    int range = max_cluster_count - min_cluster_count + 1;
    std::vector<double> distortions(range+1);
    for(int i = min_cluster_count -1; i <= max_cluster_count; ++i) {
      number_of_centers = i;
      apply_GMM_fitting_with_K_means_init();
      double distortion = fitter.calculate_distortion();
      distortions[i-(min_cluster_count -1)] = std::pow(distortion, -0.5);
    }
    double max_jump = 0.0;
    for(int i = 1; i < range+1; ++i) {
      double jump = distortions[i] - distortions[i-1];
      if(jump > max_jump) {
        max_jump = jump;
        number_of_centers = i + min_cluster_count - 1;
      }
    }
    for(int i = 0; i < range + 1; ++i) {
      std::cout << "d: " << distortions[i] << std::endl;
    }
    std::cout << "noc: " << number_of_centers << std::endl;
    //apply_GMM_fitting_and_K_means_init();
  }

  /**
   * Going to be removed
   */
  void write_sdf_values(const char* file_name) {
    std::ofstream output(file_name);
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      output << sdf_values[facet_it] << std::endl;
    }
    output.close();
  }
  /**
   * Going to be removed
   */
  void read_sdf_values(const char* file_name) {
    std::ifstream input(file_name);
    sdf_values.clear();
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      double sdf_value;
      input >> sdf_value;
      sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf_value));
    }
  }
  /**
   * Going to be removed
   */
  void read_center_ids(const char* file_name) {
    std::ifstream input(file_name);
    centers.clear();
    int max_center = 0;
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      int center_id;
      input >> center_id;
      centers.insert(std::pair<Facet_handle, int>(facet_it, center_id));
      if(center_id > max_center) {
        max_center = center_id;
      }
    }
    number_of_centers = max_center + 1;
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
    for(Facet_iterator facet_it = mesh->facets_begin();
        facet_it != mesh->facets_end(); ++facet_it) {
      output << segments[facet_it] << std::endl;
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

      CGAL::Timer t;
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

      CGAL::Timer t;
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

#ifdef SEG_DEBUG
#undef SEG_DEBUG
#endif
#endif //CGAL_SURFACE_MESH_SEGMENTATION_H
