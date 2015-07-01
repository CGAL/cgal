#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Timer.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

template<class ValueType, class Polyhedron>
struct Facet_with_id_pmap
  : public boost::put_get_helper<ValueType&,
  Facet_with_id_pmap<ValueType, Polyhedron> >
{
  typedef typename Polyhedron::Facet_const_handle key_type;
  typedef ValueType value_type;
  typedef value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  Facet_with_id_pmap(
    std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }

  reference operator[](key_type key) const
  { return internal_vector[key->id()]; }
private:
  std::vector<ValueType>& internal_vector;
};

template<class Polyhedron>
double compute_and_time_sdf(const Polyhedron& mesh, Facet_with_id_pmap<double, Polyhedron> sdf_values) 
{
  CGAL::Timer timer;
  timer.start();
  std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_values);
  timer.stop();
  std::cout << "  minimum SDF: " << min_max_sdf.first << " maximum SDF: " << min_max_sdf.second << std::endl;
  std::cout << "  SDF Time: " << timer.time() << std::endl;
  return timer.time();
}

template<class Polyhedron>
double compute_and_time_segmentation(const Polyhedron& mesh, Facet_with_id_pmap<double, Polyhedron> sdf_values, std::size_t cluster) {
  // Segment
  std::vector<std::size_t> internal_segment_map(mesh.size_of_facets());
  Facet_with_id_pmap<std::size_t, Polyhedron> segment_property_map(internal_segment_map);

  CGAL::Timer timer;
  timer.start();
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_values, segment_property_map, cluster);
  timer.stop();
  std::cout << "  number of segments: " << number_of_segments << std::endl;
  std::cout << "  Segmentation Time with " << cluster << " clusters : " << timer.time() << std::endl;
  return timer.time();
}

template<class Polyhedron>
std::vector<double>
read_and_run(const std::string& file_name) {
  std::cout << " File name " << file_name;
  std::vector<double> timings;
  
  Polyhedron mesh;
  std::ifstream input(file_name.c_str());
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    assert(false);
  }
  std::cout << " Number of Facets: " << mesh.size_of_facets() << std::endl;
  timings.push_back( mesh.size_of_facets());
  // assign id field for each facet
  std::size_t facet_id = 0;
  for(typename Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
    facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
      facet_it->id() = facet_id;
  }

  std::vector<double> sdf_values_internal(mesh.size_of_facets());
  Facet_with_id_pmap<double, Polyhedron> sdf_values(sdf_values_internal);
  double sdf_time = compute_and_time_sdf(mesh, sdf_values);
  timings.push_back(sdf_time);

  std::vector<std::size_t> cluster_counts;
  cluster_counts.push_back(2);
  cluster_counts.push_back(5);
  cluster_counts.push_back(10);
  cluster_counts.push_back(15);

  for(std::size_t i = 0; i < cluster_counts.size(); ++i) {
    double segment_time = compute_and_time_segmentation(mesh, sdf_values, cluster_counts[i]);
    timings.push_back(segment_time);
  }
  return timings;
}

int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;
  typedef CGAL::Simple_cartesian<double>                      Simple_K;
  typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>    EPICK_Polyhedron;
  typedef CGAL::Polyhedron_3<Simple_K, CGAL::Polyhedron_items_with_id_3> Simple_K_Polyhedron;

  std::vector<std::string> files;
  files.push_back("data/dino.off");
  files.push_back("data/bear.off");
  files.push_back("data/refined_elephant.off");

  std::vector<std::vector<double> > timings_simple;
  std::vector<std::vector<double> > timings_epick;

  for(std::size_t i = 0; i < files.size(); ++i) {
    std::cout << "Simple_cartesian ";
    timings_simple.push_back( read_and_run<Simple_K_Polyhedron>(files[i]) );
    std::cout << "EPICK ";
    timings_epick.push_back( read_and_run<EPICK_Polyhedron>(files[i]) );
  }
  // print timings
  std::cout << std::endl << std::endl;
  std::cout << "*** SDF Timings ****" << std::endl;
  std::cout << "+----------+------------------+-------+" << std::endl;
  std::cout << "| # Facets | Simple Cartesian | EPICK |" << std::endl;
  std::cout << "+----------+------------------+-------+" << std::endl;
  for(std::size_t i = 0; i < timings_simple.size(); ++i) {
    std::cout << "| "
              << std::setw(9) << timings_simple[i][0] << "|"
              << std::setw(18) << timings_simple[i][1] << "|" 
              << std::setw(7)  << timings_epick[i][1]  << "|" << std::endl; 
  }
  std::cout << "+----------+------------------+-------+" << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "*** Segmentation Timings (Simple Cartesian) ****" << std::endl;
  std::cout << "+----------+------------+------------+-------------+-------------+" << std::endl;
  std::cout << "| # Facets | #Cluster=2 | #Cluster=5 | #Cluster=10 | #Cluster=15 |" << std::endl;
  std::cout << "+----------+------------+------------+-------------+-------------+" << std::endl;
  for(std::size_t i = 0; i < timings_simple.size(); ++i) {
    std::cout << "| " 
              << std::setw(9) << timings_simple[i][0] << "|"
              << std::setw(12) << timings_simple[i][2] << "|"
              << std::setw(12) << timings_simple[i][3] << "|"
              << std::setw(13) << timings_simple[i][4] << "|"
              << std::setw(13) << timings_simple[i][5] << "|" << std::endl;
  }
  std::cout << "+----------+------------+------------+-------------+-------------+" << std::endl;
}
