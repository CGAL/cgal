#include <CGAL/internal/Surface_mesh_segmentation/Filters.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <cmath>
#include <map>
#include <CGAL/property_map.h>

#include "Utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

/**
 * Uses bilateral and median filtering to smooth associated values with facets.
 * It also prints average differences between them (smoothed values with two methods).
 * Note that it always return EXIT_SUCCESS if .off file is read successfully.
 */
int main(void)
{	
	Polyhedron mesh;
	if( !read_to_polyhedron("./data/cactus.off", mesh) ) { return 1; }
		
	typedef std::map< Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_1;
    boost::associative_property_map<Facet_double_map> value_pmap_1(internal_1);
	
	Facet_double_map internal_2;
    boost::associative_property_map<Facet_double_map> value_pmap_2(internal_2);
	
	double counter = 0.0;
	for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it, ++counter)   
    {
        value_pmap_1[facet_it] = value_pmap_2[facet_it] = counter;                               
    }
	const double average_value = counter / 2.0;
	
	CGAL::internal::Bilateral_filtering<Polyhedron, CGAL::internal::Neighbor_selector_by_edge<Polyhedron> > filter_1;
	CGAL::internal::Bilateral_filtering<Polyhedron, CGAL::internal::Neighbor_selector_by_vertex<Polyhedron> > filter_2;
	CGAL::internal::Median_filtering<Polyhedron, CGAL::internal::Neighbor_selector_by_edge<Polyhedron> > filter_3;
	CGAL::internal::Median_filtering<Polyhedron, CGAL::internal::Neighbor_selector_by_vertex<Polyhedron> > filter_4;
	
	filter_1(mesh, 2, value_pmap_1);
	filter_3(mesh, 2, value_pmap_2);
	
	double average_dif = 0.0;
	for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
		average_dif += std::abs(value_pmap_1[facet_it] - value_pmap_2[facet_it]);
	}		
	average_dif = (average_dif / mesh.size_of_facets()) / average_value;
	std::cout << "average differences between bilateral and median filters: " << average_dif << std::endl;
	
	filter_2(mesh, 2, value_pmap_1);
	filter_4(mesh, 2, value_pmap_2);
	
	average_dif = 0.0;
	for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
		average_dif += std::abs(value_pmap_1[facet_it] - value_pmap_2[facet_it]);
	}		
	average_dif = (average_dif / mesh.size_of_facets()) / average_value;
	std::cout << "average differences between bilateral and median filters (2) :" << average_dif << std::endl;	
}
