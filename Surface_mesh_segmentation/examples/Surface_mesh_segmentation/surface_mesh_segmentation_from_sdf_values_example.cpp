#include <iostream>
#include <fstream>
#include <cstdlib>

#include <CGAL/mesh_segmentation.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/property_map/property_map.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char **argv)
{
    // create and read Polyhedron
    Polyhedron mesh; 	
    std::ifstream(argv[1]) >> mesh;

    // create a property-map for sdf values (it is an adaptor for this case) 
    typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_sdf_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);
	
    // compute sdf values using default parameters for number of rays, and cone angle
    CGAL::sdf_values_computation(mesh, sdf_property_map);
	
    // create a property-map for segment-ids (it is an adaptor for this case)
    typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
	
    // segment the mesh using default parameters for number of levels, and smoothing lambda	
    CGAL::surface_mesh_segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
        std::cout << segment_property_map[facet_it] << std::endl;                                 
    }

    int number_of_levels = 4;       // use 4 clusters in soft clustering
    double smoothing_lambda = 20.0; // importance of surface features 

    // Note that we can use same sdf values (sdf_property_map) over and over again for segmentation.
    // This feature becomes important when we want to segment the mesh several times with different parameters.
    CGAL::surface_mesh_segmentation_from_sdf_values(
	mesh, sdf_property_map, segment_property_map, number_of_levels, smoothing_lambda);
}