#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char **argv)
{
    if (argc !=2){
        std::cerr << "Usage: " << argv[0] << " input.OFF" << std::endl;
        return 1;
    }
  
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input(argv[1]);
    
    if ( !input || !(input >> mesh) || mesh.empty() ){
        std::cerr << argv[1] << " is not a valid off file." << std::endl;
        return 1;
    }

    // create a property-map for sdf values (it is an adaptor for this case) 
    typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_sdf_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

    // compute sdf values using default parameters for number of rays, and cone angle
    CGAL::compute_sdf_values(mesh, sdf_property_map);

    // create a property-map for segment-ids (it is an adaptor for this case)
    typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

    // segment the mesh using default parameters for number of levels, and smoothing lambda
    // Note that you can use your own scalar value, instead of using SDF calculation computed using the CGAL function
    int number_of_segments = CGAL::segment_from_sdf_values(mesh, sdf_property_map, segment_property_map);

    std::cout << "Number of segments: " << number_of_segments << std::endl;
    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
        // ids are between [0, number_of_segments -1]
        std::cout << segment_property_map[facet_it] << std::endl;                                 
    }

    const int number_of_levels = 4;       // use 4 clusters in soft clustering
    const double smoothing_lambda = 0.3;  // importance of surface features, between [0,1]

    // Note that we can use same sdf values (sdf_property_map) over and over again for segmentation.
    // This feature becomes important when we want to segment the mesh several times with different parameters.
    CGAL::segment_from_sdf_values(
      mesh, sdf_property_map, segment_property_map, number_of_levels, smoothing_lambda);
}