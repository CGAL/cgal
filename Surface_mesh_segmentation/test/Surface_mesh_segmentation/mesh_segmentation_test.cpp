#include <CGAL/mesh_segmentation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <cmath>
#include <map>
#include <boost/property_map/property_map.hpp>

#include "Utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

/**
 * Note that it always return EXIT_SUCCESS if .off file is read successfully.
 */
int main(void)
{	
    Polyhedron mesh;
    if( !read_to_polyhedron("./data/cactus.off", mesh) ) { return 1; }
  
    typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);
  
    std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);
    std::cout << "minimum sdf: " << min_max_sdf.first << " maximum sdf: " << min_max_sdf.second << std::endl;
  
    typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
  
    int nb_segments = CGAL::segmentation_from_sdf_values(
      mesh, sdf_property_map, segment_property_map);
    
    if(nb_segments != 3)
    {
        std::cerr << "Number of segments should be 3 for cactus model (since it is pretty easy model to segment)" << std::endl;
        return EXIT_FAILURE;
    }
}