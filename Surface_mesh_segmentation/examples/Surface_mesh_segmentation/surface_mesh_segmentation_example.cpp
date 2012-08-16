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

    // create a property-map for segment-ids (it is an adaptor for this case)
    typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
	
    // calculate SDF values and segment the mesh using default parameters.
    CGAL::surface_mesh_segmentation(mesh, segment_property_map);

    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); 
        facet_it != mesh.facets_end(); ++facet_it)   
    {
        std::cout << segment_property_map[facet_it] << std::endl;                                 
    }
}