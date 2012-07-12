#include <iostream>
#include <fstream>
#include <cstdlib>

#include <CGAL/Surface_mesh_segmentation.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh_segmentation<Polyhedron> Segmentation;

int main(int argc, char **argv)
{
    std::ifstream meshFile("camel.off");
    if(!meshFile) 
    { 
        std::cerr << "Could not open the file!" << std::endl;
        return EXIT_FAILURE;
    }
    
    Polyhedron mesh;   
    meshFile >> mesh;
    Segmentation segmentation(&mesh);
	
	#if 1 
	//Store your own sdf values.
	typedef Polyhedron::Facet_iterator Facet_iterator;
	typedef Polyhedron::Facet_handle   Facet_handle;
	for(Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it)
    {
        double sdf = 0.0; // YOUR SDF VALUES - GOES HERE...
        segmentation.sdf_values.insert(std::pair<Facet_handle, double>(facet_it, sdf));
    }
	// If your sdf values are normalized then comment two lines below.
	segmentation.smooth_sdf_values_with_bilateral();	
	segmentation.normalize_sdf_values();
	
	#else
	segmentation.calculate_sdf_values();
	#endif
	
	// Optional changes
	// segmentation.number_of_centers = 5;   // You can change number of levels.
	// segmentation.smoothing_lambda  = 23.0 // You can change smoothness.
	
	segmentation.apply_GMM_fitting_with_K_means_init();
	segmentation.apply_graph_cut();
	segmentation.assign_segments();
	
	// Print segment-id for each facet (segment is a disconnected component)
	for(Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it)
    {
		std::cout << segmentation.segments[facet_it] << std::endl;
	}

}