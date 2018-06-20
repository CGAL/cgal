#define CGAL_APPROX_DECOMPOSITION_VERBOSE

#include <CGAL/approx_decomposition.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

int main()
{
    // read mesh
    Polyhedron mesh;
   
    std::ifstream input("data/cube.off"); 
//    std::ifstream input("data/elephant.off");
    
    if (!input || !(input >> mesh))
    {
        std::cout << "Failed to read mesh" << std::endl;
        return EXIT_FAILURE;
    }

    if (CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cout << "Input mesh is invalid" << std::endl;
        return EXIT_FAILURE;
    }

    // create property map for cluster-ids
    typedef boost::property_map<Polyhedron, CGAL::dynamic_face_property_t<int> >::type Facet_property_map;
    Facet_property_map facet_property_map = get(CGAL::dynamic_face_property_t<int>(), mesh);

    // decompose mesh with default parameters
    auto start = std::chrono::system_clock::now();
    
    std::size_t clusters_num = CGAL::convex_decomposition(mesh, facet_property_map);

    auto end = std::chrono::system_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000. << " seconds" << std::endl;

     // write cluster-ids for each facet
    std::cout << "Number of clusters: " << clusters_num << std::endl;
    BOOST_FOREACH(face_descriptor face, faces(mesh))
    {
        std::cout << get(facet_property_map, face) << " ";
    }
    std::cout << std::endl;

    // write concavity values for all clusters
//    for (std::size_t i = 0; i < clusters_num; ++i)
//    {
//        std::cout << "#" << i << ": " << CGAL::concavity_value(mesh, facet_property_map, i) << std::endl;
//    }

    return EXIT_SUCCESS;
}

