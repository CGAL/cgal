#include <CGAL/approximate_convex_decomposition.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef CGAL::Face_filtered_graph<Mesh> Filtered_graph;

typedef CGAL::Real_timer Timer;

#ifndef CGAL_LINKED_WITH_TBB
typedef CGAL::Sequential_tag Concurrency_tag;
#else
typedef CGAL::Parallel_tag Concurrency_tag;
#endif

int main()
{
    // read mesh
    Mesh mesh;
    
    std::ifstream input("data/sword.off");
    
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
    typedef Mesh::Property_map<face_descriptor, int> Clusters_pmap;
    Clusters_pmap clusters_pmap = mesh.add_property_map<face_descriptor, int>("f:cluster").first;

    // decompose mesh
    Timer timer;
    
    timer.start();
    std::size_t clusters_num = CGAL::approximate_convex_decomposition<Concurrency_tag>(mesh, clusters_pmap, 0.3, 1);
    timer.stop();

    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;

    // write cluster-ids for each facet
    std::cout << "Number of clusters: " << clusters_num << std::endl;
    BOOST_FOREACH(face_descriptor fd, faces(mesh))
    {
        std::cout << clusters_pmap[fd] << " ";
    }
    std::cout << std::endl;

    // write concavity values for all clusters
    for (std::size_t i = 0; i < clusters_num; ++i)
    {
        std::cout << "Concavity value of #" << i << " cluster: " << CGAL::concavity_values<Concurrency_tag>(mesh, clusters_pmap, i) << std::endl;
    }

    // compute convex hulls
    for (std::size_t i = 0; i < clusters_num; ++i)
    {
        Filtered_graph filtered_mesh(mesh, i, clusters_pmap);
        Mesh cluster;
        CGAL::copy_face_graph(filtered_mesh, cluster);
   
        Mesh conv_hull;
        std::vector<Point_3> pts;

        if (CGAL::num_vertices(cluster) > 3)
        {
            BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices(cluster))
            {
                pts.push_back(cluster.point(vert));
            }

            CGAL::convex_hull_3(pts.begin(), pts.end(), conv_hull);
        }
        else             
        {                
            conv_hull = cluster;             
        }            
    
        // use convex hull
        // ...    
    }

    return EXIT_SUCCESS;
}

