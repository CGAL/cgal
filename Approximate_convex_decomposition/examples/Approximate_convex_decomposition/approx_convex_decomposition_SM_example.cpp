#define CGAL_APPROX_DECOMPOSITION_VERBOSE

#include <CGAL/approx_decomposition.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <chrono>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

typedef CGAL::Face_filtered_graph<Mesh> Filtered_graph;

int main()
{
    // read mesh
    Mesh mesh;
    
    std::ifstream input("data/cube.off");
//    std::ifstream input("data/sword.off");
//    std::ifstream input("data/cactus.off");
//    std::ifstream input("data/cheese.off");
//    std::ifstream input("data/lion-head.off");
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
    Mesh::Property_map<face_descriptor, int> facet_property_map;
    facet_property_map = mesh.add_property_map<face_descriptor, int>("f:cluster").first;

    // decompose mesh
    auto start = std::chrono::system_clock::now();

    std::size_t clusters_num = CGAL::convex_decomposition(mesh, facet_property_map, 0.3, 1);

    auto end = std::chrono::system_clock::now();
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000. << " seconds" << std::endl;

    // write cluster-ids for each facet
    std::cout << "Number of clusters: " << clusters_num << std::endl;
    BOOST_FOREACH(face_descriptor fd, faces(mesh))
    {
        std::cout << facet_property_map[fd] << " ";
    }
    std::cout << std::endl;

    // write concavity values for all clusters
    for (std::size_t i = 0; i < clusters_num; ++i)
    {
        std::cout << "Concavity value of #" << i << " cluster: " << CGAL::concavity_value(mesh, facet_property_map, i) << std::endl;
    }

    // write clusters to .off files
    for (std::size_t i = 0; i < clusters_num; ++i)
    {
        Filtered_graph filtered_mesh(mesh, i, facet_property_map);
        Mesh cluster;
        CGAL::copy_face_graph(filtered_mesh, cluster);
   
        {
            std::ofstream os("cluster_" + std::to_string(i) + ".off");
            os << cluster;
        }
        {
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
            std::ofstream os("ch_cluster_" + std::to_string(i) + ".off");             
            os << conv_hull;
        }
    }

    return EXIT_SUCCESS;
}

