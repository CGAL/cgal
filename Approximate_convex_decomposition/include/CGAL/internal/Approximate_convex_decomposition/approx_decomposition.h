#ifndef CGAL_APPROX_DECOMPOSITION_H
#define CGAL_APPROX_DECOMPOSITION_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>

namespace CGAL
{
namespace internal
{

template < class TriangleMesh,
           class GeomTraits
           >
class Approx_decomposition
{
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;
 
    typedef CGAL::Surface_mesh<Point_3> SM;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    template <class Graph>
    struct Noborder_predicate;

    struct Cluster_properties;

public:
    Approx_decomposition(const TriangleMesh& mesh, const GeomTraits& traits)
    : m_mesh(mesh)
    , m_traits(traits)
    {}

    template <class FacePropertyMap, class PointPropertyMap>
    std::size_t decompose(FacePropertyMap face_ids, PointPropertyMap point_ids, double concavity_threshold, std::size_t min_number_of_clusters)
    {
        typedef CGAL::Dual<TriangleMesh> DualGraph;
        typedef boost::filtered_graph<DualGraph, Noborder_predicate<TriangleMesh>> FilteredDualGraph;

        BOOST_FOREACH(face_descriptor face, CGAL::faces(m_mesh))
        {
            face_ids[face] = -1;
        }

        DualGraph dual(m_mesh);
        FilteredDualGraph filtered_dual(dual, Noborder_predicate<TriangleMesh>(m_mesh));

//        BOOST_FOREACH(face_descriptor face, CGAL::vertices(filtered_dual))
//        {
//            std::cout << face << std::endl;
//        }

//        BOOST_FOREACH(edge_descriptor edge, CGAL::edges(filtered_dual))
//        {
//            std::cout << edge << ": " << CGAL::source(edge, m_mesh) << " -> " << CGAL::target(edge, m_mesh) << " " << CGAL::source(edge, filtered_dual) << " -> " << CGAL::target(edge, filtered_dual) << std::endl;
//        }

        typedef std::map<face_descriptor, Cluster_properties> Clusters_map;
        typedef std::pair<face_descriptor, Cluster_properties> Clusters_map_pair;
        
        Clusters_map clusters_map;

        BOOST_FOREACH(face_descriptor face, CGAL::faces(m_mesh))
        {
            clusters_map[face].concavity = 0;

            clusters_map[face].faces.push_back(face);

            BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices_around_face(CGAL::halfedge(face, m_mesh), m_mesh))
            {
//                std::cout << vert << " ";
                clusters_map[face].conv_hull_pts.push_back(m_mesh.point(vert));
            }
//            std::cout << std::endl;
        }

        bool found_edge = true;
        while (found_edge &&
               clusters_map.size() > min_number_of_clusters)
        {
            found_edge = false;

            //TODO: implement the main body
        }

        int cluster_id = 0;
        BOOST_FOREACH(Clusters_map_pair pair, clusters_map)
        {
            Cluster_properties& cluster_props = pair.second;

            if (cluster_props.faces.size() <= 1) continue;

            BOOST_FOREACH(face_descriptor face, cluster_props.faces)
            {
                face_ids[face] = cluster_id;
                
//                point_ids[] = ?; // TODO: find out which cluster_id to assign
            }

            ++cluster_id;
        }

        return cluster_id;
    }

private:
    const TriangleMesh& m_mesh;
    const GeomTraits& m_traits;

    template <class Graph>
    struct Noborder_predicate
    {
        Noborder_predicate() : m_graph(NULL) {}
        Noborder_predicate(const Graph& graph) : m_graph(graph) {}
        
        bool operator()(const edge_descriptor& edge) const { return !CGAL::is_border(edge, m_graph); }
        
        const Graph& m_graph;
    };

    struct Cluster_properties
    {
        double concavity;
        
        std::vector<face_descriptor> faces;
        
        std::vector<Point_3> conv_hull_pts;
    };
};        

}
}

#endif // CGAL_APPROX_DECOMPOSITION_H
