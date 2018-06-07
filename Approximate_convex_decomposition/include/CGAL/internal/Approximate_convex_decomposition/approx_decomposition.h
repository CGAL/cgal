#ifndef CGAL_APPROX_DECOMPOSITION_H
#define CGAL_APPROX_DECOMPOSITION_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Bbox_3.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/copy.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>

#include "concavity.h"

namespace CGAL
{
namespace internal
{

template < class TriangleMesh,
           class GeomTraits
           >
class Approx_decomposition
{
    // property tags
    struct cluster_props_t
    {
        typedef boost::vertex_property_tag kind;
    };

    struct decimation_props_t
    {
        typedef boost::edge_property_tag kind;
    };

    // predefined structs
    template <class Graph>
    struct Noborder_predicate;

    struct Cluster_properties;
    struct Decimation_properties;
    
    // typedefs
    typedef typename GeomTraits::Point_3 Point_3;
    typedef typename GeomTraits::Vector_3 Vector_3;
 
    typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

    typedef CGAL::Face_filtered_graph<Surface_mesh> Face_filtered_mesh;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    typedef CGAL::Dual<TriangleMesh> Dual_graph;
    typedef boost::filtered_graph<Dual_graph, Noborder_predicate<TriangleMesh>> Filtered_dual_graph;
    
    typedef boost::property<cluster_props_t, Cluster_properties> VertexProperty;
    typedef boost::property<decimation_props_t, Decimation_properties> EdgeProperty;

    typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, VertexProperty, EdgeProperty> Graph;
    
    typedef typename boost::graph_traits<Graph>::vertex_descriptor graph_vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::edge_descriptor graph_edge_descriptor;

    typedef typename boost::property_map<Graph, cluster_props_t>::type Graph_cluster_map;
    typedef typename boost::property_map<Graph, decimation_props_t>::type Graph_decimation_map;

    // constants
    const double CONCAVITY_FACTOR = 0.1;

public:
    Approx_decomposition(const TriangleMesh& mesh, const GeomTraits& traits)
    : m_mesh(mesh)
    , m_traits(traits)
    {}

    ~Approx_decomposition()
    {}

    template <class FacePropertyMap>
    std::size_t decompose(FacePropertyMap face_ids, double concavity_threshold, std::size_t min_number_of_clusters)
    {
        BOOST_FOREACH(face_descriptor face, CGAL::faces(m_mesh))
        {
            face_ids[face] = -1;
        }

        Dual_graph dual(m_mesh);
        Filtered_dual_graph filtered_dual(dual, Noborder_predicate<TriangleMesh>(m_mesh));
        
        setup_graph(filtered_dual);

//        BOOST_FOREACH(edge_descriptor edge, edges(orig_dual))
//        {
//            std::cout << edge << ": " << source(edge, m_mesh) << " -> " << target(edge, m_mesh) << " " << source(edge, orig_dual) << " -> " << target(edge, orig_dual) << std::endl;
//            CGAL_assertion(!CGAL::is_border(edge, m_mesh));
//        }

//        BOOST_FOREACH(face_descriptor face, vertices(m_graph))
//        {
//            std::cout << face << std::endl;
//        }

//        BOOST_FOREACH(edge_descriptor edge, edges(m_graph))
//        {
//            std::cout << edge << ": " << source(edge, m_mesh) << " -> " << target(edge, m_mesh) << " " << source(edge, m_graph) << " -> " << target(edge, m_graph) << std::endl;
//            CGAL_assertion(!CGAL::is_border(edge, m_mesh));
//        }

        BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
        {
            update_edge(edge, concavity_threshold, CONCAVITY_FACTOR);
        }

        while (num_vertices(m_graph) > min_number_of_clusters)
        {
            bool found_edge = false;

            graph_edge_descriptor optimal_edge;


            BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
            {
                Decimation_properties& decimation_props = m_decimation_map[edge];

                if (decimation_props.new_cluster_props.concavity > concavity_threshold) continue;

                if (!found_edge || decimation_props.decimation_cost < m_decimation_map[optimal_edge].decimation_cost)
                {
                    optimal_edge = edge;
                    found_edge = true;
                }
            }

            if (!found_edge) break;

            std::cout << std::endl << "Optimal edge for decimation: " << optimal_edge << " " << m_decimation_map[optimal_edge].decimation_cost << " " << m_decimation_map[optimal_edge].new_cluster_props.concavity << std::endl;
            std::cout << "Vertices: " << num_vertices(m_graph) << " Edges: " << num_edges(m_graph) << std::endl;
//            BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
//            {
//                std::cout << source(edge, m_graph) << " " << target(edge, m_graph) << std::endl;
////                std::cout << edge << std::endl;
//            }

            decimate_edge(optimal_edge, concavity_threshold, CONCAVITY_FACTOR);
        }

        int cluster_id = 0;
        BOOST_FOREACH(graph_vertex_descriptor vert, vertices(m_graph))
        {
            Cluster_properties& cluster_props = m_cluster_map[vert];

//            if (cluster_props.faces.size() <= 1) continue;

            BOOST_FOREACH(face_descriptor face, cluster_props.faces)
            {
                face_ids[face] = cluster_id;
            }

            ++cluster_id;
        }

        return cluster_id;
    }

private:
    const TriangleMesh& m_mesh;
    const GeomTraits& m_traits;
    
    Graph m_graph;

    Graph_cluster_map m_cluster_map;
    Graph_decimation_map m_decimation_map;

    std::priority_queue<graph_edge_descriptor> m_candidates;

    template <class Mesh>
    struct Noborder_predicate
    {
        Noborder_predicate() : m_mesh(NULL) {}
        Noborder_predicate(const Mesh& mesh) : m_mesh(mesh) {}
        
        bool operator()(const edge_descriptor& edge) const { return !CGAL::is_border(edge, m_mesh); }
        
        const Mesh& m_mesh;
    };

    struct Cluster_properties
    {
        double concavity;
        std::vector<face_descriptor> faces;
        std::vector<Point_3> conv_hull_pts;

//        void operator=(const Cluster_properties&& r)
//        {
//            concavity = r.concavity;
//            faces = std::move(r.faces);
//            conv_hull_pts = std::move(r.conv_hull_pts);
//        }
    };

    struct Decimation_properties
    {
        double decimation_cost;
        Cluster_properties new_cluster_props;
    };

    void setup_graph(const Filtered_dual_graph& dual)
    {
        std::map<face_descriptor, graph_vertex_descriptor> face_graph_map;

        m_cluster_map = boost::get(cluster_props_t(), m_graph);
        m_decimation_map = boost::get(decimation_props_t(), m_graph);

        BOOST_FOREACH(face_descriptor face, vertices(dual))
        {
            Cluster_properties props;
            
            props.concavity = 0;
            props.faces.push_back(face);

            BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices_around_face(CGAL::halfedge(face, m_mesh), m_mesh))
            {
                props.conv_hull_pts.push_back(m_mesh.point(vert));
            }
            CGAL_assertion(props.conv_hull_pts.size() >= 3);

            face_graph_map[face] = add_vertex(props, m_graph);
        }

        BOOST_FOREACH(edge_descriptor edge, edges(dual))
        {
            Decimation_properties props;

            add_edge(face_graph_map[source(edge, dual)], face_graph_map[target(edge, dual)], props, m_graph);
        }
    }

    void update_edge(graph_edge_descriptor edge, double concavity_threshold, double alpha_factor)
    {
//        std::cout << std::endl << "Update edge: " << edge << std::endl;

        Decimation_properties& decimation_props = m_decimation_map[edge];

        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        Cluster_properties& cluster_1_props = m_cluster_map[vert_1];
        Cluster_properties& cluster_2_props = m_cluster_map[vert_2];

        std::vector<Point_3> common_hull_pts;
//        std::cout << cluster_1_props.conv_hull_pts.size() << " " << cluster_2_props.conv_hull_pts.size() << std::endl;
        CGAL_assertion(cluster_1_props.conv_hull_pts.size() >= 3);
        CGAL_assertion(cluster_2_props.conv_hull_pts.size() >= 3);
        BOOST_FOREACH(Point_3 p, cluster_1_props.conv_hull_pts)
        {
            common_hull_pts.push_back(p);
        }
        BOOST_FOREACH(Point_3 p, cluster_2_props.conv_hull_pts)
        {
            common_hull_pts.push_back(p);
        }

        CGAL_assertion(common_hull_pts.size() > 3);

        Surface_mesh conv_hull;
        CGAL::convex_hull_3(common_hull_pts.begin(), common_hull_pts.end(), conv_hull);

        decimation_props.new_cluster_props.conv_hull_pts.clear();
        BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices(conv_hull))
        {
            decimation_props.new_cluster_props.conv_hull_pts.push_back(conv_hull.point(vert));
        }
        CGAL_warning(decimation_props.new_cluster_props.conv_hull_pts.size() > 3);

        decimation_props.new_cluster_props.faces.clear();
        BOOST_FOREACH(face_descriptor face, cluster_1_props.faces)
        {
            decimation_props.new_cluster_props.faces.push_back(face);
        }
        BOOST_FOREACH(face_descriptor face, cluster_2_props.faces)
        {
            decimation_props.new_cluster_props.faces.push_back(face);
        }

        Surface_mesh new_cluster;
        Face_filtered_mesh selected_faces(m_mesh, decimation_props.new_cluster_props.faces);
        CGAL::copy_face_graph(selected_faces, new_cluster);
//        {
//            std::ofstream os("new_cluster_" + std::to_string(edge) + ".off");
//            os << new_cluster;
//        }

        Concavity<Surface_mesh, GeomTraits> concavity_calc(new_cluster, m_traits);
        decimation_props.new_cluster_props.concavity = concavity_calc.compute(conv_hull);

//        std::cout << "Concavity: " << decimation_props.new_cluster_props.concavity << std::endl;

        double aspect_ratio = compute_aspect_ratio(new_cluster);
        double d = compute_normalization_factor(new_cluster);

//        std::cout << "Aspect ratio & normalization factor: " << aspect_ratio << " " << d << std::endl;

        double alpha = alpha_factor * concavity_threshold / d;

        decimation_props.decimation_cost = decimation_props.new_cluster_props.concavity / d + alpha * aspect_ratio;

//        std::cout << "Decimation cost: " << decimation_props.decimation_cost << std::endl;
    }

    void decimate_edge(graph_edge_descriptor edge, double concavity_threshold, double alpha_factor)
    {
//        std::cout << "Decimate edge: " << edge << std::endl;

        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        CGAL_assertion(vert_1 != vert_2);

        CGAL_assertion(m_decimation_map[edge].new_cluster_props.conv_hull_pts.size() > 3);
        m_cluster_map[vert_1] = m_decimation_map[edge].new_cluster_props;
        CGAL_assertion(m_cluster_map[vert_1].conv_hull_pts.size() > 3);

        BOOST_FOREACH(graph_vertex_descriptor vert, boost::adjacent_vertices(vert_2, m_graph))
        {
            if (vert == vert_1) continue;
            std::pair<graph_edge_descriptor, bool> result = add_edge(vert_1, vert, m_graph);
//            if (!result.second) std::cout << "!!!!" << std::endl;
        }

        clear_vertex(vert_2, m_graph);
        remove_vertex(vert_2, m_graph);

        CGAL_assertion(m_cluster_map[vert_1].conv_hull_pts.size() > 3);
       
        int cnt = 0;
        BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_1, m_graph))
        {
            update_edge(edge, concavity_threshold, alpha_factor);
            ++cnt;
        }
        std::cout << "Updated edges: " << cnt << std::endl;
    }

    double compute_aspect_ratio(const Surface_mesh& mesh)
    {
        double perimeter = 0;
        double area = 0;

        BOOST_FOREACH(edge_descriptor edge, edges(mesh))
        {
            if (!CGAL::is_border(edge, mesh)) continue;
            
//            const Point_3& a = mesh.point(source(edge, mesh));
//            const Point_3& b = mesh.point(target(edge, mesh));

//            perimeter += CGAL::sqrt(CGAL::squared_distance(a, b));
            perimeter += CGAL::Polygon_mesh_processing::edge_length(edge, mesh);
        }

        BOOST_FOREACH(face_descriptor face, faces(mesh))
        {
            area += CGAL::Polygon_mesh_processing::face_area(face, mesh); 
        }

//        std::cout << "Perimeter & area: " << perimeter << " " << area << std::endl;

        return perimeter * perimeter / (4. * boost::math::constants::pi<double>() * area);
    }

    double compute_normalization_factor(const Surface_mesh& mesh)
    {
        CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);

        Point_3 min_p(bbox.xmin(), bbox.ymin(), bbox.zmin());
        Point_3 max_p(bbox.xmax(), bbox.ymax(), bbox.zmax());

        return CGAL::sqrt(CGAL::squared_distance(min_p, max_p));
    }
};    

}
}

#endif // CGAL_APPROX_DECOMPOSITION_H
