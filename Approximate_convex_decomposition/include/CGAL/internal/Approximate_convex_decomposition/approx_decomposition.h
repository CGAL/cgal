#ifndef CGAL_APPROX_DECOMPOSITION_H
#define CGAL_APPROX_DECOMPOSITION_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/properties.h>
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

#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>

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
 
    typedef CGAL::Face_filtered_graph<TriangleMesh> Face_filtered_mesh;

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

    typedef typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type Vertex_point_map;

    // constants
    const double CONCAVITY_FACTOR = 0.1;

public:
    Approx_decomposition(const TriangleMesh& mesh, const GeomTraits& traits)
    : m_mesh(mesh)
    , m_traits(traits)
    , m_concavity_calc(mesh, traits)
    {}

    template <class FacePropertyMap>
    std::size_t decompose(FacePropertyMap face_ids, double concavity_threshold, std::size_t min_number_of_clusters)
    {
        Dual_graph dual(m_mesh);
        Filtered_dual_graph filtered_dual(dual, Noborder_predicate<TriangleMesh>(m_mesh));

//        m_points_map = get(CGAL::vertex_point, m_mesh);
        auto m = get(CGAL::vertex_point, m_mesh);

        setup_graph(filtered_dual);

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

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE            
//            std::cout << std::endl;
            std::cout << "#" << num_vertices(m_graph) << " Optimal edge for decimation: " << optimal_edge << std::endl;
            std::cout << "Decimation cost: " << m_decimation_map[optimal_edge].decimation_cost << ", Concavity value: " << m_decimation_map[optimal_edge].new_cluster_props.concavity << std::endl;
            std::cout << "Total edges: " << num_edges(m_graph) << std::endl;
//            BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
//            {
//                std::cout << source(edge, m_graph) << " " << target(edge, m_graph) << std::endl;
//            }
#endif

            decimate_edge(optimal_edge, concavity_threshold, CONCAVITY_FACTOR);
        }

        std::size_t num_clusters = num_vertices(m_graph);
        
        int cluster_id = 0;
        BOOST_FOREACH(graph_vertex_descriptor vert, vertices(m_graph))
        {
            Cluster_properties& cluster_props = m_cluster_map[vert];

            BOOST_FOREACH(face_descriptor face, cluster_props.faces)
            {
                face_ids[face] = cluster_id;
            }

            ++cluster_id;
        }

        BOOST_FOREACH(face_descriptor face, faces(m_mesh))
        {
            CGAL_assertion(face_ids[face] >= 0 && face_ids[face] < num_clusters);
        }

        return num_clusters;
    }

private:
    const TriangleMesh& m_mesh;
    const GeomTraits& m_traits;
    
    Graph m_graph;

    Graph_cluster_map m_cluster_map;
    Graph_decimation_map m_decimation_map;

    Vertex_point_map m_points_map;

//    std::priority_queue<graph_edge_descriptor> m_candidates;

    Concavity<TriangleMesh, GeomTraits> m_concavity_calc;

    template <class Mesh>
    struct Noborder_predicate
    {
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

            BOOST_FOREACH(vertex_descriptor vert, vertices_around_face(halfedge(face, m_mesh), m_mesh))
            {
//                props.conv_hull_pts.push_back(m_points_map[vert]);
                props.conv_hull_pts.push_back(get(CGAL::vertex_point, m_mesh)[vert]);
            }

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
        Decimation_properties& decimation_props = m_decimation_map[edge];

        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        Cluster_properties& cluster_1_props = m_cluster_map[vert_1];
        Cluster_properties& cluster_2_props = m_cluster_map[vert_2];

        std::vector<Point_3> common_hull_pts;
        
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

        TriangleMesh conv_hull;
        CGAL::convex_hull_3(common_hull_pts.begin(), common_hull_pts.end(), conv_hull);

        decimation_props.new_cluster_props.conv_hull_pts.clear();
        BOOST_FOREACH(vertex_descriptor vert, vertices(conv_hull))
        {
            decimation_props.new_cluster_props.conv_hull_pts.push_back(get(CGAL::vertex_point, conv_hull)[vert]);
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

        TriangleMesh new_cluster;
        Face_filtered_mesh selected_faces(m_mesh, decimation_props.new_cluster_props.faces);
        CGAL::copy_face_graph(selected_faces, new_cluster);

//        Concavity<TriangleMesh, GeomTraits> concavity_calc(new_cluster, m_traits);
        decimation_props.new_cluster_props.concavity = m_concavity_calc.compute(vertices(new_cluster), conv_hull);

//        decimation_props.new_cluster_props.concavity = m_concavity_calc.compute(decimation_props.new_cluster_props.faces, conv_hull);

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Concavity: " << decimation_props.new_cluster_props.concavity << std::endl;
#endif

        double aspect_ratio = compute_aspect_ratio(new_cluster);
        double d = compute_normalization_factor(new_cluster);

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Aspect ratio & normalization factor: " << aspect_ratio << " " << d << std::endl;
#endif

        double alpha = alpha_factor * concavity_threshold / d;

        decimation_props.decimation_cost = decimation_props.new_cluster_props.concavity / d + alpha * aspect_ratio;

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Decimation cost: " << decimation_props.decimation_cost << std::endl;
#endif
    }

    void decimate_edge(graph_edge_descriptor edge, double concavity_threshold, double alpha_factor)
    {
        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        CGAL_assertion(m_decimation_map[edge].new_cluster_props.conv_hull_pts.size() > 3);
        m_cluster_map[vert_1] = m_decimation_map[edge].new_cluster_props;

        BOOST_FOREACH(graph_vertex_descriptor vert, boost::adjacent_vertices(vert_2, m_graph))
        {
            if (vert == vert_1) continue;
            std::pair<graph_edge_descriptor, bool> result = add_edge(vert_1, vert, m_graph);
        }

        clear_vertex(vert_2, m_graph);
        remove_vertex(vert_2, m_graph);

        int cnt = 0;
        BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_1, m_graph))
        {
            update_edge(edge, concavity_threshold, alpha_factor);
            ++cnt;
        }
#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
        std::cout << "Updated edges: " << cnt << std::endl;
#endif
    }

    double compute_aspect_ratio(const TriangleMesh& mesh)
    {
        double perimeter = 0;
        double area = 0;

        BOOST_FOREACH(edge_descriptor edge, edges(mesh))
        {
            if (!CGAL::is_border(edge, mesh)) continue;
            
            perimeter += CGAL::Polygon_mesh_processing::edge_length(edge, mesh);
        }

        BOOST_FOREACH(face_descriptor face, faces(mesh))
        {
            area += CGAL::Polygon_mesh_processing::face_area(face, mesh); 
        }

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Perimeter & area: " << perimeter << " " << area << std::endl;
#endif

        return perimeter * perimeter / (4. * boost::math::constants::pi<double>() * area);
    }

    double compute_normalization_factor(const TriangleMesh& mesh)
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
