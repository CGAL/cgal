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
    , m_concavity_calc(mesh, traits)
    {}

    /**
     * Computes approximate convex decomposition of a triangle mesh and fills up face property map with cluster-ids.
     * @param face_ids which associates each face of a triangle mesh to a cluster-id [0, 'number_of_clusters'-1]
     * @param concavity_threshold concavity value each cluster must satisfy
     * @param min_number_of_clusters minimal number of cluster that can be produced, except the case when number of faces is less
     */
    template <class FacePropertyMap>
    std::size_t decompose(FacePropertyMap face_ids, double concavity_threshold, std::size_t min_number_of_clusters)
    {
        // create filtered dual graph without border edges (null source or target vertex)
        Dual_graph dual(m_mesh);
        Filtered_dual_graph filtered_dual(dual, Noborder_predicate<TriangleMesh>(m_mesh));

        // fill up adjacency list and its properties
        setup_graph(filtered_dual);

        // compute initial decimation costs and concavity values for all edges
        BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
        {
            update_edge(edge, concavity_threshold, CONCAVITY_FACTOR);
        }

        // main loop that stops when all edges with concavities lower that `concavity_threshold` (valid edges)
        // are decimated or the constraint of minimum number of cluster is met
        while (num_vertices(m_graph) > min_number_of_clusters)
        {
            bool found_edge = false;

            graph_edge_descriptor optimal_edge;

            BOOST_FOREACH(graph_edge_descriptor edge, edges(m_graph))
            {
                Decimation_properties& decimation_props = m_decimation_map[edge];

                // concavity of new cluster too large (not valid)
                if (decimation_props.new_cluster_props.concavity > concavity_threshold) continue; 

                // found first valid edge or the decimation cost of new valid edge is better than the decimation cost of current optimal valid edge
                if (!found_edge || decimation_props.decimation_cost < m_decimation_map[optimal_edge].decimation_cost)
                {
                    optimal_edge = edge;
                    found_edge = true;
                }
            }

            // if no valid edges then stop
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

            // decimate optimal valid edge and update the decimation costs of modified edges
            decimate_edge(optimal_edge, concavity_threshold, CONCAVITY_FACTOR);
        }

        // resulting number of produced clusters
        std::size_t num_clusters = num_vertices(m_graph);    

        // current cluster's id
        int cluster_id = 0;

        BOOST_FOREACH(graph_vertex_descriptor vert, vertices(m_graph))
        {
            Cluster_properties& cluster_props = m_cluster_map[vert];

            // assign id to all current cluster's faces
            BOOST_FOREACH(face_descriptor face, cluster_props.faces)
            {
                face_ids[face] = cluster_id;
            }

            ++cluster_id;
        }

        // post check if all faces are assigned to any cluster
        BOOST_FOREACH(face_descriptor face, faces(m_mesh))
        {
            CGAL_assertion(face_ids[face] >= 0 && face_ids[face] < num_clusters);
        }

        return num_clusters;
    }

private:
    const TriangleMesh& m_mesh;
    const GeomTraits& m_traits;
    
    Graph m_graph; // adjacency list

    Graph_cluster_map m_cluster_map; // vertex property map of the adjacency list
    Graph_decimation_map m_decimation_map; // edge property map of the adjacency list

//    std::priority_queue<graph_edge_descriptor> m_candidates;

    Concavity<TriangleMesh, GeomTraits> m_concavity_calc; // concavity calculator that computes concavity values for any subset of faces of the input mesh

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
        std::vector<face_descriptor> faces; // list of faces of the input mesh
        std::vector<Point_3> conv_hull_pts; // list of points on the convex hull of the cluster
        CGAL::Bbox_3 bbox; // bounding box of the cluster

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
        Cluster_properties new_cluster_props; // properties of the cluster that is produced after decimation of the edge
    };

    /**
     * Fills up adjacency list.
     */
    void setup_graph(const Filtered_dual_graph& dual)
    {
        std::map<face_descriptor, graph_vertex_descriptor> face_graph_map; // maps faces of the input mesh to vetices of the adjacency list 
        
        // extract property maps of the adjacency list
        m_cluster_map = boost::get(cluster_props_t(), m_graph);
        m_decimation_map = boost::get(decimation_props_t(), m_graph);

        // add vertices
        // fill up cluster properties (single face) for all vertices
        BOOST_FOREACH(face_descriptor face, vertices(dual))
        {
            Cluster_properties props;
            
            props.concavity = 0;
            props.faces.push_back(face);

            BOOST_FOREACH(vertex_descriptor vert, vertices_around_face(halfedge(face, m_mesh), m_mesh))
            {
                props.conv_hull_pts.push_back(get(CGAL::vertex_point, m_mesh)[vert]);
            }

            props.bbox = CGAL::Polygon_mesh_processing::face_bbox(face, m_mesh);

            face_graph_map[face] = add_vertex(props, m_graph);
        }

        // add edges
        // decimation properties are filled up in the later calls to update_edge
        BOOST_FOREACH(edge_descriptor edge, edges(dual))
        {
            Decimation_properties props;

            add_edge(face_graph_map[source(edge, dual)], face_graph_map[target(edge, dual)], props, m_graph);
        }
    }

    /**
     * Constructs cluster of two smaller clusters and computes decimation cost.
     */
    void update_edge(graph_edge_descriptor edge, double concavity_threshold, double alpha_factor)
    {
        Decimation_properties& decimation_props = m_decimation_map[edge];

        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        Cluster_properties& cluster_1_props = m_cluster_map[vert_1];
        Cluster_properties& cluster_2_props = m_cluster_map[vert_2];

        // convex hull points
        std::vector<Point_3> common_hull_pts; // merged list of the lists of the convex hull points of two clusters
        
        CGAL_assertion(cluster_1_props.conv_hull_pts.size() >= 3);
        CGAL_assertion(cluster_2_props.conv_hull_pts.size() >= 3);
        
        // add the convex hull points of the first cluster
        BOOST_FOREACH(Point_3 p, cluster_1_props.conv_hull_pts)
        {
            common_hull_pts.push_back(p);
        }
        // add the convex hull points of the second cluster
        BOOST_FOREACH(Point_3 p, cluster_2_props.conv_hull_pts)
        {
            common_hull_pts.push_back(p);
        }

        // compute convex hull on the merged list of points
        TriangleMesh conv_hull;
        CGAL::convex_hull_3(common_hull_pts.begin(), common_hull_pts.end(), conv_hull);

        // fill up the list of the convex hull points of produced cluster after decimation
        decimation_props.new_cluster_props.conv_hull_pts.clear();
        BOOST_FOREACH(vertex_descriptor vert, vertices(conv_hull))
        {
            decimation_props.new_cluster_props.conv_hull_pts.push_back(get(CGAL::vertex_point, conv_hull)[vert]);
        }
        CGAL_warning(decimation_props.new_cluster_props.conv_hull_pts.size() > 3);

        // faces
        decimation_props.new_cluster_props.faces.clear();

        // add faces from the first cluster
        BOOST_FOREACH(face_descriptor face, cluster_1_props.faces)
        {
            decimation_props.new_cluster_props.faces.push_back(face);
        }
        // add faces from the second cluster
        BOOST_FOREACH(face_descriptor face, cluster_2_props.faces)
        {
            decimation_props.new_cluster_props.faces.push_back(face);
        }

        // compute concavity value of the produced cluster using concavity calculator of the input mesh, faces, and convex hull of the produced cluster
        decimation_props.new_cluster_props.concavity = m_concavity_calc.compute(decimation_props.new_cluster_props.faces, conv_hull);

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Concavity: " << decimation_props.new_cluster_props.concavity << std::endl;
#endif

        // compute bounding box of two clusters
        decimation_props.new_cluster_props.bbox = cluster_1_props.bbox + cluster_2_props.bbox;

        // compute decimation cost
        double aspect_ratio = compute_aspect_ratio(decimation_props.new_cluster_props.faces);
        double d = compute_normalization_factor(decimation_props.new_cluster_props.bbox);

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Aspect ratio & normalization factor: " << aspect_ratio << " " << d << std::endl;
#endif

        double alpha = alpha_factor * concavity_threshold / d;

        decimation_props.decimation_cost = decimation_props.new_cluster_props.concavity / d + alpha * aspect_ratio;

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Decimation cost: " << decimation_props.decimation_cost << std::endl;
#endif
    }

    /**
     * Decimates an edge in adjacency list. The second vertex of an edge merges into the first edge.
     */
    void decimate_edge(graph_edge_descriptor edge, double concavity_threshold, double alpha_factor)
    {
        graph_vertex_descriptor vert_1 = source(edge, m_graph), vert_2 = target(edge, m_graph);

        CGAL_assertion(m_decimation_map[edge].new_cluster_props.conv_hull_pts.size() > 3);
        
        // assign cluster properties of the produced cluster to the first vertex of the edge
        m_cluster_map[vert_1] = m_decimation_map[edge].new_cluster_props;

        // connect the first vertex to all adjacent vertices of the second one except self
        BOOST_FOREACH(graph_vertex_descriptor vert, boost::adjacent_vertices(vert_2, m_graph))
        {
            if (vert == vert_1) continue;
            std::pair<graph_edge_descriptor, bool> result = add_edge(vert_1, vert, m_graph);
        }

        // remove adjacent edges incident to the second vertex and remove the vertex
        clear_vertex(vert_2, m_graph);
        remove_vertex(vert_2, m_graph);

        // update decimation costs of all modified edges (only adjacent edges to the first vertex)
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

    /**
     * Compute aspect ratio of a cluster.
     */
    double compute_aspect_ratio(const std::vector<face_descriptor>& faces)
    {
        CGAL_assertion(!faces.empty());

        double perimeter = 0;
        double area = 0;

        std::set<edge_descriptor> used_edges;

        BOOST_FOREACH(face_descriptor face, faces)
        {
            // sum up area
            area += CGAL::Polygon_mesh_processing::face_area(face, m_mesh);

            // sum up lengths of edges except internal edges (they are processed twice)
            BOOST_FOREACH(edge_descriptor edge, edges_around_face(halfedge(face, m_mesh), m_mesh))
            {
                int found = used_edges.find(edge) == used_edges.end() ? 1 : -1;
                perimeter += found * CGAL::Polygon_mesh_processing::edge_length(edge, m_mesh);
                used_edges.insert(edge);    
            }
        }

#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
//        std::cout << "Perimeter & area: " << perimeter << " " << area << std::endl;
#endif

        return perimeter * perimeter / (4. * boost::math::constants::pi<double>() * area);
    }

    /**
     * Computes length of the diagonal of a bounding box (normalization factor).
     */
    double compute_normalization_factor(const CGAL::Bbox_3& bbox)
    {
        Point_3 min_p(bbox.xmin(), bbox.ymin(), bbox.zmin());
        Point_3 max_p(bbox.xmax(), bbox.ymax(), bbox.zmax());

        return CGAL::sqrt(CGAL::squared_distance(min_p, max_p));
    }
};    

}
}

#endif // CGAL_APPROX_DECOMPOSITION_H
