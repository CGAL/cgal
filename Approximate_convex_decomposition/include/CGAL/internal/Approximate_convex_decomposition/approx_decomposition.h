#ifndef CGAL_APPROX_DECOMPOSITION_H
#define CGAL_APPROX_DECOMPOSITION_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/move/move.hpp>
#include <boost/unordered_map.hpp>
#include <tbb/parallel_for_each.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <utility>
#include <functional>

#include <CGAL/internal/Approximate_convex_decomposition/concavity.h>

namespace CGAL
{
namespace internal
{

template < class TriangleMesh,
           class Vpm,
           class GeomTraits,
           class ConcurrencyTag
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
    struct CandidateComparator;
    
    template <class Graph>
    struct Noborder_predicate;

    struct Cluster_properties;
    struct Decimation_properties;
    
    // typedefs
    typedef typename GeomTraits::Point_3 Point_3;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    typedef CGAL::Dual<TriangleMesh> Dual_graph;
    typedef boost::filtered_graph<Dual_graph, Noborder_predicate<TriangleMesh> > Filtered_dual_graph;
    
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
    Approx_decomposition(const TriangleMesh& mesh, const Vpm& vpm, const GeomTraits& traits)
    : m_mesh(mesh)
    , m_vpm(vpm)
    , m_traits(traits)
    , m_candidates(CandidateComparator(*this))
    , m_concavity_calc(mesh, vpm, traits)
    {}

    /**
     * Computes approximate convex decomposition of a triangle mesh and fills up face property map with cluster-ids.
     * @param face_ids which associates each face of a triangle mesh to a cluster-id [0, 'number_of_clusters'-1]
     * @param concavity_threshold concavity value each cluster must satisfy
     * @param min_number_of_clusters minimal number of cluster that can be produced
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
            add_candidate(edge, concavity_threshold);
        }
        
        // remove edges that can not become candidates
        remove_invalid_edges();

        // main loop that stops when all edges with concavities lower that `concavity_threshold` (valid edges)
        // are decimated or the constraint of minimum number of cluster is met
        while (num_vertices(m_graph) > min_number_of_clusters)
        {
            CGAL_assertion(m_candidates.size() == num_edges(m_graph));
  
            // if no valid edges left than stop
            if (m_candidates.empty()) break;

            // pop edge with the lowest decimation cost (optimal edge)
            graph_edge_descriptor optimal_edge = *m_candidates.begin();
            m_candidates.erase(m_candidates.begin());
            
#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE            
            std::cout << "#" << num_vertices(m_graph) << " Optimal edge for decimation: " << optimal_edge << std::endl;
            std::cout << "Decimation cost: " << m_decimation_map[optimal_edge].decimation_cost << ", Concavity value: " << m_decimation_map[optimal_edge].new_cluster_props.concavity << std::endl;
            std::cout << "Total edges: " << num_edges(m_graph) << std::endl;
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
                put(face_ids, face, cluster_id);
            }

            ++cluster_id;
        }

        // post check if all faces are assigned to any cluster
        BOOST_FOREACH(face_descriptor face, faces(m_mesh))
        {
            CGAL_assertion(get(face_ids, face) >= 0 && get(face_ids, face) < num_clusters);
        }

        return num_clusters;
    }

private:
    const TriangleMesh& m_mesh;
    Vpm m_vpm;
    const GeomTraits& m_traits;
    
    Graph m_graph; // adjacency list

    Graph_cluster_map m_cluster_map; // vertex property map of the adjacency list
    Graph_decimation_map m_decimation_map; // edge property map of the adjacency list
       
    std::vector<graph_edge_descriptor> m_invalid_edges; // list of edges to be removed

    // comparator that orders candidate edges by their decimation costs
    struct CandidateComparator
    {
        CandidateComparator(Approx_decomposition& alg) : m_alg(alg) {}

        bool operator() (const graph_edge_descriptor& a, const graph_edge_descriptor& b) const
        {
            if (m_alg.m_decimation_map[a].decimation_cost != m_alg.m_decimation_map[b].decimation_cost)
            {
                return m_alg.m_decimation_map[a].decimation_cost < m_alg.m_decimation_map[b].decimation_cost;
            }
            return a < b;
        }

    private:
        Approx_decomposition& m_alg;
    };
    
    std::set<graph_edge_descriptor, CandidateComparator> m_candidates; // ordered by decimation cost list of edges

    Concavity<TriangleMesh, Vpm, GeomTraits, ConcurrencyTag> m_concavity_calc; // concavity calculator that computes concavity values for any subset of faces of the input mesh

    // a predicate for border edges removal
    template <class Mesh>
    struct Noborder_predicate
    {
        Noborder_predicate(const Mesh& mesh) : m_mesh(mesh) {}
        
        bool operator()(const edge_descriptor& edge) const { return !CGAL::is_border(edge, m_mesh); }
        
        const Mesh& m_mesh;
    };

    // all the necessary information to describe a cluster 
    struct Cluster_properties
    {
        int id; // needed to prevent dupblications of edges
        double concavity;
        std::vector<face_descriptor> faces; // list of faces of the input mesh
        std::vector<Point_3> conv_hull_pts; // list of points on the convex hull of the cluster
        CGAL::Bbox_3 bbox; // bounding box of the cluster
    };

    // all the necessary information to describe a decimation operation of two clusters
    // while computing the decimation_cost, information about the merged cluster is computed and is stored in the new_cluster_props
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
        boost::unordered_map<face_descriptor, graph_vertex_descriptor> face_graph_map; // maps faces of the input mesh to vetices of the adjacency list 
        
        // extract property maps of the adjacency list
        m_cluster_map = boost::get(cluster_props_t(), m_graph);
        m_decimation_map = boost::get(decimation_props_t(), m_graph);

        // add vertices
        // fill up cluster properties (single face) for all vertices
        int id = 0;
        BOOST_FOREACH(face_descriptor face, vertices(dual))
        {
            Cluster_properties props;
           
            props.id = id++;
            props.concavity = 0;
            props.faces.push_back(face);

            BOOST_FOREACH(vertex_descriptor vert, vertices_around_face(halfedge(face, m_mesh), m_mesh))
            {
                props.conv_hull_pts.push_back(get(m_vpm, vert));
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

        const Cluster_properties& cluster_1_props = m_cluster_map[vert_1];
        const Cluster_properties& cluster_2_props = m_cluster_map[vert_2];

        decimation_props.new_cluster_props.id = -1;

        // faces 
        decimation_props.new_cluster_props.faces.resize(cluster_1_props.faces.size() + cluster_2_props.faces.size());

        // add faces from the first cluster
        std::copy(cluster_1_props.faces.begin(), cluster_1_props.faces.end(), decimation_props.new_cluster_props.faces.begin());
        // add faces from the second cluster
        std::copy(cluster_2_props.faces.begin(), cluster_2_props.faces.end(), decimation_props.new_cluster_props.faces.begin() + cluster_1_props.faces.size());

        // convex hull points
        std::vector<Point_3> common_hull_pts(cluster_1_props.conv_hull_pts.size() + cluster_2_props.conv_hull_pts.size()); // merged list of the lists of the convex hull points of two clusters
        
        // add the convex hull points of the first cluster
        std::copy(cluster_1_props.conv_hull_pts.begin(), cluster_1_props.conv_hull_pts.end(), common_hull_pts.begin());
        // add the convex hull points of the second cluster
        std::copy(cluster_2_props.conv_hull_pts.begin(), cluster_2_props.conv_hull_pts.end(), common_hull_pts.begin() + cluster_1_props.conv_hull_pts.size());

        // if can construct convex hull
        if (common_hull_pts.size() > 3) //TODO: add collinearity check
        {
            // compute convex hull on the merged list of points
            TriangleMesh conv_hull;
            CGAL::convex_hull_3(common_hull_pts.begin(), common_hull_pts.end(), conv_hull);

            // fill up the list of the convex hull points of produced cluster after decimation
            decimation_props.new_cluster_props.conv_hull_pts.clear();
            BOOST_FOREACH(vertex_descriptor vert, vertices(conv_hull))
            {
                decimation_props.new_cluster_props.conv_hull_pts.push_back(get(CGAL::vertex_point, conv_hull, vert));
            }

            // compute concavity value of the produced cluster using concavity calculator of the input mesh, faces, and convex hull of the produced cluster
            decimation_props.new_cluster_props.concavity = m_concavity_calc.compute(decimation_props.new_cluster_props.faces, conv_hull);
        }
        else
        {
            decimation_props.new_cluster_props.conv_hull_pts = boost::move(common_hull_pts);
            decimation_props.new_cluster_props.concavity = 0;
        }

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

        CGAL_assertion(vert_1 != vert_2);

        // id of the produced cluster after decimation
        int result_id = m_cluster_map[vert_1].id;

        // assign cluster properties of the produced cluster to the first vertex of the edge
        m_cluster_map[vert_1] = m_decimation_map[edge].new_cluster_props;
        m_cluster_map[vert_1].id = result_id;

        // connect the first vertex to all adjacent vertices of the second one except self
        BOOST_FOREACH(graph_vertex_descriptor vert, boost::adjacent_vertices(vert_2, m_graph))
        {
            if (vert == vert_1) continue; // backward edge
            if (m_cluster_map[vert].id == result_id) continue; // no need to add edge between the same cluster

            add_edge(vert_1, vert, m_graph);
        }

        // remove edges from candidates that are going to be removed
        BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_2, m_graph))
        {
            remove_candidate(edge);
        }
        
        // remove adjacent edges incident to the second vertex and remove the vertex
        clear_vertex(vert_2, m_graph);
        remove_vertex(vert_2, m_graph);

        // remove candidate edges with old value of decimation cost
#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
        int cnt = 0;
#endif
        BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_1, m_graph))
        {
            remove_candidate(edge); 
#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
                ++cnt;
#endif
        }

        // update decimation costs of all modified edges (only adjacent edges to the first vertex)

#ifdef CGAL_LINKED_WITH_TBB

        // functor that calls update_edge method
        struct Update_edge_functor
        {
            Update_edge_functor(Approx_decomposition& alg, double& concavity_threshold, double& alpha_factor)
            : m_alg(alg), m_concavity_threshold(concavity_threshold), m_alpha_factor(alpha_factor) {}

            void operator() (const graph_edge_descriptor& edge) const
            {
                m_alg.update_edge(edge, m_concavity_threshold, m_alpha_factor);
            }

        private:
            Approx_decomposition& m_alg;
            double& m_concavity_threshold;
            double& m_alpha_factor;
        };

        if (boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
        {
            Update_edge_functor update_functor(*this, concavity_threshold, alpha_factor);
            
            tbb::parallel_for_each(boost::out_edges(vert_1, m_graph).first, boost::out_edges(vert_1, m_graph).second, update_functor);
        }
        else
#endif
        {
            BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_1, m_graph))
            {
                update_edge(edge, concavity_threshold, alpha_factor);
            }
        }

        // add candidate edges with new value of decimation cost
        BOOST_FOREACH(graph_edge_descriptor edge, boost::out_edges(vert_1, m_graph))
        {
            add_candidate(edge, concavity_threshold);
        }
        
#ifdef CGAL_APPROX_DECOMPOSITION_VERBOSE
        std::cout << "Updated edges: " << cnt << std::endl;
#endif

        // remove edges that can not become candidates
        remove_invalid_edges();
    }

    /**
     * Compute aspect ratio of a cluster.
     */
    double compute_aspect_ratio(const std::vector<face_descriptor>& faces)
    {
        CGAL_assertion(!faces.empty());

        double perimeter = 0;
        double area = 0;

        // sum up area
        BOOST_FOREACH(face_descriptor face, faces)
        {
            area += CGAL::Polygon_mesh_processing::face_area(face, m_mesh);
        }

        // sum up lengths of border edges
        struct Perimeter_calculator 
        {
            Perimeter_calculator(const TriangleMesh& mesh, double& perimeter)
            : m_mesh(mesh), m_perimeter(perimeter)
            {}

            void operator() (const halfedge_descriptor& h) const
            {
                m_perimeter += CGAL::Polygon_mesh_processing::edge_length(edge(h, m_mesh), m_mesh);
            }
            
            const TriangleMesh& m_mesh;
            double& m_perimeter;
        };

        CGAL::Polygon_mesh_processing::border_halfedges(faces, m_mesh, boost::make_function_output_iterator(Perimeter_calculator(m_mesh, perimeter)));

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

    /**
     * Adds an edge to the candidates list if its concavity value satisfies the threshold, otherwise adds it to the list of edges to be removed
     */
    void add_candidate(graph_edge_descriptor edge, double concavity_threshold)
    {
        // if concavity value of the produced cluster doesn't satisfy the threshold then mark the edge as invalid (for further removal)
        if (m_decimation_map[edge].new_cluster_props.concavity > concavity_threshold)
        {
            m_invalid_edges.push_back(edge);
            return;
        }
        
        // otherwise the edge is a candidate for decimation operator
        m_candidates.insert(edge);
    }

    /**
     * Removes an edge from the candidates list
     */
    void remove_candidate(graph_edge_descriptor edge)
    {
        m_candidates.erase(edge);
    }

    /**
     * Removes edges from the adjacency list that were marked to be removed
     */
    void remove_invalid_edges()
    {
        BOOST_FOREACH(graph_edge_descriptor edge, m_invalid_edges)
        {
            remove_edge(edge, m_graph);
        }
        m_invalid_edges.clear();
    }
};    

}
}

#endif // CGAL_APPROX_DECOMPOSITION_H
