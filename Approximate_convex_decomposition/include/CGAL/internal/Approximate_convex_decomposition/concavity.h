#ifndef CGAL_CONCAVITY_H
#define CGAL_CONCAVITY_H

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <tbb/parallel_for_each.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

namespace CGAL
{
namespace internal
{

    template <class TriangleMesh, class GeomTraits, class ConcurrencyTag, class Mesh = TriangleMesh>
    class Concavity
    {
        // predefined structs
        struct Intersection_functor;

        // typedefs
        typedef typename GeomTraits::Point_3 Point_3;
        typedef typename GeomTraits::Vector_3 Vector_3;
        typedef typename GeomTraits::Ray_3 Ray_3;
        
        typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
        typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

        typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;
        
        typedef std::map<vertex_descriptor, Vector_3> Normals_map;

        typedef CGAL::Face_filtered_graph<TriangleMesh> Filtered_graph;

        typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> AABB_primitive;
        typedef CGAL::AABB_tree<CGAL::AABB_traits<GeomTraits, AABB_primitive>> AABB_tree;

        typedef boost::optional<typename AABB_tree::template Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection;
    
    public:
        Concavity(const TriangleMesh& mesh, const GeomTraits& traits)
        : m_mesh(mesh)
        , m_traits(traits)
        {}

        /**
         * Computes concavity value of a cluster of the mesh which id is specified.
         */
        template <class FacetPropertyMap>
        double compute(FacetPropertyMap facet_ids, std::size_t cluster_id)
        {
            Filtered_graph filtered_mesh(m_mesh, cluster_id, facet_ids);

            Concavity<Filtered_graph, GeomTraits, ConcurrencyTag, TriangleMesh> concavity(filtered_mesh, m_traits);
            return concavity.compute();
        }

        /**
         * Computes concavity value of the whole mesh.
         */
        double compute()
        {
            CGAL_assertion(!CGAL::is_empty(m_mesh));

            Mesh conv_hull;
            std::vector<Point_3> pts;

            if (num_vertices(m_mesh) <= 3) return 0;

            // extract the list points of the mesh
            BOOST_FOREACH(vertex_descriptor vert, vertices(m_mesh))
            {
                pts.push_back(get(CGAL::vertex_point, m_mesh)[vert]);
            }

            // compute convex hull
            CGAL::convex_hull_3(pts.begin(), pts.end(), conv_hull); 
            
            return compute(vertices(m_mesh), conv_hull);
        }

        /**
         * Constructs list of vertices from the list of faces and computes concavity value with the convex hull provided.
         * Faces list is a subset of all faces in the mesh.
         */
        double compute(const std::vector<face_descriptor>& faces, const Mesh& conv_hull)
        {
            boost::unordered_set<vertex_descriptor> pts;

            BOOST_FOREACH(face_descriptor face, faces)
            {
                BOOST_FOREACH(vertex_descriptor vert, vertices_around_face(halfedge(face, m_mesh), m_mesh))
                {
                    pts.insert(vert);
                }
            }

            return compute(std::make_pair(pts.begin(), pts.end()), conv_hull);
        }

        /**
         * Computes concavity value projecting vertices from a list onto a convex hull.
         * Vertices list a subset of all vertices in the mesh.
         */
        template <class iterator>
        double compute(const std::pair<iterator, iterator>& verts, const Mesh& conv_hull)
        {
            // compute normals if normals are not computed
            compute_normals();

            // construct AABB for fast computations of intersections between ray and convex hull
            AABB_tree tree(faces(conv_hull).begin(), faces(conv_hull).end(), conv_hull);

            // compute intersections and select the largest projection length
            double result = 0;

            // functor that computes intersection, its projection length from a vertex and maximizes the result variable
            struct Intersection_functor
            {
                Intersection_functor(const TriangleMesh& mesh, const Normals_map& normals_map, const AABB_tree& tree, double& result)
                : m_mesh(mesh), m_normals_map(normals_map), m_tree(tree), m_result(result) {}

                void operator() (const vertex_descriptor& vert) const
                {
                    Point_3 origin = get(CGAL::vertex_point, m_mesh)[vert];
                    Ray_3 ray(origin, m_normals_map.at(vert));
                    
                    Ray_intersection intersection = m_tree.first_intersection(ray);
                    if (intersection)
                    {
                        const Point_3* intersection_point =  boost::get<Point_3>(&(intersection->first));
                        if (intersection_point)
                        {
                            m_result = std::max(m_result, CGAL::squared_distance(origin, *intersection_point));
                        }
                    }
                }

            private:
                const TriangleMesh& m_mesh;
                const Normals_map& m_normals_map;
                const AABB_tree& m_tree;
                double& m_result;
            };

            Intersection_functor intersection_functor(m_mesh, m_normals_map, tree, result);

#ifdef CGAL_LINKED_WITH_TBB
            if (boost::is_convertible<ConcurrencyTag, Parallel_tag>::value)
            {
                tbb::parallel_for_each(verts.first, verts.second, intersection_functor);
            }
            else
#endif
            {
                BOOST_FOREACH(vertex_descriptor vert, verts)
                {
                    intersection_functor(vert);
                }
            }

            return CGAL::sqrt(result);
        }

    private:
        const TriangleMesh& m_mesh;
        const GeomTraits& m_traits;

        Normals_map m_normals_map;
        bool m_normals_computed = false;

        void compute_normals()
        {
            if (m_normals_computed) return; // if the normals are already computed, then skip

            CGAL::Polygon_mesh_processing::compute_vertex_normals(m_mesh, boost::associative_property_map<Normals_map>(m_normals_map));
            m_normals_computed = true;
        }
    };
}
}

#endif // CGAL_CONCAVITY_H
