#ifndef CGAL_CONCAVITY_H
#define CGAL_CONCAVITY_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <boost/foreach.hpp>

#include <fstream>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>

namespace CGAL
{
namespace internal
{

    template <class TriangleMesh, class GeomTraits>
    class Concavity
    {
        typedef typename GeomTraits::Point_3 Point_3;
        typedef typename GeomTraits::Vector_3 Vector_3;
        typedef typename GeomTraits::Ray_3 Ray_3;
        
        typedef CGAL::Surface_mesh<Point_3> SM;
        
        typedef typename boost::graph_traits<SM>::vertex_descriptor vertex_descriptor;

        typedef CGAL::Face_filtered_graph<SM> Filtered_graph;
        typedef CGAL::AABB_face_graph_triangle_primitive<SM> AABB_primitive;
        typedef CGAL::AABB_tree<CGAL::AABB_traits<GeomTraits, AABB_primitive>> AABB_tree;

        typedef boost::optional<AABB_tree::Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection; // TODO: fix the compilation error
    
    public:
        Concavity(const TriangleMesh& mesh, const GeomTraits& traits)
        : m_mesh(mesh)
        , m_traits(traits)
        {}

        template <class FacetPropertyMap>
        double compute(FacetPropertyMap facet_ids, std::size_t cluster_id)
        {
            Filtered_graph filtered_mesh(m_mesh, cluster_id, facet_ids);

            SM cluster;
            CGAL::copy_face_graph(filtered_mesh, cluster);
//            {
//                std::ofstream os("cluster_" + std::to_string(cluster_id) + ".off");
//                os << cluster;
//            }

            Concavity concavity(cluster, m_traits);
            return concavity.compute();
        }

        double compute()
        {
            CGAL_assertion(!CGAL::is_empty(m_mesh));

            SM conv_hull;
            std::vector<Point_3> pts;

            BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices(m_mesh))
            {
//                std::cout << typeid(vert).name() << " " << vert << " " << m_mesh.point(vert) << std::endl;
                pts.push_back(m_mesh.point(vert));
            }

            CGAL::convex_hull_3(pts.begin(), pts.end(), conv_hull); 
//            {
//                std::ofstream os("ch.off");
//                os << conv_hull;
//            }
            
            return compute(conv_hull);
        }

        double compute(const SM& conv_hull)
        {
            typedef std::map<vertex_descriptor, Vector_3> Normals_map;
            Normals_map normals_map;

            CGAL::Polygon_mesh_processing::compute_vertex_normals(m_mesh, boost::associative_property_map<Normals_map>(normals_map)); 
            
            AABB_tree tree(CGAL::faces(conv_hull).begin(), CGAL::faces(conv_hull).end(), conv_hull);

            double result = 0;

            BOOST_FOREACH(vertex_descriptor vert, CGAL::vertices(m_mesh))
            {
//                std::cout << m_mesh.point(vert) << " " << normals_map[vert] << " " << normals_map[vert].squared_length() << " -> ";
                Ray_3 ray(m_mesh.point(vert), normals_map[vert]);
                
                auto intersection = tree.first_intersection(ray);
                if (intersection)
                {
                    const Point_3* p =  boost::get<Point_3>(&(intersection->first));
                    if (p)
                    {
//                        std::cout << *p << std::endl;
                        result = std::max(result, CGAL::squared_distance(m_mesh.point(vert), *p));
                    }
                }
//                std::cout << std::endl;
            }

            return result;
        }

    private:
        const TriangleMesh& m_mesh;
        const GeomTraits& m_traits;
    };

}
}

#endif // CGAL_CONCAVITY_H
