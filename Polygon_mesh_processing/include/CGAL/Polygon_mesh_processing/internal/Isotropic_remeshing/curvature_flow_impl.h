#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Bbox_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <utility>
#include <math.h>


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {



template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
class Curvature_flow
{

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator vertex_iterator;

    typedef typename GeomTraits::Point_3  Point;
    typedef typename GeomTraits::Vector_3 Vector;
    typedef typename GeomTraits::FT       FT;
    typedef typename GeomTraits::Triangle_3 Triangle;
    typedef std::vector<Triangle> Triangle_list;

    // <one halfedge around v, pair of incident halfedges to this halfedge around v>
    typedef std::pair<halfedge_descriptor, halfedge_descriptor> He_pair;
    typedef std::map<halfedge_descriptor, He_pair> Edges_around_map;

    typedef CGAL::AABB_triangle_primitive<GeomTraits, typename Triangle_list::iterator> AABB_Primitive;
    typedef CGAL::AABB_traits<GeomTraits, AABB_Primitive> AABB_Traits;
    typedef CGAL::AABB_tree<AABB_Traits> Tree;



public:

    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap), cot_calculator_(pmesh, vpmap){}

    template<typename FaceRange>
    void init_remeshing(const FaceRange& face_range)
    {
        for (face_descriptor f : face_range)
        {
            input_triangles_.push_back(triangle(f));
        }

        tree_ptr_ = new Tree(input_triangles_.begin(), input_triangles_.end());
        tree_ptr_->accelerate_distance_queries();

        //TODO: update constrained edges
        //check_constrained_edges();

    }

    std::size_t remove_degenerate_faces()
    {
        std::size_t nb_removed_faces = 0;

        // from repair.h
        nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
        std::cout<<"nb_collapsed_faces: "<<nb_removed_faces<<std::endl;
#endif

        return nb_removed_faces;
    }

    void curvature_smoothing()
    {
        std::map<vertex_descriptor, Point> barycenters;
        std::map<vertex_descriptor, Vector> n_map;

        for(vertex_descriptor v : vertices(mesh_))
        {
            if(!is_border(v, mesh_)) // add && !is_constrained
            {
                // normals
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(traits_));
                n_map[v] = vn;

                // find incident halfedges
                Edges_around_map he_map;
                typename Edges_around_map::iterator it;
                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
                    he_map[hi] = He_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );

                // maybe use a seperate function for this
                // with cot weights
                Vector weighted_move = CGAL::NULL_VECTOR;
                double sum_cot_weights = 0;
                for(it = he_map.begin(); it!= he_map.end(); ++it)
                {
                    halfedge_descriptor hi = it->first; //main_he

                    // weight
                    double weight = cot_angles(hi, it->second); // incd_edges
                    sum_cot_weights += weight;

                    // displacement vector
                    Point Xi = get(vpmap_, source(hi, mesh_));
                    Point Xj = get(vpmap_, target(hi, mesh_));
                    Vector vec(Xi, Xj);

                    // add weight. If weight is 0, then vec becomes 0.
                    vec *= weight;
                    // sum vecs
                    weighted_move += vec;
                }

                // divide with total weight - if there is actually weight
                if(sum_cot_weights != 0)
                    weighted_move /= sum_cot_weights;

                Point weighted_barycenter = get(vpmap_, v) + weighted_move;
                barycenters[v] = weighted_barycenter;

            } // not on border

        } // all vertices


        // compute locations on tangent plane
        typedef typename std::map<vertex_descriptor, Point>::value_type VP;
        std::map<vertex_descriptor, Point> new_locations;
        for(const VP& vp: barycenters)
        {
            Point p = get(vpmap_, vp.first);
            Point q = vp.second;
            Vector n = n_map[vp.first];

            new_locations[vp.first] = q + ( n * Vector(q, p) ) * n ;
            //new_locations[vp.first] = q;
        }

        // update location
        for(const VP& vp : new_locations)
        {

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
            std::cout << "from: "<< get(vpmap_, vp.first);
#endif
            put(vpmap_, vp.first, vp.second);

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
            std::cout<<" moved at: "<< vp.second << std::endl;
#endif
        }

    }

    void project_to_surface()
    {
        for( vertex_descriptor v : vertices(mesh_))
        {
            if(!is_border(v, mesh_) ) // todo: && !is_constrained(v)
            {
                Point p_query = get(vpmap_, v);
                Point projected = tree_ptr_->closest_point(p_query);
                put(vpmap_, v, projected);
            }
        }
    }


private:
    Triangle triangle(face_descriptor f) const
    {
        halfedge_descriptor h = halfedge(f, mesh_);
        vertex_descriptor v1 = target(h, mesh_);
        vertex_descriptor v2 = target(next(h, mesh_), mesh_);
        vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
        return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
    }

    double sqlength(const vertex_descriptor& v1, const vertex_descriptor& v2) const
    {
        return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
    }

    double sqlength(const halfedge_descriptor& h) const
    {
      vertex_descriptor v1 = target(h, mesh_);
      vertex_descriptor v2 = source(h, mesh_);
      return sqlength(v1, v2);
    }

    double sqlength(const edge_descriptor& e) const
    {
      return sqlength(halfedge(e, mesh_));
    }

    double cot_angles(const halfedge_descriptor& main_he, const He_pair& incd_edges)
    {
        vertex_descriptor vs = source(main_he, mesh_);
        vertex_descriptor vt = target(main_he, mesh_);
        vertex_descriptor v1 = target(incd_edges.first, mesh_);
        vertex_descriptor v2 = source(incd_edges.second, mesh_);

        CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

        Point p1 = get(vpmap_, v1);
        Point p2 = get(vpmap_, v2);
        Point pt = get(vpmap_, vt);
        Point ps = get(vpmap_, vs);
        Vector edge1(pt, p1);
        Vector edge2(pt, p2);
        Vector vec_main_he(pt, ps);

        double tolerance = 1e-3;
        if ( edge1.squared_length()           < tolerance ||
             edge2.squared_length()           < tolerance ||
             sqlength(main_he)                < tolerance ||
             (edge1 - vec_main_he).squared_length() < tolerance ||
             (edge2 - vec_main_he).squared_length() < tolerance   )
        {
            return 0;
            // zero means 90 degrees angle (also means no weight)
        }

        CGAL_assertion(vec_main_he.squared_length() > tolerance);

        double a1 = cot_calculator_(vs, v1, vt);
        double a2 = cot_calculator_(vs, v2, vt);

        return a1 + a2;
    }

    // data members
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    //CGAL::internal::Cotangent_value_Meyer_secure<PolygonMesh, VertexPointMap> cot_calculator_;
    CGAL::internal::Cotangent_value_clamped_2<PolygonMesh, VertexPointMap> cot_calculator_;
    //CGAL::internal::Cotangent_value_clamped<PolygonMesh, VertexPointMap> cot_calculator_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;
    GeomTraits traits_;
    double min_sq_edge_len_;


};




} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
