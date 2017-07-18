#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Monge_via_jet_fitting.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <utility>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {



template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
class Curvature_flow
{

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;


    typedef typename CGAL::Monge_via_jet_fitting<GeomTraits> Monge_via_jet_fitting;
    typedef typename Monge_via_jet_fitting::Monge_form Monge_form;

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


    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap),
                                                                cot_calculator_(mesh_, vpmap_)
    {

        //std::vector<vertex_descriptor> points
        std::size_t num_points = vertices(mesh_).size();
        std::size_t min_num_of_points = 6; // (d+1)(d+2)/2, for d=2
        if(num_points < min_num_of_points)
        {
            CGAL_error_msg("Not enough points in the mesh.");
        }


    }


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



    void curvature_smoothing()
    {

        std::map<vertex_descriptor, Point> barycenters;
        std::map<vertex_descriptor, Vector> n_map;

        for(vertex_descriptor v : vertices(mesh_))
        {

            if(!is_border(v, mesh_)) // add && !is_constrained
            {


                // mean curvature
                double k_mean = compute_mean_curvature();


                // normals
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(GeomTraits()));
                n_map[v] = vn;



                // find incident halfedges
                Edges_around_map he_map;
                typename Edges_around_map::iterator it;
                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
                    he_map[hi] = He_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );



                /*
                // find barycenter - without cot weights
                Vector displacement = CGAL::NULL_VECTOR;
                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
                {
                    displacement += Vector(get(vpmap_, target(hi, mesh_)) - get(vpmap_, source(hi, mesh_)));
                }
                barycenters[v] = get(vpmap_, v) + (displacement / halfedges_around_source(v, mesh_).size()) ;
                Point new_location = barycenters[v] + k_mean * n_map[v]; // point + vector
                */




                // with cotangent weights
                Vector displacement = CGAL::NULL_VECTOR;
                double sum_c = 0;

                for(it = he_map.begin(); it != he_map.end(); ++it)
                {
                    halfedge_descriptor main_he = it->first;
                    He_pair inc_pair = it->second;

                    std::pair<double, double> a1a2 = cot_angles(main_he, inc_pair);
                    double weight = a1a2.first + a1a2.second;
                    sum_c += weight;

                    displacement += Vector(get(vpmap_, source(main_he, mesh_)) - get(vpmap_, target(main_he, mesh_)));
                    displacement *= weight;
                }

                displacement /= sum_c;

                barycenters[v] = get(vpmap_, v) + (displacement / halfedges_around_source(v, mesh_).size()) ;
                Point new_location = barycenters[v] + k_mean * n_map[v]; // point + vector



                // update location
                put(vpmap_, v, new_location);




            } // in not on border


        } // all vertices




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


    double compute_mean_curvature()
    {
        // gather points around v
        std::vector<Point> incident_points = gather_all_points();

        Monge_form monge_form;
        Monge_via_jet_fitting monge_fit;

        std::size_t d_fit = 2; // d_fit >= d_monge
        std::size_t d_monge = 2; // need 2 principal coeeficients
        std::size_t Nd = (d_fit + 1)*(d_fit + 1) / 2.0;
        CGAL_assertion(incident_points.size() >= Nd);

        monge_form = monge_fit(incident_points.begin(), incident_points.end(), d_fit, d_monge);
        const FT k1 = monge_form.principal_curvatures(0);
        const FT k2 = monge_form.principal_curvatures(1);

        return (k1 + k2) / 2.0;
    }

    std::pair<double, double> cot_angles(const halfedge_descriptor& main_he, const He_pair& incd_edges)
    {
        vertex_descriptor v0 = source(main_he, mesh_);
        vertex_descriptor vs = target(main_he, mesh_);
        vertex_descriptor v1 = target(incd_edges.first, mesh_);
        vertex_descriptor v2 = source(incd_edges.second, mesh_);
        CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

        double a1 = cot_calculator_(v0, v1, vs);
        double a2 = cot_calculator_(v0, v2, vs);

        // to check degeneracies

        return std::make_pair(a1, a2);

    }

    std::vector<Point> points_around_vertex(vertex_descriptor v)
    {
        std::vector<Point> incident_vertices;
        for(halfedge_descriptor h : halfedges_around_target(v, mesh_))
        {
            vertex_descriptor vs = source(h, mesh_);
            incident_vertices.push_back(get(vpmap_, vs));
        }

        // temp assertion
        std::vector<Point> incident_vertices2;
        for(vertex_descriptor vi : vertices_around_target(v, mesh_))
        {
           incident_vertices2.push_back(get(vpmap_, vi));
        }
        CGAL_assertion(incident_vertices.size() == incident_vertices2.size());

        return incident_vertices;
    }


    std::vector<Point> gather_all_points()
    {
        std::vector<Point> points; // todo: preallocate it and fill it by pointing
        for(vertex_descriptor v : vertices(mesh_))
        {
            points.push_back(get(vpmap_, v));
        }

        return points;
    }

    Triangle triangle(face_descriptor f) const
    {
        halfedge_descriptor h = halfedge(f, mesh_);
        vertex_descriptor v1 = target(h, mesh_);
        vertex_descriptor v2 = target(next(h, mesh_), mesh_);
        vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
        return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
    }


    // data members
    // ------------
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    // Cotagent calculator class
    CGAL::internal::Cotangent_value_Meyer<
        PolygonMesh,
        typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type > cot_calculator_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;


};







} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
