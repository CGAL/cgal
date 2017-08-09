#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Monge_via_jet_fitting.h>
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


    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap), cot_calculator_(pmesh, vpmap)
    {

        //std::vector<vertex_descriptor> points
       /* std::size_t num_points = vertices(mesh_).size();
        std::size_t min_num_of_points = 6; // (d+1)(d+2)/2, for d=2
        if(num_points < min_num_of_points)
        {
            CGAL_error_msg("Not enough points in the mesh.");
        }

        */

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


        //min_sq_edge_len_ = init_min_edge_length();

        //TODO: update constrained edges
        //check_constrained_edges();

    }


    // fixing degenerates

    std::size_t remove_degenerate_faces()
    {
        std::size_t nb_removed_faces = 0;

        // from repair.h
        nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);

        // debug
        std::cout<<"nb_collapsed_faces: "<<nb_removed_faces<<std::endl;

        return nb_removed_faces;
    }

    /*
    std::size_t collapse_short_edges()
    {
        std::size_t nb_collapsed_edges = 0;
        std::set<edge_descriptor> edges_to_collapse, non_topologically_valid_collapses;

        // collect short edges
        for(edge_descriptor e : edges(mesh_))
        {
            if(edge_should_collapse(e))
                edges_to_collapse.insert(e);
        }

        // collapse
        while(!edges_to_collapse.empty())
        {
            while(!edges_to_collapse.empty())
            {
                edge_descriptor e = *edges_to_collapse.begin();
                edges_to_collapse.erase(edges_to_collapse.begin());

                // link condition
                if(!Euler::does_satisfy_link_condition(e, mesh_))
                {
                    non_topologically_valid_collapses.insert(e);
                    continue;
                }

                // take out from short_edges prev and prev(opposite), which will be collapsed
                halfedge_descriptor he = halfedge(e, mesh_);

                // verbose - todo
                vertex_descriptor vs = source(he, mesh_);
                vertex_descriptor vt = target(he, mesh_);
                //

                edges_to_collapse.erase(edge(prev(he, mesh_), mesh_));
                edges_to_collapse.erase(edge(prev(opposite(he, mesh_), mesh_), mesh_));

                // shoot out edge to be collapsed from topologically_non_valid
                non_topologically_valid_collapses.erase(e);
                non_topologically_valid_collapses.erase(edge(prev(he, mesh_), mesh_));
                non_topologically_valid_collapses.erase(edge(prev(opposite(he, mesh_), mesh_), mesh_));

                vertex_descriptor v = Euler::collapse_edge(e, mesh_);

                // check if an out_edge now is too short
                for (edge_descriptor out_e : out_edges(v, mesh_))
                {
                    if(edge_should_collapse(out_e))
                        edges_to_collapse.insert(out_e);
                }

                nb_collapsed_edges++;
            }

            edges_to_collapse.swap(non_topologically_valid_collapses);
        }

        // debug
        std::cout<<"nb_collapsed_edges: "<<nb_collapsed_edges<<std::endl;

        return nb_collapsed_edges;
    }
    */

    void curvature_smoothing()
    {

        std::map<vertex_descriptor, Point> barycenters;
        std::map<vertex_descriptor, Vector> n_map;

        for(vertex_descriptor v : vertices(mesh_))
        {

            if(!is_border(v, mesh_)) // add && !is_constrained
            {


                // mean curvature
                //double k_mean = compute_mean_curvature();


                // normals

                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(GeomTraits())); //traits_



                /*
                // normalize (kn)
                /*
                Vector kn = k_mean * vn;
                normalize(kn, GeomTraits()); //traits_
                */

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

                    /*
                    // check if main_he is degenerate
                    double tolerance = 1e-3;
                    if(sqlength(hi)< tolerance)
                    {
                        continue;
                        // think of somehting better
                    }
                    */

                    // weight
                    std::pair<double, double> a1a2 = cot_angles(hi, it->second); // incd_edges
                    double weight = a1a2.first + a1a2.second;
                    sum_cot_weights += weight;

                    // if I return 0 weight, then vec becomes 0.

                    // displacement vector
                    Point Xi = get(vpmap_, source(hi, mesh_));
                    Point Xj = get(vpmap_, target(hi, mesh_));
                    Vector vec(Xj, Xi); // pointing to the point to be moved
                    // add weight
                    vec *= weight;
                    // sum vecs
                    weighted_move += vec;
                }

                // divide with total weight - be careful of a bloody division with zero.
                if(sum_cot_weights != 0) // risky business
                {
                    weighted_move /= sum_cot_weights;
                }


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
/*
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
*/
    std::pair<double, double> cot_angles(const halfedge_descriptor& main_he, const He_pair& incd_edges)
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
            return std::make_pair(0, 0);
            // zero means 90 degrees angle (and means also no weight)
        }

        CGAL_assertion(vec_main_he.squared_length() > tolerance);

        double a1 = cot_calculator_(vs, v1, vt);
        double a2 = cot_calculator_(vs, v2, vt);

        return std::make_pair(a1, a2);
    }
/*
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
*/
    // degenerate fixing functions
   /*
    bool edge_should_collapse(edge_descriptor e)
    {
        halfedge_descriptor he = halfedge(e, mesh_);
        Point s = get(vpmap_, source(he, mesh_));
        Point t = get(vpmap_, target(he, mesh_));

        double sq_length = traits_.compute_squared_distance_3_object()(s, t);

        if(sq_length < min_sq_edge_len_)
            return true;
        else
            return false;
    }

    struct Vertex_to_point
    {
        Vertex_to_point(VertexPointMap vpmap) : vpmap(vpmap){}

        typedef typename boost::property_traits<VertexPointMap>::reference output_type;

        output_type operator()(vertex_descriptor vd) const
        {
            return get(vpmap, vd);
        }

        VertexPointMap vpmap;
    };

    double init_min_edge_length()
    {
        vertex_iterator vi, ve;
        boost::tie(vi, ve) = vertices(mesh_);
        Vertex_to_point v_to_p(vpmap_);

        Bbox_3 bbox= CGAL::bbox_3(
                    boost::make_transform_iterator(vi, v_to_p),
                    boost::make_transform_iterator(ve, v_to_p));

        return 0.002 * diagonal_length(bbox);
    }

    double diagonal_length(const Bbox_3& bbox)
    {
      double dx = bbox.xmax() - bbox.xmin();
      double dy = bbox.ymax() - bbox.ymin();
      double dz = bbox.zmax() - bbox.zmin();

      double diag = dx * dx + dy * dy + dz * dz;
      return std::sqrt(diag);
    }

    */

    // data members
    // ------------
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    CGAL::internal::Cotangent_value_Meyer_secure<PolygonMesh, VertexPointMap> cot_calculator_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;
    GeomTraits traits_;
    double min_sq_edge_len_;


};




} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
