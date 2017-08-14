#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H


#include <fstream>
#include <math.h>
#include <utility>

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {



template<typename PolygonMesh, typename VertexPointMap, typename VertexConstraintMap, typename EdgeConstraintMap, typename GeomTraits>
class Compatible_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

    typedef typename GeomTraits::Point_3 Point;
    typedef typename GeomTraits::Vector_3 Vector;
    typedef typename GeomTraits::Segment_3 Segment;
    typedef typename GeomTraits::Triangle_3 Triangle;
    typedef std::vector<Triangle> Triangle_list;

    typedef CGAL::AABB_triangle_primitive<GeomTraits, typename Triangle_list::iterator> AABB_Primitive;
    typedef CGAL::AABB_traits<GeomTraits, AABB_Primitive> AABB_Traits;
    typedef CGAL::AABB_tree<AABB_Traits> Tree;

    // <for each halfedge around v, pair of incident halfedges to this halfedge around v>
    typedef std::pair<halfedge_descriptor, halfedge_descriptor> he_pair;
    typedef std::map<halfedge_descriptor, he_pair> Edges_around_map;


public:
    Compatible_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap, VertexConstraintMap& vcmap, EdgeConstraintMap& ecmap) :
        mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap), ecmap_(ecmap)
    {}

    ~Compatible_remesher()
    {
        delete tree_ptr_;
    }

    template<typename FaceRange>
    void init_remeshing(const FaceRange& face_range)
    {        
        check_constraints();

        BOOST_FOREACH(face_descriptor f, face_range)
        {
            // todo: avoid push back and reserve the space
            input_triangles_.push_back(triangle(f));
        }

        tree_ptr_ = new Tree(input_triangles_.begin(), input_triangles_.end());
        tree_ptr_->accelerate_distance_queries();

    }

    std::size_t remove_degenerate_faces()
    {
        std::size_t nb_removed_faces = 0;

        // from repair.h
        nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
        std::cout<<"nb_collapsed_faces: "<<nb_removed_faces<<std::endl;
#endif

        return nb_removed_faces;
    }

    void angle_relaxation(bool use_weights)
    {
        std::map<vertex_descriptor, Point> barycenters;
        std::map<vertex_descriptor, Vector> n_map;
        BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
        {

            if(!is_border(v, mesh_) && !is_constrained(v))
            {

                // compute normal to v
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(traits_));
                n_map[v] = vn;

                Edges_around_map he_map;
                typename Edges_around_map::iterator it;
                BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_source(v, mesh_))
                    he_map[hi] = he_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );

                // calculate movement
                Vector move = CGAL::NULL_VECTOR;
                double weights_sum = 0;
                for(it = he_map.begin(); it != he_map.end(); ++it)
                {
                    halfedge_descriptor main_he = it->first;
                    he_pair incident_pair = it->second;

                    Vector rotated_edge = rotate_edge(main_he, incident_pair);

                    // calculate angle
                    double angle = get_angle(incident_pair, vn);
                    if(angle < 1e-5)
                        continue;

                    // small angles carry more weight
                    double weight = 1.0 / (angle*angle);
                    weights_sum += weight;

                    if(use_weights)
                        move += weight * rotated_edge;
                    else
                        move += rotated_edge;
                }

                // if at least 1 angle was summed
                if(use_weights && weights_sum != 0)
                    move /= weights_sum;
                else
                    move /= CGAL::to_double(he_map.size());

                barycenters[v] = (get(vpmap_, v) + move) ;

            } // not on border
        } // for each v

        // compute locations on tangent plane
        typedef typename std::map<vertex_descriptor, Point>::value_type VP;
        std::map<vertex_descriptor, Point> new_locations;
        BOOST_FOREACH(const VP& vp, barycenters)
        {
            Point p = get(vpmap_, vp.first);
            Point q = vp.second;
            Vector n = n_map[vp.first];

            new_locations[vp.first] = q + ( n * Vector(q, p) ) * n ;
        }

        // update location
        BOOST_FOREACH(const VP& vp, new_locations)
        {

#ifdef CGAL_PMP_SMOOTHING_DEBUG
            std::cout << "from: "<< get(vpmap_, vp.first);
#endif
            put(vpmap_, vp.first, vp.second);

#ifdef CGAL_PMP_SMOOTHING_DEBUG
            std::cout<<" moved at: "<< vp.second << std::endl;
#endif
        }

    }

    void area_relaxation(const double& precision)
    {

#ifdef CGAL_PMP_SMOOTHING_DEBUG
        count_non_convex_energy_ = 0;
#endif

        unsigned int moved_points = 0;
        BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
        {
             if(!is_border(v, mesh_) && !is_constrained(v))
             {
                 if (gradient_descent(v, precision))
                 {
                     moved_points++;
                 }

             }

        }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
        std::cout<<"moved points: "<<moved_points<<" times"<<std::endl;
        std::cout<<"non_convex_energy found: "<<count_non_convex_energy_<<" times"<<std::endl;
#endif

    }

    void project_to_surface()
    {
        BOOST_FOREACH( vertex_descriptor v, vertices(mesh_))
        {
            if(!is_border(v, mesh_) && !is_constrained(v))
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

    double element_area(const vertex_descriptor& p1,
                        const vertex_descriptor& p2,
                        const vertex_descriptor& p3) const
    {
        return to_double(CGAL::approximate_sqrt(
                             traits_.compute_squared_area_3_object()(
                                   get(vpmap_, p1),
                                   get(vpmap_, p2),
                                   get(vpmap_, p3))));
    }

    double element_area(const Point& P,
                        const vertex_descriptor& p2,
                        const vertex_descriptor& p3) const
    {
        return to_double(CGAL::approximate_sqrt(
                             traits_.compute_squared_area_3_object()(
                                   P,
                                   get(vpmap_, p2),
                                   get(vpmap_, p3))));
    }

    double compute_average_area_around(const vertex_descriptor& v)
    {
        double sum_areas = 0;
        unsigned int number_of_edges = 0;

        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
        {
            // opposite vertices
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);

            double S = element_area(v, pi, pi1);
            sum_areas += S;
            number_of_edges++;
        }

        return sum_areas / number_of_edges;
    }

    double measure_energy(const vertex_descriptor& v, const double& S_av)
    {
        double energy = 0;
        unsigned int number_of_edges = 0;

        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);

            energy += (S - S_av)*(S - S_av);
            number_of_edges++;
        }

        return to_double( energy / number_of_edges );
    }

    double measure_energy(const vertex_descriptor& v, const double& S_av, const Point& new_P)
    {
        double energy = 0;
        unsigned int number_of_edges = 0;
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(new_P, pi, pi1);

            energy += (S - S_av)*(S - S_av);
            number_of_edges++;
        }

        return to_double( energy / (2 * number_of_edges) );
    }

    std::vector<double> calc_areas(const vertex_descriptor& v)
    {
        std::vector<double> areas;
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);
            areas.push_back(S);
        }

        return areas;
    }

    void compute_derivatives(double& drdx, double& drdy, double& drdz, const vertex_descriptor& v, const double& S_av)
    {
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);

            Vector vec(get(vpmap_, pi), get(vpmap_, pi1));

            // minimize r:
            // r = Σ(S-S_av)^2
            // dr/dx = 2 Σ(S - S_av) dS/dx
            // area of triangle with respect to (x_a, y_a, z_a) =
            // (1/2) [(v_z - v_y)x_a + (v_x - v_z)y_a + (v_y - v_x)z_a + constants]
            // vector v is (x_c - x_b, y_c - y_b, z_c - z_b)
            drdx += (S - S_av) * 0.5 * (vec.z() - vec.y());
            drdy += (S - S_av) * 0.5 * (vec.x() - vec.z());
            drdz += (S - S_av) * 0.5 * (vec.y() - vec.x());
        }

        drdx *= 2;
        drdy *= 2;
        drdz *= 2;
    }

    bool gradient_descent(const vertex_descriptor& v, const double& precision)
    {
        bool move_flag;
        double x, y, z, x_new, y_new, z_new, drdx, drdy, drdz;
        x = get(vpmap_, v).x();
        y = get(vpmap_, v).y();
        z = get(vpmap_, v).z();

        double S_av = compute_average_area_around(v);
        double energy = measure_energy(v, S_av);

        // if the adjacent areas are absolutely equal
        if(energy == 0)
            return false;

        double energy_new = 0;
        double relative_energy = precision + 1;
        unsigned int t = 1;
        double eta0 = 0.01;
        //double power_t = 0.25;
        double t0 = 0.001;
        double eta = eta0 / (1 + t0*t);

        while(relative_energy > precision)
        {
            drdx=0, drdy=0, drdz=0;
            compute_derivatives(drdx, drdy, drdz, v, S_av);

            x_new = x - eta * drdx;
            y_new = y - eta * drdy;
            z_new = z - eta * drdz;

            Point moved(x_new, y_new, z_new);
            energy_new = measure_energy(v, S_av, moved);

            if(energy_new < energy)
            {
                put(vpmap_, v, moved);
                move_flag = true;
            }
            else
            {

#ifdef CGAL_PMP_SMOOTHING_DEBUG
                count_non_convex_energy_++;
#endif
                return false;
            }

            relative_energy = CGAL::to_double( (energy - energy_new) / energy );

            // update
            x = x_new;
            y = y_new;
            z = z_new;
            energy = energy_new;
            t++;

            //eta = eta0 / pow(t, power_t);
            eta = eta0 / (1 + t0 * t);

        }

        return move_flag;
    }

    Vector rotate_edge(const halfedge_descriptor& main_he, const he_pair& incd_edges)
    {
        // get common vertex around which the edge is rotated
        Point pt = get(vpmap_, target(main_he, mesh_)); // CORRECT SYNATX vt

        // ps is the vertex that is being moved
        Point ps = get(vpmap_, source(main_he, mesh_));

        // get "equidistant" points - in fact they are at equal angles
        Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
        Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));
        CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

        Vector edge1(pt, equidistant_p1);
        Vector edge2(pt, equidistant_p2);
        Vector vec_main_he(pt, ps);

        // check degenerate cases
        double tolerance = 1e-3;

        if ( edge1.squared_length()           < tolerance ||
             edge2.squared_length()           < tolerance ||
             sqlength(main_he)                < tolerance ||
             (edge1 - vec_main_he).squared_length() < tolerance ||
             (edge2 - vec_main_he).squared_length() < tolerance   )
        {
            return CGAL::NULL_VECTOR;
        }

        CGAL_assertion(vec_main_he.squared_length() > tolerance);

        // find bisector
        Vector bisector = CGAL::NULL_VECTOR;
        internal::normalize(edge1, traits_);
        internal::normalize(edge2, traits_);
        bisector = edge1 + edge2;

        // under about 0.5 degrees deviation consider it flat
        if( bisector.squared_length() < 0.001 ) // -> length = 0.01 -> sin(theta) = 0.01 -> theta = 0.57 degrees
        {
            // angle is (almost) 180 degrees, take the perpendicular
            Vector normal_vec = find_perpendicular(edge1, pt, ps); // normal to edge and found on (vec_main_he)'s plane

            CGAL_assertion(normal_vec != CGAL::NULL_VECTOR);
            CGAL_assertion(CGAL::scalar_product(edge1, normal_vec) < tolerance);

            Segment b_segment(pt, pt + normal_vec);
            Point b_segment_end = b_segment.target();

            if(CGAL::angle(b_segment_end, pt, ps) == CGAL::OBTUSE)
            {
                b_segment = b_segment.opposite();
            }

            bisector = Vector(b_segment);
        }

        correct_bisector(bisector, main_he);

        double target_length = CGAL::sqrt(sqlength(main_he));
        double bisector_length = CGAL::sqrt(bisector.squared_length());

        CGAL_assertion(   ( target_length - tolerance    <   bisector_length     ) &&
                          ( bisector_length        <   target_length + tolerance )    );

        return bisector;
    }

    double get_angle(const he_pair& incd_edges, const Vector& vn)
    {
        Point s = get(vpmap_, source(incd_edges.first, mesh_));
        Point p1 = get(vpmap_, target(incd_edges.first, mesh_));
        Point p2 = get(vpmap_, source(incd_edges.second, mesh_));
        CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

        Vector v1(s, p1);
        Vector v2(s, p2);

        Vector cp = CGAL::cross_product(v1, v2);
        double det = CGAL::scalar_product(vn, cp);
        double dot = CGAL::scalar_product(v1, v2);

        // transform to range [0, 2pi]
        double res = atan2(-det, -dot) + CGAL_PI;

        return res;
    }


    Vector find_perpendicular(const Vector& input_vec, const Point& s, const Point& pv)
    {
        Vector s_pv(s, pv);
        Vector aux_normal = CGAL::cross_product(input_vec, s_pv);

        return CGAL::cross_product(aux_normal, input_vec);
    }

    Vector max_vector(const Vector& vec1, const Vector& vec2)
    {
        if (vec1.squared_length() > vec2.squared_length())
            return vec1;
        else
            return vec2;
    }

    void correct_bisector(Vector& bisector_vec, const halfedge_descriptor& main_he)
    {
        // get common vertex around which the edge is rotated
        Point pt = get(vpmap_, target(main_he, mesh_));

        // create a segment to be able to translate
        Segment bisector(pt, pt + bisector_vec);

        // scale
        double scale_factor = CGAL::sqrt(  sqlength(main_he) / bisector.squared_length() );
        typename GeomTraits::Aff_transformation_3 t_scale(CGAL::SCALING, scale_factor);
        bisector = bisector.transform(t_scale);

        // translate
        Vector vec(bisector.source(), pt);
        typename GeomTraits::Aff_transformation_3 t_translate(CGAL::TRANSLATION, vec);
        bisector = bisector.transform(t_translate);

        // take the opposite so that their sum is the overall displacement
        bisector_vec = -Vector(bisector);
    }

    bool is_constrained(const edge_descriptor& e)
    {
        return get(ecmap_, e);
    }

    bool is_constrained(const vertex_descriptor& v)
    {
        return get(vcmap_, v);
    }

    void check_constraints()
    {
        BOOST_FOREACH(edge_descriptor e, edges(mesh_))
        {
            if (is_constrained(e))
            {
                vertex_descriptor vs = source(e, mesh_);
                vertex_descriptor vt = target(e, mesh_);
                put(vcmap_, vs, true);
                put(vcmap_, vt, true);
            }
        }
    }




private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    VertexConstraintMap vcmap_;
    EdgeConstraintMap ecmap_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;
    GeomTraits traits_;

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    unsigned int count_non_convex_energy_;
#endif


};




} //namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL




#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H


