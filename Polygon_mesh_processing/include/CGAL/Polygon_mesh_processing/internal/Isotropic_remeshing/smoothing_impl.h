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


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {




template<typename PolygonMesh, typename EdgeRange>
struct Edge_constraint_map
{

    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

    boost::shared_ptr<std::set<edge_descriptor>> constrained_edges_ptr;

public:
    typedef edge_descriptor                     key_type;
    typedef bool                                value_type;
    typedef value_type&                         reference;
    typedef boost::read_write_property_map_tag  category;

    Edge_constraint_map() : constrained_edges_ptr(new std::set<edge_descriptor>() )
    {}

    Edge_constraint_map(const EdgeRange& edges) : constrained_edges_ptr(new std::set<edge_descriptor>() )
    {
        for (edge_descriptor e : edges)
        {
            constrained_edges_ptr->insert(e);
        }
    }

    friend bool get(const Edge_constraint_map<PolygonMesh, EdgeRange>& map, const edge_descriptor& e)
    {
        // assertion on pmesh
        return map.constrained_edges_ptr->find(e) != map.constrained_edges_ptr->end();
    }

    friend void put(const Edge_constraint_map<PolygonMesh, EdgeRange>& map, const edge_descriptor& e)
    {
        map.constrained_edges_ptr->insert(e);
    }

    friend void remove(const Edge_constraint_map<PolygonMesh, EdgeRange>& map, const edge_descriptor& e)
    {
        map.constrained_edges_ptr->erase(e);
    }


};



template<typename PolygonMesh, typename VertexPointMap, typename VertexConstraitMap, typename EdgeConstraintMap, typename GeomTraits>
class Compatible_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator vertex_iterator;
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

    // <one halfedge around v, pair of incident halfedges to this halfedge around v>
    typedef std::pair<halfedge_descriptor, halfedge_descriptor> He_pair;
    typedef std::map<halfedge_descriptor, He_pair> Edges_around_map;


public:
    Compatible_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap, VertexConstraitMap& vcmap, EdgeConstraintMap& ecmap) :
        mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap), ecmap_(ecmap)
    {}

    ~Compatible_remesher()
    {
        delete tree_ptr_;
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

        //update constrained edges
        check_constrained_edges();

        min_edge_len_ = init_min_edge_length();
    }

    std::size_t remove_degenerate_faces()
    {
        std::size_t nb_removed_faces = 0;

        /*
        for(face_descriptor f : faces(mesh_))
        {
            Triangle tr = triangle(f);
            if(tr.is_degenerate())
            {
                halfedge_descriptor he = halfedge(f, mesh_);
                Euler::remove_face(he, mesh_);
                nb_removed_faces++;
            }
        }
        */

        // from repair.h
        nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);

        std::cout<<"nb_collapsed_faces: "<<nb_removed_faces<<std::endl;

        return nb_removed_faces;
    }

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

    void angle_relaxation(bool use_weights)
    {
        std::map<vertex_descriptor, Point> barycenters;
        //boost::vector_property_map<Vector> n_map;
        std::map<vertex_descriptor, Vector> n_map;

        for(vertex_descriptor v : vertices(mesh_))
        {

            if(!is_border(v, mesh_) && !is_constrained(v))
            {

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
std::cout<<"processing vertex: "<< v << std::endl;
#endif

                // compute normal to v
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(traits_));
                n_map[v] = vn;

                Edges_around_map he_map;
                typename Edges_around_map::iterator it;

                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_)) // or make it around target
                    he_map[hi] = He_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
for(it = he_map.begin(); it!=he_map.end(); ++it)
{
    halfedge_descriptor main_he = it->first;
    He_pair he_pair = it->second;
    std::cout<< "main: " << main_he;
    std::cout<<" ("<<source(main_he, mesh_)<<"->"<< target(main_he, mesh_)<<")";
    std::cout<< " - incident: "<< he_pair.first;
    std::cout<<" ("<<source(he_pair.first, mesh_)<<"->"<< target(he_pair.first, mesh_)<<")";
    std::cout<<" and " << he_pair.second;
    std::cout<<" ("<<source(he_pair.second, mesh_)<<"->"<< target(he_pair.second, mesh_)<<")"<<std::endl;
}
#endif

                // calculate movement
                Vector move = CGAL::NULL_VECTOR;
                double opposite_weight_factor = 0;

                for(it = he_map.begin(); it != he_map.end(); ++it)
                {
                    halfedge_descriptor main_he = it->first;
                    He_pair incident_pair = it->second;

                    Vector rotated_edge = rotate_edge(main_he, incident_pair);

                    // calculate angle
                    double angle = get_angle(incident_pair, vn);
                    if(angle < 1e-5) // temp
                        continue;

                    //small angles have more weight
                    double weight = 1.0 / (angle*angle);
                    opposite_weight_factor += weight;

                    if(use_weights)
                        move += weight * rotated_edge;
                    else
                        move += rotated_edge;
                }

                // if at least 1 angle was summed
                if(use_weights && opposite_weight_factor != 0)
                    move /= opposite_weight_factor;
                else
                    move /= CGAL::to_double(he_map.size());


                barycenters[v] = (get(vpmap_, v) + move) ;


            } // not on border
        } // for each v



        // compute locations on tangent plane
        typedef typename std::map<vertex_descriptor, Point>::value_type VP;
        std::map<vertex_descriptor, Point> new_locations;
        for(const VP& vp: barycenters)
        {
            Point q = vp.second;
            Point p = get(vpmap_, vp.first);
            Vector n = n_map[vp.first];

            new_locations[vp.first] = q + ( n * Vector(q, p) ) * n ;
        }


        // perform moves
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

    void area_relaxation(const double& precision)
    {

        count_non_convex_energy_ = 0; //temp;
        unsigned int moved_points = 0;

        for(vertex_descriptor v : vertices(mesh_))
        {

             if(!is_border(v, mesh_) && !is_constrained(v))
             {

                 if (gradient_descent(v, precision))
                 {
                     moved_points++;
                 }

             } // not on border


        }

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
        std::cout<<"moved points: "<<moved_points<<" times"<<std::endl;
        std::cout<<"non_convex_energy found: "<<count_non_convex_energy_<<" times"<<std::endl;
#endif

    }

    void project_to_surface()
    {
        for( vertex_descriptor v : vertices(mesh_))
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

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
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

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
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

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
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
        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
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
        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {
            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);

            Vector vec(get(vpmap_, pi), get(vpmap_, pi1));

            // minimize r:
            // r = Σ(S-S_av)^2
            // dr/dx = 2 Σ(S - S_av) dS/dx
            // area of triangle with respect to (x_a, y_a, z_a) =
            // (1/2) [(v_z - v_y)x_a + (v_x - v_z)y_a + (y_y - v_x)z_a + constants]
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

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
        std::ofstream out("data/energy.txt");
        std::ofstream out("data/coords.txt");
        std::ofstream out("data/areas.txt");
#endif
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

#ifdef CGAL_PMP_COMPATIBLE_REMESHING_DEBUG
            std::vector<double> areas = calc_areas(v);
            std::cout<<"v= "<<v<<std::endl;
            for(unsigned int i=0; i<areas.size(); ++i)
            {
                std::cout<<areas[i]<<"\t";
            }
            std::cout<<std::endl;
#endif

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
                count_non_convex_energy_++;
                return false;
            }

            relative_energy = to_double( (energy - energy_new) / energy );

            // update
            x = x_new;
            y = y_new;
            z = z_new;
            energy = energy_new;
            t++;

            //eta = eta0 / pow(t, power_t);
            eta = eta0 / (1 + t0*t);

        }

    return move_flag;

    }

    Vector rotate_edge(const halfedge_descriptor& main_he, const He_pair& incd_edges)
    {

        // get common vertex around which the edge is rotated
        Point s = get(vpmap_, target(main_he, mesh_));

        // pv is the vertex that is being moved
        Point pv = get(vpmap_, source(main_he, mesh_));

        // get "equidistant" points - in fact they are at equal angles
        Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
        Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));
        CGAL_assertion(target(incd_edges.second, mesh_) == source(incd_edges.first, mesh_));

        Vector edge1(s, equidistant_p1);
        Vector edge2(s, equidistant_p2);
        Vector s_pv(s, pv);

        // check degenerate cases
        double tolerance = 1e-3; // to think about it

        if ( edge1.squared_length()          < tolerance ||
             edge2.squared_length()          < tolerance ||
             sqlength(main_he)               < tolerance ||
             (edge1 - s_pv).squared_length() < tolerance ||
             (edge2 - s_pv).squared_length() < tolerance   )
        {
            return CGAL::NULL_VECTOR;
        }

        CGAL_assertion(s_pv.squared_length() > tolerance);

        // get bisector
        Vector bisector = CGAL::NULL_VECTOR;
        internal::normalize(edge1, traits_);
        internal::normalize(edge2, traits_);
        bisector = edge1 + edge2;


        // under about 0.5 degrees deviation consider it flat
        if( bisector.squared_length() < 0.001 ) // -> length = 0.01 -> sin(theta) = 0.01 -> theta = 0.57 degrees
        {

            // angle is (almost) 180 degrees, take the perpendicular
            Vector normal_vec = find_perpendicular(edge1, s, pv); // normal to edge and found on (s-pv)'s plane

            CGAL_assertion(normal_vec != CGAL::NULL_VECTOR);
            CGAL_assertion(CGAL::scalar_product(edge1, normal_vec) < tolerance);

            Segment b_segment(s, s + normal_vec);
            Point b_segment_end = b_segment.target();

            if(CGAL::angle(b_segment_end, s, pv) == CGAL::OBTUSE)
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

    double get_angle(const He_pair& incd_edges, const Vector& vn)
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
        double res = atan2(-det, -dot) + M_PI;

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
        Point s = get(vpmap_, target(main_he, mesh_));

        // create a segment to be able to translate
        Segment bisector(s, s + bisector_vec);

        // scale
        double scale_factor = CGAL::sqrt(  sqlength(main_he) / bisector.squared_length() );
        typename GeomTraits::Aff_transformation_3 t_scale(CGAL::SCALING, scale_factor);
        bisector = bisector.transform(t_scale);

        // translate
        Vector vec(bisector.source(), s);
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

    void check_constrained_edges()
    {
        for(edge_descriptor e : edges(mesh_))
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

    // degenerate fixing functions

    bool edge_should_collapse(edge_descriptor e)
    {
        halfedge_descriptor he = halfedge(e, mesh_);
        Point s = get(vpmap_, source(he, mesh_));
        Point t = get(vpmap_, target(he, mesh_));

        double sq_length = traits_.compute_squared_distance_3_object()(s, t);

        if(sq_length < min_edge_len_)
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

        return 0.01 * diagonal_length(bbox);
    }

    double diagonal_length(const Bbox_3& bbox)
    {
      double dx = bbox.xmax() - bbox.xmin();
      double dy = bbox.ymax() - bbox.ymin();
      double dz = bbox.zmax() - bbox.zmin();

      double diag = dx * dx + dy * dy + dz * dz;
      return std::sqrt(diag);
    }



private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    VertexConstraitMap vcmap_;
    EdgeConstraintMap ecmap_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;
    unsigned int count_non_convex_energy_;
    GeomTraits traits_;
    double min_edge_len_;



};












} //namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL




#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H


