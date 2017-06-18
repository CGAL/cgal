#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H


#include <fstream>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {






template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
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

    // <one halfedge around v, pair of incident halfedges to this halfedge around v>
    typedef std::pair<halfedge_descriptor, halfedge_descriptor> He_pair;
    typedef std::map<halfedge_descriptor, He_pair> Edges_around_map;


public:
    Compatible_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap)
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
    }

    void angle_relaxation()
    {

        std::map<vertex_descriptor, Point> barycenters;
        boost::vector_property_map<Vector> n_map;

        for(vertex_descriptor v : vertices(mesh_))
        {

            if(!is_border(v, mesh_))
            {

#ifdef CGAL_ANGLE_BASED_SMOOTHING_DEBUG
std::cout<<"processing vertex: "<< v << std::endl;
#endif

                // compute normal to v
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(GeomTraits()));
                n_map[v] = vn;



                Edges_around_map he_map;
                typename Edges_around_map::iterator it;

                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_)) // or make it around target
                    he_map[hi] = He_pair(next(hi, mesh_), prev(opposite(hi, mesh_) ,mesh_));



#ifdef CGAL_ANGLE_BASED_SMOOTHING_DEBUG
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

                for(it = he_map.begin(); it != he_map.end(); ++it)
                {
                    halfedge_descriptor main_he = it->first;
                    He_pair incident_pair = it->second;

                    Vector rotated_edge = rotate_edge(main_he, incident_pair);

                    move += rotated_edge;

                }

                barycenters[v] = get(vpmap_, v) + (move / (double)he_map.size());


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

#ifdef CGAL_ANGLE_BASED_SMOOTHING_DEBUG
std::cout << "from: "<< get(vpmap_, vp.first);
#endif
            put(vpmap_, vp.first, vp.second);

#ifdef CGAL_ANGLE_BASED_SMOOTHING_DEBUG
std::cout<<" moved at: "<< vp.second << std::endl;
#endif
        }

    }

    void area_relaxation()
    {

        count_non_convex_energy_ = 0; //temp;
        unsigned int moved_points = 0;

        for(vertex_descriptor v : vertices(mesh_))
        {

             if(!is_border(v, mesh_))
             {

                 if (gradient_descent(v))
                 {

#ifdef CGAL_AREA_BASED_REMESHING_DEBUG
std::cout<<"point moved from: "<<get(vpmap_, v);
#endif
                     Point p_query = get(vpmap_, v);
                     Point projected = tree_ptr_->closest_point(p_query);

                     //do the projection
                     put(vpmap_, v, projected);

#ifdef CGAL_AREA_BASED_REMESHING_DEBUG
std::cout<<" to after projection: "<<get(vpmap_, v)<<std::endl;
#endif
                     moved_points++;

                 }

             } // not on border


        }

        //std::cout<<"moved points: "<<moved_points<<" times"<<std::endl;
        //std::cout<<"non_convex_energy found: "<<count_non_convex_energy_<<" times"<<std::endl;

    }

    void project_to_surface()
    {
        for( vertex_descriptor v : vertices(mesh_))
        {

            // to check if is constrained

            if(!is_border(v, mesh_))
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
                             GeomTraits().compute_squared_area_3_object()(
                                   get(vpmap_, p1),
                                   get(vpmap_, p2),
                                   get(vpmap_, p3))));
    }

    double element_area(const Point& P,
                        const vertex_descriptor& p2,
                        const vertex_descriptor& p3) const
    {
        return to_double(CGAL::approximate_sqrt(
                             GeomTraits().compute_squared_area_3_object()(
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

    void compute_derivatives(double& dFdx, double& dFdy, double& dFdz, const vertex_descriptor& v, const double& S_av)
    {

        for(halfedge_descriptor h : halfedges_around_source(v, mesh_))
        {

            vertex_descriptor pi = source(next(h, mesh_), mesh_);
            vertex_descriptor pi1 = target(next(h, mesh_), mesh_);
            double S = element_area(v, pi, pi1);

            Vector vec(get(vpmap_, pi), get(vpmap_, pi1));

            dFdx += (S - S_av) * 0.5 * (vec.z() - vec.y());
            dFdy += (S - S_av) * 0.5 * (vec.x() - vec.z());
            dFdz += (S - S_av) * 0.5 * (vec.y() - vec.x());

        }

        dFdx *= 2;
        dFdy *= 2;
        dFdz *= 2;

    }

    bool gradient_descent(const vertex_descriptor& v)
    {

        double eta = 0.01; //learning rate
        bool move_flag;
        double x, y, z, x_new, y_new, z_new, dFdx, dFdy, dFdz;
        x = get(vpmap_, v).x();
        y = get(vpmap_, v).y();
        z = get(vpmap_, v).z();

        double S_av = compute_average_area_around(v);
        double energy = measure_energy(v, S_av);

        // if the adjacent areas are absolutely equal
        if(energy == 0)
            return false;
        double energy_new = 0;

        //std::ofstream out("data/energy.txt");
        //std::ofstream out("data/areas.txt");

        double criterion = to_double( (energy - energy_new) / energy );

        while( criterion  > 0.001 ) // make it a named parameter
        {

            dFdx=0, dFdy=0, dFdz=0;
            compute_derivatives(dFdx, dFdy, dFdz, v, S_av);

            /*
            std::vector<double> areas = calc_areas(v);
            for(unsigned int i=0; i<areas.size(); ++i)
            {
                out<<areas[i]<<"\t";
            }
            out<<std::endl;
            */


            x_new = x - eta * dFdx;
            y_new = y - eta * dFdy;
            z_new = z - eta * dFdz;

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

            criterion = to_double( (energy - energy_new) / energy );

            // update
            x = x_new;
            y = y_new;
            z = z_new;
            energy = energy_new;

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

        Vector edge1(s, equidistant_p1);
        Vector edge2(s, equidistant_p2);
        Vector s_pv(s, pv);


        // check degenerate cases
        double tolerance = 1e-5; // to think about it

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
        internal::normalize(edge1, GeomTraits());
        internal::normalize(edge2, GeomTraits());
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





private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;
    Triangle_list input_triangles_;
    Tree* tree_ptr_;
    unsigned int count_non_convex_energy_;



};












} //namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL




#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_IMPL_H


