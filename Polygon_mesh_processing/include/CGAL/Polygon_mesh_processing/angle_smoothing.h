#ifndef ANGLE_SMOOTHING_H
#define ANGLE_SMOOTHING_H


#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {

namespace Polygon_mesh_processing {





template<typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
class Angle_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef CGAL::Halfedge_around_source_iterator<PolygonMesh> halfedge_around_source_iterator;


    typedef typename GeomTraits::Point_3 Point;
    typedef typename GeomTraits::Vector_3 Vector;
    typedef typename GeomTraits::Segment_3 Segment;



    // gather incident lines to a map
    // <one halfedge around v, pair of incident halfedges to this halfedge around v>
    typedef std::pair<halfedge_descriptor, halfedge_descriptor> He_pair;
    typedef std::map<halfedge_descriptor, He_pair> Edges_around_map;


public:
    Angle_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap)
    {}




    void angle_relaxation()
    {

        std::map<vertex_descriptor, Point> barycenters;
        boost::vector_property_map<Vector> n_map;

        unsigned int count = 0;
        for(vertex_descriptor v : vertices(mesh_))
        {


            if(!is_border(v, mesh_))
            {
                //std::cout<<"processing vertex: "<< v << std::endl;


                // compute normal to v
                Vector vn = compute_vertex_normal(v, mesh_,
                                                  Polygon_mesh_processing::parameters::vertex_point_map(vpmap_)
                                                  .geom_traits(GeomTraits()));
                n_map[v] = vn;



                Edges_around_map he_map;
                typename Edges_around_map::iterator it;

                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_)) // or make it around target
                    he_map[hi] = He_pair(next(hi, mesh_), prev(opposite(hi, mesh_) ,mesh_));


                // take a look
                /*
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
                */


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
            std::cout << "from: "<< get(vpmap_, vp.first);
            put(vpmap_, vp.first, vp.second);
            std::cout<<" moved at: "<< vp.second << std::endl;
        }




    }






private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;



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


    Vector rotate_edge(const halfedge_descriptor& main_he, const He_pair& incd_edges)
    {

        double tol = 1e-14;
        double precision = 1e-3;
        double magnifier = 1e+3;

        // get common vertex around which the edge is rotated
        Point s = get(vpmap_, target(main_he, mesh_));

        // get "equidistant" points - actualy ther are at equal angles
        Point equidistant_p1 = get(vpmap_, target(incd_edges.first, mesh_));
        Point equidistant_p2 = get(vpmap_, source(incd_edges.second, mesh_));

        Vector edge1(s, equidistant_p1);
        Vector edge2(s, equidistant_p2);

        if(edge1 != CGAL::NULL_VECTOR)
            internal::normalize(edge1, GeomTraits());
        if(edge2 != CGAL::NULL_VECTOR)
            internal::normalize(edge2, GeomTraits());

        // get bisector!
        Vector bisector = CGAL::NULL_VECTOR;

        // but be aware of stupid numerical errors - check possible overflow?
        bisector = edge1 * magnifier + edge2 * magnifier;


        // handle degenerate cases
        if(bisector.squared_length() < precision)
        {
            // pv is the vertex that is being moved
            Point pv = get(vpmap_, source(main_he, mesh_));

            if(s == pv)
            {
                // s, end_v, and equidistant points are almcollinear.
                CGAL_assertion(CGAL::collinear(s, pv, equidistant_p1));
                CGAL_assertion(CGAL::collinear(s, pv, equidistant_p2));
                return Vector(CGAL::NULL_VECTOR);
            }
            else
            {
                // s and equidistant points are collinear
                //CGAL_assertion(CGAL::collinear(s, equidistant_p1, equidistant_p2)); // any overload with some tolerance?
                Vector n(s, pv);
                bisector = n;
            }

        }

        correct_bisector(bisector, main_he);


        double target_length = CGAL::sqrt(sqlength(main_he));
        double bisector_length = CGAL::sqrt(bisector.squared_length());

        if( ! ( ( target_length - tol    <   bisector_length     ) &&
                ( bisector_length        <   target_length + tol ) ) )
        {
            std::cerr<<"problem";
        }

        CGAL_assertion(   ( target_length - tol    <   bisector_length     ) &&
                          ( bisector_length        <   target_length + tol )    );


        return bisector;


    }



    void correct_bisector(Vector& bisector_vec, const halfedge_descriptor& main_he)
    {

        // get common vertex around which the edge is rotated
        Point s = get(vpmap_, target(main_he, mesh_));

        // create a segment to be able to translate - avoid it?
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






};









template<typename PolygonMesh, typename NamedParameters>
void angle_remeshing(PolygonMesh& pmesh, const NamedParameters& np)
{



    //CGAL_PMP_NP_CLASS np;
    //np = CGAL::Polygon_mesh_processing::parameters::all_default();

    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Traits;


    CGAL::Polygon_mesh_processing::Angle_remesher<PolygonMesh, VertexPointMap, Traits> remesher(pmesh, vpmap);
    remesher.angle_relaxation();




}











} // namespace Polygon_mesh_processing
} //namespace CGAL






#endif // ANGLE_SMOOTHING_H
