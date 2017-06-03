#ifndef ANGLE_SMOOTHING_H
#define ANGLE_SMOOTHING_H


#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>


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


public:
    Angle_remesher(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap)
    {}




    void angle_relaxation()
    {

        std::map<vertex_descriptor, Point> barycenters;

        BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
        {

            if(!is_border(v, mesh_))
            {
                std::cout<<"processing vertex: "<< v << std::endl;

                // gather incident lines to a map
                // <halfedge around v, <incident halfedges to halfedge around v>>
                typedef std::pair<halfedge_descriptor, halfedge_descriptor> He_pair;
                typedef std::map<halfedge_descriptor, He_pair> Edges_around_map;
                Edges_around_map he_map;


                for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
                {
                    he_map[hi] = He_pair(next(hi, mesh_), prev(opposite(hi, mesh_) ,mesh_));
                    //std::cout<<"edges: "<<next(hi, mesh_)<<std::endl;
                }


                // take a look
                typename Edges_around_map::iterator it;
                for(it = he_map.begin(); it!=he_map.end(); ++it)
                {
                    halfedge_descriptor main_he = it->first;
                    He_pair he_pair = it->second;
                    std::cout<< "main: " << main_he;
                    std::cout<< " - incident: "<< he_pair.first <<" and " << he_pair.second <<std::endl;
                }



                Vector move = CGAL::NULL_VECTOR;

                for(it = he_map.begin(); it != he_map.end(); ++it)
                {

                    He_pair he_pair = it->second;
                    Point eq_d_p1 = get(vpmap_, target(he_pair.first, mesh_));
                    Point eq_d_p2 = get(vpmap_, source(he_pair.second, mesh_));
                    Point m = CGAL::midpoint(eq_d_p1, eq_d_p2);

                    Point s = get(vpmap_, source(he_pair.first, mesh_));
                    typename GeomTraits::Segment_3 bisector(s, m);

                    // scale
                    halfedge_descriptor main_he = it->first;
                    double len = CGAL::sqrt(sqlength(main_he));
                    typename GeomTraits::Aff_transformation_3 t_scale(CGAL::SCALING, len);
                    bisector = bisector.transform(t_scale);

                    // translate
                    Vector vec(bisector.source(), s);
                    typename GeomTraits::Aff_transformation_3 t_translate(CGAL::TRANSLATION, vec);
                    bisector = bisector.transform(t_translate);

                    move += Vector(bisector.target(), bisector.source());

                }


                barycenters[v] = get(vpmap_, v) + (move / (double)he_map.size());


            } // if is not on border
        } // for each v



        // perform moves
        typedef typename std::map<vertex_descriptor, Point>::value_type VP;
        for(const VP& vp : barycenters)
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
