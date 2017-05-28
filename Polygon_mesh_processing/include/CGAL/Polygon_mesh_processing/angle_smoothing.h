#ifndef ANGLE_SMOOTHING_H
#define ANGLE_SMOOTHING_H

#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Line_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_2.h>

//#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
//#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>




namespace CGAL {

namespace Polygon_mesh_processing {



template<typename PolygonMesh>
class Angle_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    typedef CGAL::Halfedge_around_source_iterator<PolygonMesh> halfedge_around_source_iterator;

    typedef CGAL::Line_2<PolygonMesh> Line;
    typedef CGAL::Point_2<PolygonMesh> Point;






public:
    Angle_remesher(PolygonMesh& pmesh) : mesh_(pmesh)
    {}




    void angle_relaxation()
    {

        BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
        {

            if(!is_border(v, mesh_))
            {
                std::cout<<"processing vertex: "<< v << std::endl;

                // gather lines adjacent to angles
                halfedge_around_source_iterator hi, he;
                halfedge_descriptor hnext;
                std::vector<halfedge_descriptor> he_lines; //lines stored here
                for(boost::tie(hi, he) = halfedges_around_source(v, mesh_); hi != he; ++hi)
                {
                    hnext = next(*hi, mesh_);
                    he_lines.push_back(hnext); //TODO: avoid push_back
                }

                // take a look
                vertex_descriptor vs, vt;
                for(int i=0; i<he_lines.size(); ++i)
                {
                    std::cout<<he_lines[i]<<std::endl;
                    vs = source(he_lines[i], mesh_);
                    vt = target(he_lines[i], mesh_);
                    std::cout<<"vs= "<<vs<<"   vt= "<<vt<<std::endl;
                    std::cout<<"Pvs= "<< mesh_.point(vs) <<"   Pvt= "<< mesh_.point(vt) <<std::endl;
                }

                // TODO: run it for each pair of halfedges
                Line bs = bisector(he_lines[0], he_lines[1]);

                double radius = sqlength(v, source(he_lines[0], mesh_));





            }

        }


    }






private:
    PolygonMesh& mesh_;
    //VertexPointMap& vpmap_;



    Line bisector(halfedge_descriptor h1, halfedge_descriptor h2) const
    {
        Line l1 (mesh_.point(source(h1, mesh_)), mesh_.point(target(h1, mesh_)));
        Line l2 (mesh_.point(source(h2, mesh_)), mesh_.point(target(h2, mesh_)));
        return CGAL::bisector(l1, l2);
    }


    double sqlength(const vertex_descriptor& v1, const vertex_descriptor& v2) const
    {
        // TO FIX with vpmap
        //return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
        return to_double(CGAL::squared_distance(mesh_.point(v1), mesh_.point(v2)));
    }



};





















} // namespace Polygon_mesh_processing
} //namespace CGAL






#endif // ANGLE_SMOOTHING_H
