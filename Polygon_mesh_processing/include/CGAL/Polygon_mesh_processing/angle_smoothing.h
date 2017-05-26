#ifndef ANGLE_SMOOTHING_H
#define ANGLE_SMOOTHING_H


#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/Euler_operations.h>

//#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


namespace CGAL {

namespace Polygon_mesh_processing {



template<typename PolygonMesh>
class Angle_remesher
{

    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;


public:
    Angle_remesher(PolygonMesh& pmesh) : mesh_(pmesh)
    {}




    void angle_relaxation()
    {

        BOOST_FOREACH(vertex_descriptor v, vertices(mesh_))
        {

            if(!is_border(v, mesh_))
            {
                std::cout<<"vertex inside: ";
                std::cout<< v << std::endl;
            }

        }


    }






private:
    PolygonMesh& mesh_;

};





















} // namespace Polygon_mesh_processing
} //namespace CGAL






#endif // ANGLE_SMOOTHING_H
