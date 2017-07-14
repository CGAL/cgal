#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/Weights.h>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {



template<typename PolygonMesh, typename VertexPointMap>
class Curvature_flow
{










public:


    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap),
                                                                cot_calculator_(pmesh_, vpmap_)
    {}





    // k = div n






private:
    PolygonMesh& pmesh_;
    VertexPointMap& vpmap_;

    // Cotagent calculator class
    CGAL::internal::Cotangent_value_Meyer<
        PolygonMesh,
        typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type >   cot_calculator_;

};







} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
