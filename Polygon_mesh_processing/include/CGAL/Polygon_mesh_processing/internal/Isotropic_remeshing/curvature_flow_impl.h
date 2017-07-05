#ifndef CURVATURE_FLOW_IMPL_H
#define CURVATURE_FLOW_IMPL_H



namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {




template<typename PolygonMesh, typename VertexPointMap>
class Curvature_flow
{





public:


    Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap) : mesh_(pmesh), vpmap_(vpmap)
    {}




    // k = div n


    // cot a = v1 * v2 / v1 x v2




private:
    PolygonMesh& mesh_;
    VertexPointMap& vpmap_;

};







} // internal
} // Polygon_mesh_processing
} // CGAL





#endif // CURVATURE_FLOW_IMPL_H
