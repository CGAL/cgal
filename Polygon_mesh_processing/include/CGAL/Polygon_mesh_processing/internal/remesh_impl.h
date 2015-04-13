
#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

  template<typename PolygonMesh>
  class Incremental_remesher
  {
  public:
    Incremental_remesher(PolygonMesh& pmesh)
      : mesh_(pmesh)
    {}

    void split_long_edges(const double& high)
    {
      ;
    }
    void collapse_short_edges(const double& low, const double& high)
    {
      ;
    }
    void equalize_valences()
    {
      ;
    }
    void tangential_relaxation()
    {
      ;
    }
    void project_to_surface()
    {
      ;
    }

  private:
    PolygonMesh& mesh_;

  };//end class Incremenal_remesher
}//end namespace internal
}//end namesapce Polygon_mesh_processing
}//end namesapce CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REMESH_IMPL_H
