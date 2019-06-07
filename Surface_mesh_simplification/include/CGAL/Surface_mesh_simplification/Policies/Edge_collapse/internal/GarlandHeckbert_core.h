#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Cartesian/MatrixC33.h>

#include <Eigen/Dense>
#include <limits>
#include <vector>
namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template<class TM_>
struct GarlandHeckbertCore
{
  typedef TM_                                                                   TM;
  typedef boost::graph_traits<TM>                                               GraphTraits;

  typedef typename GraphTraits::edges_size_type                                 size_type;
  typedef typename GraphTraits::vertex_descriptor                               vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                             halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor                                 face_descriptor;
  typedef typename boost::property_map<TM, CGAL::vertex_point_t>::type          Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type        Point;
  typedef typename Kernel_traits<Point>::Kernel                                 Kernel;
  typedef typename Kernel::Plane_3                                              Plane_3;
  typedef typename Kernel::FT                                                   FT;
  typedef typename Kernel::RT                                                   RT;

  typedef typename Eigen::Matrix<FT, 4, 4> Matrix4x4;
  typedef typename Eigen::Matrix<FT, 1, 4> Row4;
  typedef typename Eigen::Matrix<FT, 4, 1> Col4;


  typedef std::unordered_map<vertex_descriptor, Matrix4x4> garland_heckbert_map_type;


  /**
  * Combines two Q matrices.
  * It is simply the addition of two matrices
  */
  static Matrix4x4 combine_matrices(const Matrix4x4& aFirst, const Matrix4x4& aSecond) {
    return aFirst + aSecond;
  }

  /*
  * fundamental error quidric for the target vertex of aHD in aTM
  */
  static Matrix4x4 fundamental_error_quidric(const halfedge_descriptor& aHD, const TM& aTM) {
    Matrix4x4 quidric;
    quidric.setZero();

    for(face_descriptor fd: faces_around_target(aHD, aTM)) {
      halfedge_descriptor incident_hd = halfedge(fd, aTM);
      std::vector<vertex_descriptor> vds;
      for(vertex_descriptor vd: vertices_around_face(incident_hd, aTM)) {
        vds.push_back(vd);
      }
      Plane_3 plane(get(boost::vertex_point, aTM, vds[0]),
                    get(boost::vertex_point, aTM, vds[1]),
                    get(boost::vertex_point, aTM, vds[2]));
      Row4 plane_mtr;
      plane_mtr << plane.a(), plane.b(), plane.c(), plane.d();
    //  std::cout << plane_mtr << std::endl << std::endl;
      quidric += plane_mtr.transpose() * plane_mtr;
    }
    return quidric;
  }


  /*
  * Return the point p that minimizes p' Q p
  */
  static Col4 optimal_point(const Matrix4x4& quidric) {
    Matrix4x4 X;
    X.block(0, 0, 3, 4) = quidric.block(0,0,3,4);
    X.block(3,0,1,3).setZero();
    X(3,3) = 1;

    Col4 rhs;
    rhs.setZero();
    rhs(3) = 1;

    Col4 pt = X.inverse() * rhs;

    return pt;
  }

};


} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
