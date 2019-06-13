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


  static Col4 point_to_homogenous_column(const Point& pt) {
    return Col4(pt.x(), pt.y(), pt.z(), 1);
  }

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
      if(fd != GraphTraits::null_face()) {
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
    }
    return quidric;
  }


  /*
  * Return the point p that minimizes p' Q p where p is free.
  * aP0, and aP1 are the points that are being collapsed.
  * aQuidric is the matrix that is the combination of matrices
  * of aP0 and aP1.
  */
  static Col4 optimal_point(const Matrix4x4& aQuidric, const Col4& aP0, const Col4& aP1) {
    Matrix4x4 X;

    X.block(0, 0, 3, 4) = aQuidric.block(0,0,3,4);
    X.block(3,0,1,3).setZero();
    X(3,3) = 1;

    Col4 opt_pt;

    if(X.determinant() == 0) {
      // not invertible
      Col4 p1mp0 = std::move(aP1 - aP0);
      FT a = (p1mp0.transpose() * aQuidric * p1mp0)(0,0);
      FT b = (p1mp0.transpose() * aQuidric * aP0 + aP0.transpose() * aQuidric * p1mp0)(0,0);

      if(a == 0) {
        if(b < 0) {
          opt_pt = aP1;
        } else if (b == 0) {
          opt_pt = (aP0 + aP1) / 2;
        } else {
          opt_pt = aP0;
        }
      } else {
        FT ext_t = -b/(2*a);
        if(ext_t < 0 || ext_t > 1 || a > 0) {
          // one of endpoints
          FT aP0_cost = (aP0.transpose() * aQuidric * aP0)(0,0);
          FT aP1_cost = (aP1.transpose() * aQuidric * aP1)(0,0);
          if(aP0_cost > aP1_cost) {
            opt_pt = aP1;
          } else {
            opt_pt = aP0;
          }
        } else {
          // extremum of the parabola
          opt_pt = aP0 + (aP1 - aP0) * ext_t;
        }
      }
    } else {
      // invertible
      Col4 rhs;
      rhs.setZero();
      rhs(3) = 1;
      opt_pt = X.inverse() * rhs;
    }
    return opt_pt;
  }

};


} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
