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
  typedef typename Kernel::Line_3                                               Line_3;
  typedef typename Kernel::FT                                                   FT;
  typedef typename Kernel::RT                                                   RT;
  typedef typename Kernel::Vector_3                                             Vector_3;

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
  inline static Matrix4x4 combine_matrices(const Matrix4x4& aFirst, const Matrix4x4& aSecond) {
    return aFirst + aSecond;
  }

  /*
  * Returns true if the target of aHD is a discontinuity vertex in aTM
  *
  * Currently, only checks if vertex belongs to a border.
  */
  /*static bool is_discontinuity_vertex(const halfedge_descriptor& aHD, const TM& aTM) {
    if(is_border(target(aHD, aTM), aTM)) {
      return true;
    }
    return false;
  }*/

  static bool is_discontinuity_edge(const halfedge_descriptor& aHD, const TM& aTM) {
    return is_border_edge(aHD, aTM);
  }

  /*
  * fundamental error quadric for the target vertex of aHD in aTM
  */
  static Matrix4x4 fundamental_error_quadric(const halfedge_descriptor& aHD, const TM& aTM, FT aDiscontinuityMultiplier = 1.0) {
    Matrix4x4 quadric;
    quadric.setZero();

    vertex_descriptor target_vd = target(aHD, aTM);
    const Point target_vertex_point = std::move(get(boost::vertex_point, aTM, target_vd));
    const Vector_3 target_vertex_vector(target_vertex_point.x(),
                                        target_vertex_point.y(),
                                        target_vertex_point.z());

  //std::cout << "point: " << target_vertex_point.x() << " "
  //            << target_vertex_point.y() << " "
  //            << target_vertex_point.z() << std::endl;

    //bool discontinuity_vertex = is_discontinuity_vertex(aHD, aTM);

    int i = 0;
    for(const halfedge_descriptor hd: halfedges_around_target(target_vd, aTM)) {
      face_descriptor fd = face(hd, aTM);
      if(fd != GraphTraits::null_face()) {
        //std::cout << "face" << std::endl;
        std::vector<vertex_descriptor> vds;
        for(vertex_descriptor vd: vertices_around_face(hd, aTM)) {
          vds.push_back(vd);
        }
        Plane_3 plane(get(boost::vertex_point, aTM, vds[0]),
                      get(boost::vertex_point, aTM, vds[1]),
                      get(boost::vertex_point, aTM, vds[2]));

        Row4 plane_mtr;
        FT norm = sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        plane_mtr << plane.a() / norm, plane.b() / norm, plane.c() / norm, plane.d() / norm;
        quadric += plane_mtr.transpose() * plane_mtr;
        //std::cout << plane_mtr << std::endl;

        if(is_discontinuity_edge(hd, aTM)) {
          const vertex_descriptor source_vd = source(hd, aTM);
          const Vector_3 p1p2(target_vertex_point, get(boost::vertex_point, aTM, source_vd));
          const Vector_3 normal = cross_product(p1p2, plane.orthogonal_vector());
          FT d = - normal * target_vertex_vector;
          FT norm = sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
          plane_mtr << normal.x() / norm, normal.y() / norm, normal.z() / norm, d / norm;
          quadric += plane_mtr.transpose() * plane_mtr * aDiscontinuityMultiplier;
        }

        halfedge_descriptor shd = next(hd, aTM);
        if(is_discontinuity_edge(shd, aTM)) {
          const vertex_descriptor target_vd = target(shd, aTM);
          const Vector_3 p1p2(target_vertex_point, get(boost::vertex_point, aTM, target_vd));
          const Vector_3 normal = cross_product(p1p2, plane.orthogonal_vector());
          FT d = - normal * target_vertex_vector;
          FT norm = sqrt(normal.x() * normal.x() + normal.y() * normal.y() + normal.z() * normal.z());
          plane_mtr << normal.x() / norm, normal.y() / norm, normal.z() / norm, d / norm;
          quadric += plane_mtr.transpose() * plane_mtr * aDiscontinuityMultiplier;
        }

        /*if(discontinuity_vertex) {
          const vertex_descriptor source_vd = source(hd, aTM);
          const Line_3 edge_line = Line_3(get(boost::vertex_point, aTM, source_vd),
                                          target_vertex_point);

          Plane_3 plane = edge_line.perpendicular_plane(target_vertex_point);

          Row4 plane_mtr;
          FT norm = sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
          plane_mtr << plane.a() / norm, plane.b() / norm, plane.c() / norm, plane.d() / norm;
          std::cout << plane_mtr << std::endl;
          quadric += plane_mtr.transpose() * plane_mtr;
        }*/
      } else {
        //std::cout << "null face." << std::endl;
      }
    }

    //std::cout << std::endl;
    return quadric;
  }


  /*
  * Return the point p that minimizes p' Q p where p is free.
  * aP0, and aP1 are the points that are being collapsed.
  * aQuadric is the matrix that is the combination of matrices
  * of aP0 and aP1.
  */
  static Col4 optimal_point(const Matrix4x4& aQuadric, const Col4& aP0, const Col4& aP1) {
    Matrix4x4 X;
    X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

    Col4 opt_pt;

    if(X.determinant() == 0) {
      // not invertible
      Col4 p1mp0 = std::move(aP1 - aP0);
      FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0,0);
      FT b = 2 * (aP0.transpose() * aQuadric * p1mp0)(0,0);

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
          FT aP0_cost = (aP0.transpose() * aQuadric * aP0)(0,0);
          FT aP1_cost = (aP1.transpose() * aQuadric * aP1)(0,0);
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
      //Col4 rhs;
      //rhs << 0, 0, 0, 1;
      //opt_pt = X.inverse() * rhs;
      opt_pt = X.inverse().col(3);
    }
    return opt_pt;
  }

};


} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_CORE_H
