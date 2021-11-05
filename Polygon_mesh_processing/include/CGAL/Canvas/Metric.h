// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//
// This class provides a definition of a metric by three orthogonal eigen-vectors.
// The magnitude of the vectors are the squared eigen value.
//
// For example, (3, 0, 0), (2, 0, 0), (1, 0, 0) denote that the space is elongated
// in (1, 0, 0) three times, in (2, 0, 0) two-times.
//
// Denote the three vectors V1, V2 and V3, and the eigen vectors v1, v2 and v3,
// squared eigen-values of the tensor metric M is e1, e2, e3, we have:
//      V1 = v1 * e1, V2 = v2 * e2, V3 = v3 * e3
//      R = [v1, v2, v3]
//      F_M = R diag{e1, e2, e3} R^{-1}
//      M = F_M * F_M^t

#ifndef CGAL_CANVAS_METRIC_H
#define CGAL_CANVAS_METRIC_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Triangle_3.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <Eigen/Dense>
#else
#error
#endif

#include <boost/algorithm/minmax_element.hpp>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL {
namespace Canvas {

template<typename K, typename KExact = K>
class Metric_base
{
public:
  typedef typename K::FT                      FT;
  typedef typename K::Point_3                 Point_3;
  typedef typename K::Vector_3                Vector_3;

private:
  Eigen::Matrix3d eigen_transformation, eigen_inverse_transformation;
  Eigen::Matrix3d full_matrix;
  mutable double e_max, e_min, e_n;
  Vector_3 v_max, v_min, v_n;

public:
  const Eigen::Matrix3d& get_transformation() const { return eigen_transformation; }
  const Eigen::Matrix3d& get_inverse_transformation() const { return eigen_inverse_transformation; }
  const Vector_3& get_vmin() const { return v_min; }
  const Vector_3& get_vmax() const { return v_max; }
  const Vector_3& get_vn() const { return v_n; }
  void get_min_eigenvector(Vector_3& v) const { v = v_min; }
  void get_max_eigenvector(Vector_3& v) const { v = v_max; }
  void get_third_eigenvector(Vector_3& v) const { v = v_n; }
  double get_min_eigenvalue() const { return e_min; }
  double get_max_eigenvalue() const { return e_max; }
  double get_third_eigenvalue() const { return e_n; }
  void set_min_eigenvalue(double emin) const { e_min = emin; }
  void set_max_eigenvalue(double emax) const { e_max = emax; }
  void set_third_eigenvalue(double en) const { e_n = en; }

  const Eigen::Matrix3d& get_mat() const
  {
    return full_matrix;
  }

public:
  // Transform
  template <typename Object>
  Object transform(const Object& p) const
  {
    Eigen::Vector3d ep(p.x(), p.y(), p.z());
    ep = eigen_transformation * ep;
    Object result(ep[0],ep[1],ep[2]);
    return result;
  }

  /* transformation_exact doesn't exist anymore
#ifdef ANISO_USE_EXACT
  template<typename Object>
  Object transform(const Object &p, KExact k_exact) const {
    return p.transform(transformation_exact);
  }
#endif
  */

  // Inverse transform
  template <typename Object>
  Object inverse_transform(const Object &p) const
  {
    Eigen::Vector3d ep(p.x(),p.y(),p.z());
    ep = eigen_inverse_transformation * ep;

    Object result(ep[0], ep[1], ep[2]);
    return result;
  }

  template <typename Kernel>
  CGAL::Triangle_3<Kernel> inverse_transform(const CGAL::Triangle_3<Kernel>& t) const
  {
    typename Kernel::Point_3 p = t.vertex(0);
    typename Kernel::Point_3 q = t.vertex(1);
    typename Kernel::Point_3 r = t.vertex(2);
    Eigen::Vector3d ep(p.x(), p.y(), p.z());
    ep = eigen_inverse_transformation * ep;
    p = Point_3(ep[0], ep[1], ep[2]);

    ep = Eigen::Vector3d(q.x(), q.y(), q.z());
    ep = eigen_inverse_transformation * ep;
    q = Point_3(ep[0], ep[1], ep[2]);

    ep = Eigen::Vector3d(r.x(), r.y(), r.z());
    ep = eigen_inverse_transformation * ep;
    r = Point_3(ep[0], ep[1], ep[2]);

    return CGAL::Triangle_3<Kernel>(p, q, r);
  }

  Bbox_3 inverse_transform(const Bbox_3& bb) const
  {
    double x[8],y[8],z[8];
    Eigen::Vector3d ev(bb.xmin(), bb.ymin(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[0] = ev[0]; y[0] = ev[1]; z[0] = ev[2];

    ev = Eigen::Vector3d(bb.xmin(), bb.ymin(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[1] = ev[0]; y[1] = ev[1]; z[1] = ev[2];

    ev = Eigen::Vector3d(bb.xmin(), bb.ymax(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[2] = ev[0]; y[2] = ev[1]; z[2] = ev[2];

    ev = Eigen::Vector3d(bb.xmin(), bb.ymax(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[3] = ev[0]; y[3] = ev[1]; z[3] = ev[2];

    ev = Eigen::Vector3d(bb.xmax(), bb.ymin(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[4] = ev[0]; y[4] = ev[1]; z[4] = ev[2];

    ev = Eigen::Vector3d(bb.xmax(), bb.ymin(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[5] = ev[0]; y[5] = ev[1]; z[5] = ev[2];

    ev = Eigen::Vector3d(bb.xmax(), bb.ymax(),bb.zmin());
    ev = eigen_inverse_transformation * ev;
    x[6] = ev[0]; y[6] = ev[1]; z[6] = ev[2];

    ev = Eigen::Vector3d(bb.xmax(), bb.ymax(),bb.zmax());
    ev = eigen_inverse_transformation * ev;
    x[7] = ev[0]; y[7] = ev[1]; z[7] = ev[2];

    std::pair<double*,double*> xmm = boost::minmax_element(x, x+8);
    std::pair<double*,double*> ymm = boost::minmax_element(y, y+8);
    std::pair<double*,double*> zmm = boost::minmax_element(z, z+8);
    return Bbox_3(*(xmm.first), *(ymm.first), *(zmm.first),
                  *(xmm.second), *(ymm.second), *(zmm.second));
  }

  FT compute_distortion(const Metric_base &m) const
  {
    double eigen_dis1 = (m.eigen_transformation * eigen_inverse_transformation).operatorNorm();
    double eigen_dis2 = (eigen_transformation * m.eigen_inverse_transformation).operatorNorm();

#ifdef ANISO_DEBUG_METRIC
    double max_dist = 1.001;
    if(max(eigen_dis1, eigen_dis2) > max_dist)
    {
      std::cout << "--------" << std::endl;
      std::cout << "this evs : " << get_third_eigenvalue() << " " << get_min_eigenvalue() << " " << get_max_eigenvalue() << std::endl;
      std::cout << "m evs : " << m.get_third_eigenvalue() << " " << m.get_min_eigenvalue() << " " << get_max_eigenvalue() << std::endl;
      std::cout << "--THIS-ET---" << std::endl;
      std::cout << eigen_transformation << std::endl;
      std::cout << "--THIS-EIT---" << std::endl;
      std::cout << eigen_inverse_transformation << std::endl;
      std::cout << "--M-ET---" << std::endl;
      std::cout << m.eigen_transformation << std::endl;
      std::cout << "--M-EIT---" << std::endl;
      std::cout << m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--M-ET-*-THIS-EIT-" << std::endl;
      std::cout << m.eigen_transformation * eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--THIS-ET-*-M-EIT" << std::endl;
      std::cout << eigen_transformation * m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "vals dist1,dist2 : " << eigen_dis1 << " " << eigen_dis2 << std::endl;
    }
#endif

    return max(eigen_dis1, eigen_dis2);
  }

  void construct(const Vector_3& axis_x, //normal
                 const Vector_3& axis_y,
                 const Vector_3& axis_z,
                 const double& vpn,
                 const double& vpmax, //curvature max
                 const double& vpmin, //curvature min
                 const FT& epsilon)
  {
    e_n = (std::max)(epsilon, std::abs(vpn));
    e_max = (std::max)(epsilon, std::abs(vpmax)); //vpmax
    e_min = (std::max)(epsilon, std::abs(vpmin)); //vpmin

    v_n = axis_x;
    v_max = axis_y;
    v_min = axis_z;

    Eigen::Matrix3d eigen_m;
    eigen_m(0,0) = v_n.x();  eigen_m(0,1) = v_max.x();  eigen_m(0,2) = v_min.x();
    eigen_m(1,0) = v_n.y();  eigen_m(1,1) = v_max.y();  eigen_m(1,2) = v_min.y();
    eigen_m(2,0) = v_n.z();  eigen_m(2,1) = v_max.z();  eigen_m(2,2) = v_min.z();

    Eigen::Matrix3d eigen_diag = Eigen::Matrix3d::Zero();
    eigen_diag(0,0) = e_n;  //should be max(e_max)
    eigen_diag(1,1) = e_max;
    eigen_diag(2,2) = e_min;

    Eigen::Matrix3d eigen_mtransp = eigen_m.transpose();
    eigen_transformation = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "--- eigen_m ---" << std::endl;
    std::cout << eigen_m << std::endl;
    std::cout << "--- eigen_diag ---" << std::endl;
    std::cout << eigen_diag << std::endl;
    std::cout << "-- eigen_mtransp ---" << std::endl;
    std::cout << eigen_mtransp << std::endl;
    std::cout << "-- eigen_transformation ---" << std::endl;
    std::cout << eigen_transformation << std::endl;
#endif

    eigen_diag(0,0) = 1.0/eigen_diag(0,0);
    eigen_diag(1,1) = 1.0/eigen_diag(1,1);
    eigen_diag(2,2) = 1.0/eigen_diag(2,2);
    eigen_inverse_transformation = eigen_m * eigen_diag * eigen_mtransp;

    full_matrix = eigen_transformation.transpose()*eigen_transformation;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "-- eigen_inverse_transformation ---" << std::endl;
    std::cout << eigen_inverse_transformation << std::endl;
#endif
  }

public:
  friend
  std::ostream& operator<<(std::ostream& out, const Metric_base& x)
  {
    out << "M  = " << x.eigen_transformation << std::endl;
    return out;
  }

  Metric_base()
  {
    construct(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1, 1, 1, 1e-6);
  }

  Metric_base(Eigen::Matrix3d eigen_transformation_)
  {
    //std::cout << "building a metric from a matrix3d with input : " << std::endl;
    //std::cout << eigen_transformation_ << std::endl;

    construct(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 1, 1, 1, 1e-6);
    eigen_transformation = eigen_transformation_;

    bool invertible;
    double determinant;
    eigen_transformation.computeInverseAndDetWithCheck(eigen_inverse_transformation, determinant, invertible);
/*
    if(invertible)
      std::cout << "It is invertible, and its inverse is:" << std::endl << eigen_inverse_transformation << std::endl;
    else
      std::cout << "Not invertible and determinant was : " << determinant << std::endl;
*/
    full_matrix = eigen_transformation.transpose()*eigen_transformation;
  }

  Metric_base(const Vector_3& axis_x,  //normal (unit)
              const Vector_3& axis_y,  //principal direction 1 (unit) // axis_vpmax
              const Vector_3& axis_z,  //principal direction 2 (unit) // axis_vpmin
              const double& vpn,
              const double& vpy,
              const double& vpz,
              const double& epsilon)
  {
#ifdef ANISO_DEBUG_METRIC
    double xz = axis_x*axis_z;
    double xy = axis_x*axis_y;
    double yz = axis_y*axis_z;
    double xx = axis_x*axis_x;
    double yy = axis_y*axis_y;
    double zz = axis_z*axis_z;
    std::cout << "produits scalaires : "
              << xz << " " << xy << " " << yz << std::endl;
    std::cout << "sq_norms           : "
              << xx << " " << yy << " " << zz << std::endl;
    std::cout << "valeurs propres    : " << vpn << " " << vpy << " " << vpz << std::endl;
    std::cout << std::endl;
    if(xz > 0.01 || xy > 0.01 || yz > 0.01)
      std::cout << "Warning : This dot product should be 0" << std::endl;
    if(std::abs(xx-1.) > 0.01 || std::abs(yy-1.) > 0.01 || std::abs(zz-1.) > 0.01)
      std::cout << "Warning : This sqnorm should be 1" << std::endl;
#endif
    if(std::abs(vpy) > std::abs(vpz))
      construct(axis_x, axis_y, axis_z, vpn, vpy, vpz, epsilon);
    else
      construct(axis_x, axis_z, axis_y, vpn, vpz, vpy, epsilon);
  }

  ~Metric_base() { }
};

} // namespace Canvas
} // namespace CGAL

#endif // CGAL_CANVAS_METRIC_H
