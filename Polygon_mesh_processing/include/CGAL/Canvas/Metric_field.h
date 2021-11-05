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

#ifndef CGAL_CANVAS_METRIC_FIELD_H
#define CGAL_CANVAS_METRIC_FIELD_H

#include <CGAL/Canvas/Metric.h>
#include <CGAL/Canvas/metric_helper.h>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL {
namespace Canvas {

template<typename K, typename KExact = K>// = CGAL::Exact_predicates_exact_constructions_kernel>
class Metric_field
{
public:
  typedef Metric_base<K, KExact>  Metric;
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_3;
  typedef typename K::Vector_3    Vector_3;

public:
  FT epsilon;
  double en_factor;

  virtual Metric compute_metric(const Point_3 &p) const = 0;

  Metric build_metric(const Vector_3& v0, const Vector_3& v1, const Vector_3& v2,
                      const double& e0, const double& e1, const double& e2) const
  {
    return Metric(v0, v1, v2, en_factor*e0, e1, e2, epsilon);
  }

  Metric intersection(Metric metric_p, Metric metric_q) const
  {
    Eigen::Matrix3d eigen_transformation_p = metric_p.get_transformation();
    Eigen::Matrix3d eigen_transformation_q = metric_q.get_transformation();
    Eigen::Matrix3d M_p = eigen_transformation_p.transpose()*eigen_transformation_p;
    Eigen::Matrix3d M_q = eigen_transformation_q.transpose()*eigen_transformation_q;

#ifdef ANISO_DEBUG_MATRIX_OPERATIONS
    std::cout << "now intersecting..." << std::endl;
    std::cout << "p" << std::endl << eigen_transformation_p << std::endl;
    std::cout << "p evalues : " << metric_p.get_third_eigenvalue() << " " << metric_p.get_max_eigenvalue() << " " << metric_p.get_min_eigenvalue() << std::endl;
    std::cout << "q" << std::endl << eigen_transformation_q << std::endl;
    std::cout << "q evalues : " << metric_q.get_third_eigenvalue() << " " << metric_q.get_max_eigenvalue() << " " << metric_q.get_min_eigenvalue() << std::endl;
    std::cout << "M_p" << std::endl << M_p << std::endl;
    std::cout << "M_q" << std::endl << M_q << std::endl;
#endif
    Eigen::Matrix3d intersected_metric = matrix_intersection<K>(M_p, M_q);

    Vector_3 intersected_v0, intersected_v1, intersected_v2;
    double e0, e1, e2;
    get_eigen_vecs_and_vals<K>(intersected_metric,
                               intersected_v0, intersected_v1, intersected_v2,
                               e0, e1, e2);

    //need F, not M
    e0 = std::sqrt(e0);
    e1 = std::sqrt(e1);
    e2 = std::sqrt(e2);

    Metric ret_met = build_metric(intersected_v0, intersected_v1, intersected_v2, e0, e1, e2);
    std::cout << "final intersected metric : " << std::endl << ret_met.get_transformation() << std::endl;

    return ret_met;
  }

  Metric_field(FT epsilon_ = 1.0,
               FT en_factor_ = 0.999) // to avoid the degen case of two equal evs
    :
      epsilon(epsilon_),
      en_factor(en_factor_)
  { }

  virtual ~Metric_field( ) { }
};

} // namespace Canvas
} // namespace CGAL

#endif // CGAL_CANVAS_METRIC_FIELD_H
