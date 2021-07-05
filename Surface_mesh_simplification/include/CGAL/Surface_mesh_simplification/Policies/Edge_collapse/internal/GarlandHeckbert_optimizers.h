// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_OPTIMIZERS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_OPTIMIZERS_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <Eigen/Dense>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {
template<typename GeomTraits>
class GarlandHeckbert_invertible_optimizer
{
  typedef typename GeomTraits::FT FT;
  typedef typename Eigen::Matrix<FT, 4, 1> Col_4;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;

  public:

  Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const
  {
    Mat_4 X;
    X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

    Col_4 opt_pt;

    opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)

    return opt_pt;
  }
};

template<typename GeomTraits>
class GarlandHeckbert_singular_optimizer
{
  typedef typename GeomTraits::FT FT;
  typedef typename Eigen::Matrix<FT, 4, 1> Col_4;
  typedef typename Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;

  public:
  Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const 
  {
    Mat_4 X;
    X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

    Col_4 opt_pt;

    if(X.determinant() == 0)
    {
      // not invertible
      const Col_4 p1mp0 = std::move(p1 - p0);
      const FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0, 0);
      const FT b = 2 * (p0.transpose() * aQuadric * p1mp0)(0, 0);

      if(a == 0)
      {
        if(b < 0)
          opt_pt = p1;
        else if(b == 0)
          opt_pt = 0.5 * (p0 + p1);
        else
          opt_pt = p0;
      }
      else
      {
        FT ext_t = -b/(2*a);
        if(ext_t < 0 || ext_t > 1 || a < 0)
        {
          // one of endpoints
          FT p0_cost = (p0.transpose() * aQuadric * p0)(0, 0);
          FT p1_cost = (p1.transpose() * aQuadric * p1)(0, 0);
          if(p0_cost > p1_cost)
            opt_pt = p1;
          else
            opt_pt = p0;
        }
        else
        {
          // extremum of the parabola
          opt_pt = p0 + ext_t * (p1 - p0);
        }
      }
    }
    else // invertible
    {
      opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
    }
    return opt_pt;
  }
};

} //namespace internal
} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_OPTIMIZERS_H
