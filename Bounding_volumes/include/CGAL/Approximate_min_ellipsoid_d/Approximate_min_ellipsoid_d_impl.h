// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#include <CGAL/eigen_2.h>
#include <CGAL/eigen.h>

#include <CGAL/Approximate_min_ellipsoid_d.h>

namespace CGAL {

  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  find_lower_dimensional_approximation()
  {
    // (No implementation yet, in accordance with the current
    // specification in doc_tex ...)
  }

  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  compute_center()
  {
    // According to (****), the computed ellipsoid E* has the representation
    //
    //   E* = { y in R^d | y^T M'/alpha y + 2/alpha y^Tm + (nu/alpha-1) <= 0 }
    // 
    // for
    //
    //       [ M'  m  ]
    //   M = [ m^T nu ]
    //
    // where M is the matrix defined via E->matrix(i,j).  We need
    // something in the form
    //
    //    E* = { y | (y - c)^T M'/alpha (y - c) + mu <= 0 }.
    //
    // Expanding the later and comparing with the original form we
    // obtain
    //
    //    c = - M'^{-1} m
    //
    // as the formula for the ellipsoid's center.  Comparing
    // coefficients we also get
    //
    //    mu = nu/alpha-1 - c^T M'/alpha c 
    //       = nu/alpha-1 + c^T m / alpha
    //       = (nu + c^Tm)/alpha - 1                               (********)
    //
    // In order to compute c, we will compute the inverse of M' and
    // multiply with m.

    // precondition checking:
    CGAL_APPEL_ASSERT(!has_center);

    // compute M'^{-1}:
    mi.resize(d * d);
    E->compute_inverse_of_submatrix(mi.begin());

    // compute center:
    center_.resize(d);
    for (int i=0; i<d; ++i) {
      FT ci(0);
      for (int j=0; j<d; ++j)
	ci += mi[i+d*j] * E->matrix(d,j);
      center_[i] = -ci;
    }

    // remember that we have computed the center:
    has_center = true;
  }

  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  compute_axes_2_3()
  {
    // According to (****), the computed ellipsoid E* has the representation
    //
    //   E* = { y in R^d | y^T M'/alpha y + 2/alpha y^Tm + (nu/alpha-1) <= 0 }
    // 
    // for
    //
    //       [ M'  m  ]
    //   M = [ m^T nu ]
    //
    // where M is the matrix defined via E->matrix(i,j).  After caling
    // compute_center() (see above), we have in center_ a point c such
    // that
    // 
    //   E* = { y | (y - c)^T M'/alpha (y - c) + mu <= 0 }.
    //
    // where mu = nu/alpha-1 - c^T M'/alpha c.
    //
    // Now if we can write M' = U D U^T holds for some diagonal matrix
    // D and an orthogonal matrix U then the length l_i of the ith axes
    // (corresponding to the ith "direcion" stored in the ith row of
    // U) can be obtained by plugging (0,...,0,l_i,0,...,0)U^T=y-c into
    // the above equation for E*:
    //
    //   l_i^2 d[i]/alpha = -mu, 
    //
    // which gives l_i = sqrt(-mu*alpha/d[i]).

    // precondition checking:
    CGAL_APPEL_ASSERT(!has_axes && lengths_.size() == 0 &&
		      directions_.size() == 0);

    // compute M'^{-1}, if need be:
    if (!has_center)
      compute_center();

    const double alpha = (1+achieved_epsilon()) * (d+1);

    // compute mu according to (********):
    double mu = E->matrix(d,d);
    for (int i=0; i<d; ++i)
      mu += center_[i] * E->matrix(d,i);
    mu = mu/alpha - 1.0;
    const double factor = -mu*alpha;
    
    // compute Eigendecomposition:
    if (d == 2)
      compute_axes_2(alpha, factor);
    else
      compute_axes_3(alpha, factor);

    // remember that we have computed the axes:
    has_axes = true;
  }

  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  compute_axes_2(const double /* alpha */, const double factor)
  {
    CGAL_APPEL_ASSERT(d==2);

    typedef Simple_cartesian<double> K;
    typedef Vector_2<K> Vector_2;

    // write matrix M' as [ a, b; b, c ]:
    const double matrix[3] = { E->matrix(0, 0),   // a
			       E->matrix(0, 1),   // b
			       E->matrix(1, 1) }; // c
    
    std::pair<Vector_2, Vector_2> eigenvectors; // Note: not neces. normalized.
    std::pair<double, double>     eigenvalues;  // Note: sorted descendent.
    internal::eigen_symmetric_2<K>(matrix, eigenvectors, eigenvalues);
    
    // normalize eigenvectors:
    double l1=1.0/std::sqrt(eigenvectors.first.x()*eigenvectors.first.x()+
			    eigenvectors.first.y()*eigenvectors.first.y());
    double l2=1.0/std::sqrt(eigenvectors.second.x()*eigenvectors.second.x()+
			    eigenvectors.second.y()*eigenvectors.second.y());
    
    // store axes lengths:
    lengths_.push_back(std::sqrt(factor/eigenvalues.first));
    lengths_.push_back(std::sqrt(factor/eigenvalues.second));
    
    // store directions:
    directions_.resize(2);
    directions_[0].push_back(eigenvectors.first.x()*l1);
    directions_[0].push_back(eigenvectors.first.y()*l1);
    directions_[1].push_back(eigenvectors.second.x()*l2);
    directions_[1].push_back(eigenvectors.second.y()*l2);
  }
  
  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  compute_axes_3(const double /* alpha */, const double factor)
  {
    CGAL_APPEL_ASSERT(d==3);

    // write matrix M' as
    //
    //        [ a b c ]
    //   M' = [ b d e ]
    //        [ c e f ]
    //
    const double matrix[6] = { E->matrix(0, 0),   // a
			       E->matrix(0, 1),   // b
			       E->matrix(1, 1),   // d
			       E->matrix(0, 2),   // c
			       E->matrix(1, 2),   // e
			       E->matrix(2, 2) }; // f
    
    double eigenvectors[3 * 3]; // Note: not necessarily normalized.
    double eigenvalues[3];      // Note: sorted descendent.
    internal::eigen_symmetric<double>(matrix, 3, eigenvectors, eigenvalues);
    
    // normalize eigenvectors:
    double l1 = 1.0/std::sqrt(eigenvectors[0] * eigenvectors[0]+  // x^2
			      eigenvectors[1] * eigenvectors[1]+  // y^2
			      eigenvectors[2] * eigenvectors[2]); // z^2
    double l2 = 1.0/std::sqrt(eigenvectors[3] * eigenvectors[3]+  // x^2
			      eigenvectors[4] * eigenvectors[4]+  // y^2
			      eigenvectors[5] * eigenvectors[5]); // z^2
    double l3 = 1.0/std::sqrt(eigenvectors[6] * eigenvectors[6]+  // x^2
			      eigenvectors[7] * eigenvectors[7]+  // y^2
			      eigenvectors[8] * eigenvectors[8]); // z^2
    
    // store axes lengths:
    for (int i=0; i<3; ++i)
      lengths_.push_back(std::sqrt(factor/eigenvalues[i]));
    
    // store directions:
    directions_.resize(3);
    directions_[0].push_back(eigenvectors[0]*l1);
    directions_[0].push_back(eigenvectors[1]*l1);
    directions_[0].push_back(eigenvectors[2]*l1);
    directions_[1].push_back(eigenvectors[3]*l2);
    directions_[1].push_back(eigenvectors[4]*l2);
    directions_[1].push_back(eigenvectors[5]*l2);
    directions_[2].push_back(eigenvectors[6]*l3);
    directions_[2].push_back(eigenvectors[7]*l3);
    directions_[2].push_back(eigenvectors[8]*l3);
  }

  template<class Traits>
  void Approximate_min_ellipsoid_d<Traits>::
  write_eps(const std::string& name)
  {
    CGAL_APPEL_ASSERT(d==2 && !is_degenerate());

    namespace Impl = Approximate_min_ellipsoid_d_impl;
    Impl::Eps_export_2 epsf(name,2.0);
    epsf.set_label_mode(Impl::Eps_export_2::None);

    // output the input points:
    for (unsigned int k=0; k<P.size(); ++k) {
      C_it pk = tco.cartesian_begin(P[k]);
      const double u = to_double(*pk++);
      const double v = to_double(*pk);
      epsf.write_circle(u,v,0.0);
    }

    // output center:
    Center_coordinate_iterator c = center_cartesian_begin();
    const double cx = *c;
    const double cy = *++c;
    epsf.write_circle(cx, cy, 0.0, true); // (Note: last arg says to color.)

    // output axes:
    Axes_lengths_iterator axes = axes_lengths_begin();
    const double a1 = *axes;
    const double a2 = *++axes;
    Axes_direction_coordinate_iterator d1 = axis_direction_cartesian_begin(0);
    const double d1_x = *d1;
    const double d1_y = *++d1;
    Axes_direction_coordinate_iterator d2 = axis_direction_cartesian_begin(1);
    const double d2_x = *d2;
    const double d2_y = *++d2;
    epsf.write_line(cx, cy, cx+a1*d1_x, cy+a1*d1_y, true);
    epsf.write_line(cx, cy, cx+a2*d2_x, cy+a2*d2_y, true);

    // output the inscribed ellipse:
    const double alpha_inv = 1.0/((1+achieved_epsilon())*(d+1));
    epsf.set_stroke_mode(Impl::Eps_export_2::Dashed);
    epsf.write_ellipse(to_double(defining_matrix(0,0)*alpha_inv),
		       to_double(defining_matrix(1,1)*alpha_inv),
		       to_double(defining_matrix(0,1)*alpha_inv),
		       to_double(defining_vector(0)*alpha_inv),
		       to_double(defining_vector(1)*alpha_inv),
		       to_double(defining_scalar()*alpha_inv-1.0));
  }

}
