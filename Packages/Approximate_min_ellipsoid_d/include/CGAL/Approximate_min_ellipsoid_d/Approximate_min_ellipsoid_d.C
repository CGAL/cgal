// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

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
  write_eps(const std::string& name) const
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

    // output the inscribed ellipse:
    const double alpha = 1/((1+achieved_epsilon())*(d+1));
    epsf.set_stroke_mode(Impl::Eps_export_2::Dashed);
    epsf.write_ellipse(to_double(defining_matrix(0,0)*alpha),
		       to_double(defining_matrix(1,1)*alpha),
		       to_double(defining_matrix(0,1)*alpha),
		       to_double(defining_vector(0)*alpha),
		       to_double(defining_vector(1)*alpha),
		       to_double(defining_scalar()*alpha-1.0));
  }

}
