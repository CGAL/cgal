// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_CL_RATIONAL_H
#define CGAL_CL_RATIONAL_H

#include <CGAL/CLN/common.h>
#include <cl_rational.h>
#include <cl_rational_io.h>

CGAL_BEGIN_NAMESPACE

// Requirements.

inline double	to_double	(const cl_RA &I) { return cl_double_approx(I); }

// We hope the artificially added error is enough...
inline Interval_base to_interval (const cl_RA & z)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  Interval_nt_advanced approx(cl_double_approx(z));
  FPU_set_cw(CGAL_FE_UPWARD);

  return ( (approx + Interval_base::Smallest) + Interval_base::Smallest)
         + Interval_base::Smallest;
}

// Specialized utilities.

inline bool is_negative		(const cl_RA &I) { return minusp(I); } 
inline bool is_positive		(const cl_RA &I) { return plusp(I); }
inline bool is_zero		(const cl_RA &I) { return zerop(I); }
inline Comparison_result compare (const cl_RA &I, const cl_RA &J)
{ return Comparison_result(cl_compare(I,J)); }

CGAL_END_NAMESPACE

#endif // CGAL_CL_RATIONAL_H
