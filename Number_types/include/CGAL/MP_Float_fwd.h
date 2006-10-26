// Copyright (c) 2001-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_MP_FLOAT_FWD_H
#define CGAL_MP_FLOAT_FWD_H

#include <CGAL/Quotient_fwd.h>
#include <CGAL/Root_of_2_fwd.h>

// Forward declarations

CGAL_BEGIN_NAMESPACE

class MP_Float;

namespace INTERN_MP_FLOAT {
Comparison_result compare(const MP_Float&, const MP_Float&);
MP_Float square(const MP_Float&);
double to_double(const MP_Float&);
double to_double(const Quotient<MP_Float>&);
double to_double(const Root_of_2<MP_Float> &x);
std::pair<double,double> to_interval(const MP_Float &);
std::pair<double,double> to_interval(const Quotient<MP_Float>&);
MP_Float div(const MP_Float& n1, const MP_Float& n2);
MP_Float gcd( const MP_Float& a, const MP_Float& b);
}

MP_Float approximate_sqrt(const MP_Float&);
MP_Float exact_division(const MP_Float & n, const MP_Float & d);

CGAL_END_NAMESPACE

#endif // CGAL_MP_FLOAT_FWD_H
