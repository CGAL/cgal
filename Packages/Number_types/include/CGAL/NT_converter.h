// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_NT_CONVERTER_H
#define CGAL_NT_CONVERTER_H

#include <CGAL/functional_base.h>

CGAL_BEGIN_NAMESPACE

// A number type converter usable as default, using the conversion operator.

template < class NT1, class NT2 >
struct NT_converter
  : public CGAL_STD::unary_function< NT1, NT2 >
{
    NT2
    operator()(const NT1 &a) const
    {
        return NT2(a);
    }
};

// ... and specialized for double to call to_double().

template < class NT1 >
struct NT_converter < NT1, double >
  : public CGAL_STD::unary_function< NT1, double >
{
    double
    operator()(const NT1 &a) const
    {
        return CGAL_NTS to_double(a);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_NT_CONVERTER_H
