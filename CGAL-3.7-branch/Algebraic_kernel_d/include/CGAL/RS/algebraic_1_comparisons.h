// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#ifndef CGAL_RS_ALGEBRAIC_1_COMPARISONS_H
#define CGAL_RS_ALGEBRAIC_1_COMPARISONS_H

#include <boost/config.hpp>
#include <CGAL/RS/compare_1.h>
#include <CGAL/RS/polynomial_1_utils.h>

namespace CGAL{

inline
bool operator<(const Algebraic_1 &n1,const Algebraic_1 &n2){
        typedef CGAL::Rsgcd_1  Gcd;
        return(CGAL::RS_COMPARE::compare_1<Gcd>(n1,n2)==SMALLER);
}

inline
bool operator==(const Algebraic_1 &n1,const Algebraic_1 &n2){
        typedef CGAL::Rsgcd_1  Gcd;
        return(CGAL::RS_COMPARE::compare_1<Gcd>(n1,n2)==EQUAL);
}

inline
Algebraic_1 min BOOST_PREVENT_MACRO_SUBSTITUTION (const Algebraic_1 &a,const Algebraic_1 &b){
        return (a<b?a:b);
}

inline
Algebraic_1 max BOOST_PREVENT_MACRO_SUBSTITUTION (const Algebraic_1 &a,const Algebraic_1 &b){
        return (a>b?a:b);
}

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_COMPARISONS_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
