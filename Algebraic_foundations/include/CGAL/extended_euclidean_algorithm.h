// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Michael Hemmer, Dominik Huelse
//
// ============================================================================

#ifndef CGAL_EXTENDED_EUCLIDEAN_ALGORITHM_H
#define CGAL_EXTENDED_EUCLIDEAN_ALGORITHM_H 1

#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

// EEA computing the normalized gcd
// Modern Computer Algebra (Hardcover)
// by Joachim von zur Gathen (Author), Juergen Gerhard (Author)
// Publisher: Cambridge University Press; 2 edition (September 1, 2003)
//  Language: English
// ISBN-10: 0521826462
// ISBN-13: 978-0521826464
// pp.: 55

template< class AS >
AS extended_euclidean_algorithm(const AS& f, const AS& g, AS& s_, AS& t_){
    typename Algebraic_structure_traits<AS>::Integral_division idiv;
    typename Algebraic_structure_traits<AS>::Div div;
    typename Algebraic_structure_traits<AS>::Unit_part unit_part;

    std::vector<AS> p,r,s,t,q;
    p.push_back(unit_part(f));
    r.push_back(idiv(f,p[0]));
    s.push_back(idiv(AS(1),p[0]));
    t.push_back(AS(0));
    q.push_back(AS(0));

    p.push_back(unit_part(g));
    r.push_back(idiv(g,p[1]));
    s.push_back(AS(0));
    t.push_back(idiv(AS(1),p[1]));

    int i = 1;
    while(!is_zero(r[i])){
        q.push_back(div(r[i-1],r[i]));
        r.push_back(r[i-1]-q[i]*r[i]);
        p.push_back(unit_part(r[i+1]));
        r[i+1] = idiv(r[i+1],p[i+1]);
        s.push_back(idiv(s[i-1]-q[i]*s[i],p[i+1]));
        t.push_back(idiv(t[i-1]-q[i]*t[i],p[i+1]));
        i++;
    }

    s_=s[i-1];
    t_=t[i-1];
    AS h = r[i-1];
    CGAL_precondition( h == f*s_ + g*t_);
    return h;
}

} //namespace CGAL

#endif // NiX_EXTENDED_EUCLIDEAN_ALGORITHM_H //
