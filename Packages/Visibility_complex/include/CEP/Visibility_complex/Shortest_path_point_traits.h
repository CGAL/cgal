// Copyright (c) 2001-2004  ENS of Paris (France).
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef SHORTEST_PATH_POINT_TRAITS_H
#define SHORTEST_PATH_POINT_TRAITS_H

#ifndef VISIBILITY_COMPLEX_POINT_TRAITS_H
#include <CEP/Visibility_complex/Visibility_complex_point_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------- 

template < class _R , class _E >
class Shortest_path_point_traits 
    : public Visibility_complex_point_traits<_R>
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_point_traits<_R>   Base;
    typedef typename Base::Disk                   Disk;
    typedef typename Base::Bitangent_2            Bitangent_2;
    typedef typename Base::Arc_2                  Arc_2;
    typedef typename Base::Point_2                Point_2;
    typedef _E                                    Exact_NT;       
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    typedef Simple_cartesian<Exact_NT>::Point_2   Exact_point_2;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    Disk make_convex_from_point(const Point_2& p) 
    { return p; }
    bool is_vertex(const Disk& c) { return true; }

    Exact_NT length(const Arc_2& /*a*/, 
		    const Bitangent_2& /*inf*/, const Bitangent_2& /*sup*/) const 
    { return 0; }
    Exact_NT length(const Bitangent_2& b) const 
	{ return distance(b.source(),b.target()); }
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    Exact_point_2 make_exact(const Point_2& p) const
    {
	Exact_NT px = static_cast<Exact_NT>(CGAL::to_double(p.x()));
	Exact_NT py = static_cast<Exact_NT>(CGAL::to_double(p.y()));
	return Exact_point_2(px,py);
    }
    Exact_NT distance(const Point_2& p, const Point_2& q) const 
    {
	return CGAL_NTS sqrt(squared_distance(make_exact(p),make_exact(q)));
    }
    // -------------------------------------------------------------------------
};

// ----------------------------------------------------------------------------- 

CGAL_END_NAMESPACE

#endif // SHORTEST_PATH_TRAITS_H
