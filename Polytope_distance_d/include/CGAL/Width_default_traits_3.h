// Copyright (c) 1997-2000  ETH Zurich (Switzerland).
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
// Author(s)     : Thomas Herrmann, Lutz Kettner

#ifndef CGAL_WIDTH_DEFAULT_TRAITS_3_H
#define CGAL_WIDTH_DEFAULT_TRAITS_3_H

#include <CGAL/Convex_hull_traits_3.h>

namespace CGAL {

template <class Kernel_>
class Width_default_traits_3 {
public:
    typedef Kernel_                      Kernel;
    typedef typename Kernel::RT          RT;
    typedef typename Kernel::Point_3     Point_3;
    typedef typename Kernel::Plane_3     Plane_3;
    typedef typename Kernel::Vector_3    Vector_3;
    typedef Convex_hull_traits_3<Kernel> ChullTraits;

    RT get_hx( const Point_3& p) const { return p.hx(); }
    RT get_hy( const Point_3& p) const { return p.hy(); }
    RT get_hz( const Point_3& p) const { return p.hz(); }
    RT get_hw( const Point_3& p) const { return p.hw(); }
    void get_point_coordinates( const Point_3& p, 
                                RT& px, RT& py, RT& pz, RT& ph) const {
        px = get_hx(p);
        py = get_hy(p);
        pz = get_hz(p);
        ph = get_hw(p);
    }
    RT get_a( const Plane_3& f) const { return f.a(); }
    RT get_b( const Plane_3& f) const { return f.b(); }
    RT get_c( const Plane_3& f) const { return f.c(); }
    RT get_d( const Plane_3& f) const { return f.d(); }
    void get_plane_coefficients( const Plane_3& f, 
                                 RT& a, RT& b, RT& c, RT& d) const {
        a = get_a(f);
        b = get_b(f);
        c = get_c(f);
        d = get_d(f);
    }
    
    Point_3 make_point( const RT& hx, const RT& hy, 
                        const RT& hz, const RT& hw) const  {
        return Point_3(hx,hy,hz,hw);
    }
    Plane_3 make_plane( const RT& a, const RT& b, 
                        const RT& c, const RT& d) const {
        return Plane_3(a,b,c,d);
    }
    Vector_3 make_vector( const RT& a, const RT& b, const RT& c) const {
        return Vector_3(a,b,c);
    }
};

} //namespace CGAL

#endif //CGAL_WIDTH_DEFAULT_TRAITS_3_H
