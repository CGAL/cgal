// Copyright (c) 2005  Stanford University (USA).
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
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_CARTESIAN_MOVING_ORTHOSPHERE_H
#define CGAL_KINETIC_CARTESIAN_MOVING_ORTHOSPHERE_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/determinant.h>

namespace CGAL { namespace Kinetic { namespace internal {

template <class K>
class Cartesian_moving_orthosphere_3
{
    typedef typename K::Motion_function MF;
    typedef typename K::Point_3 PT;
    public:
        template <class P>
            Cartesian_moving_orthosphere_3(const P &a,
            const P &b,
            const P &c,
        const P &d, const K &k) {
            typename K::Delaunay_lifting lift= k.Delaunay_lifting_3_object();
            typename K::Center center k.center_3_object();
            initialize(center(a), lift(a),
                center(b), lift(b),
                center(c), lift(c),
                center(d), lift(d));
        }
        template <class P>
            Cartesian_moving_orthosphere_3(const P &a,
            const P &b,
        const P &c, const K &k) {
            typename K::Delaunay_lifting lift= k.Delaunay_lifting_3_object();
            typename K::Center center k.center_3_object();
            initialize(center(a), lift(a),
                center(b), lift(b),
                center(c), lift(c));
        }

        const PT& center() const
        {
            return center_;
        }
        const MF l() const
        {
            return l_;
        };
        const MF o() const
        {
            return o_;
        };
        const PT &origin()
        protected:
        PT center_;
        PT origin_;
        MF l_, o_;
};

} } } //namespace CGAL::Kinetic::internal
#endif
