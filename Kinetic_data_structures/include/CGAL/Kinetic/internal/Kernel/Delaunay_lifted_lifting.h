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

#ifndef CGAL_KINETIC_INTERNAL_DELAUNAY_LIFTED_LIFTING_H
#define CGAL_KINETIC_INTERNAL_DELAUNAY_LIFTED_LIFTING_H
#include <CGAL/Kinetic/basic.h>

namespace CGAL { namespace Kinetic { namespace internal {

template <class K>
struct Delaunay_lifted_lifting
{
    typedef typename K::Point_3 argument_type;
    typedef typename K::Motion_function result_type;

    result_type operator()(const argument_type &p) const
    {
        return CGAL::square(p.x())+ CGAL::square(p.y()) + CGAL::square(p.z());
    }
    result_type operator()(const typename K::Weighted_point_3 &wp) const
    {
        return wp.lifted();
    }
};

} } } //namespace CGAL::Kinetic::internal
#endif
