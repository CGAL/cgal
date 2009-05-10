// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef TRIANGLE_3_LINE_3_INTERSECTION_H_
#define TRIANGLE_3_LINE_3_INTERSECTION_H_

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

namespace CGAL {
namespace CGALi {

template <class K>
Object
intersection(const typename K::Triangle_3 &t,
             const typename K::Line_3 &l,
             const K&)
{
  return CGAL::intersection(t.supporting_plane(),
                            l);
}

} // end namespace CGALi

template <class K>
inline
Object
intersection(const Triangle_3<K> &t, const Line_3<K> &l)
{
  return typename K::Intersect_3()(t,l);
}

template <class K>
inline
Object
intersection(const Line_3<K> &l, const Triangle_3<K> &t)
{
  return typename K::Intersect_3()(t,l);
}

} // end namespace CGAL

#endif // TRIANGLE_3_LINE_3_INTERSECTION_H_
