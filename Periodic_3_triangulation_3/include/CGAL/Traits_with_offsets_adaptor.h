// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H
#define CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H

namespace CGAL {

template < class K, class Functor_ >
  class Traits_with_offsets_adaptor {
  typedef K Kernel;
  typedef Functor_ Functor;

  typedef typename Kernel::Point_3       Point;
  typedef typename Kernel::Offset        Offset;

public:
  typedef typename Kernel::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename Kernel::Construct_point_3 Construct_point_3;
  typedef typename Functor::result_type result_type;

  Traits_with_offsets_adaptor(const Iso_cuboid_3 * dom) : _domain(dom) { }

  result_type operator()(const Point& p0, const Point& p1,
      const Offset& o0, const Offset& o1) const {
    return Functor()(pp(p0,o0),pp(p1,o1));
  }
  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
      const Offset& o0, const Offset& o1, const Offset& o2) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2));
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3,
      const Offset& o0, const Offset& o1,
      const Offset& o2, const Offset& o3) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2),pp(p3,o3));
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4,
      const Offset& o0, const Offset& o1, const Offset& o2,
      const Offset& o3, const Offset& o4) const {
    return Functor()(pp(p0,o0),pp(p1,o1),pp(p2,o2),
	pp(p3,o3),pp(p4,o4));
  }

  result_type operator()(const Point& p0, const Point& p1) const {
    return Functor()(p0, p1);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2) const {
    return Functor()(p0, p1, p2);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3) const {
    return Functor()(p0, p1, p2, p3);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4) const {
    return Functor()(p0, p1, p2, p3, p4);
  }

protected:
  Point pp(const Point &p, const Offset &o) const {
    return Construct_point_3(*_domain)(p,o);
  }
 public:
  const Iso_cuboid_3* _domain;
};
}  // namespace CGAL

#endif /* CGAL_TRAITS_WITH_OFFSETS_ADAPTOR_H */
