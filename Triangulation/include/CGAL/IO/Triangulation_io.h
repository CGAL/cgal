// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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
// $URL:  $
// $Id:  $
//
// Author(s)     : Clement Jamin


#ifndef CGAL_TRIANGULATION_IO_H
#define CGAL_TRIANGULATION_IO_H

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation_vertex.h>
#include <string>
#include <iostream>

namespace CGAL {

template<typename K>
std::ostream &
operator<<(std::ostream & os, const typename Wrap::Point_d<K> & p)
{
  typename K::Cartesian_const_iterator_d it = p.cartesian_begin();
  os << "(" << *it;
  ++it;
  for ( ; it != p.cartesian_end() ; ++it)
  {
    os << ", " << *it;
  }
  os << ")";
  return os;
}

template<typename K>
std::ostream &
operator<<(std::ostream & os, const typename Wrap::Weighted_point_d<K> & p)
{
  return os << p.point();
}

/*template< class A, class B >
std::ostream &
operator<<(std::ostream & os, const Triangulation_vertex<A, Data, B> & v)
{
  os << v.point();
  return os;
}*/

} //namespace CGAL

#endif // CGAL_TRIANGULATION_IO_H
