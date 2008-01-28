// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_TRIANGLE_WITH_SURFACE_INDEX_H
#define CGAL_TRIANGLE_WITH_SURFACE_INDEX_H

#include <string>

namespace CGAL {

template <class Triangle>
class Triangle_with_surface_index : public Triangle
{
public:
  Triangle_with_surface_index() : Triangle(), index(0) {}

  Triangle_with_surface_index(const Triangle& t) : Triangle(t), index(0) {}

  Triangle_with_surface_index(const Triangle_with_surface_index& ti)
    : Triangle(ti), index(ti.surface_index()) {}

  int surface_index() const
  {
    return index;
  }

  void set_surface_index(const int i)
  {
    index = i;
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Triangle>()() + "+i";
  }
#endif
private:
  int index;
}; // end class Triangle_with_surface_index

template <class Triangle>
std::ostream&
operator<<(std::ostream &os, const Triangle_with_surface_index<Triangle>& t)
{
  os << static_cast<const Triangle&>(t);
  if(is_ascii(os))
    os << ' ' << t.surface_index();
  else
    write(os, t.surface_index());
  return os;
}

template <class Triangle>
std::istream&
operator>>(std::istream &is, Triangle_with_surface_index<Triangle>& t)
{
  is >>  static_cast<Triangle&>(t);
  int index;
  if(is_ascii(is))
    is >> index;
  else
    read(is, index);
  t.set_surface_index(index);
  return is;
}

} // end namespace CGAL

#endif
