// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_SEGMENT_WITH_SURFACE_INDEX_H
#define CGAL_SEGMENT_WITH_SURFACE_INDEX_H

#include <CGAL/Point_traits.h>

#include <string>

namespace CGAL {

template <class Segment>
class Segment_with_surface_index : public Segment
{
public:
  Segment_with_surface_index() : Segment(), index(0) {}

  Segment_with_surface_index(const Segment& s) : Segment(s), index(0) {}

  Segment_with_surface_index(const Segment_with_surface_index& si)
    : Segment(si), index(si.surface_index()) {}

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
    return Get_io_signature<Segment>()() + "+i";
  }
#endif
private:
  int index;
}; // end class Segment_with_surface_index

template <class Segment>
std::ostream&
operator<<(std::ostream &os, const Segment_with_surface_index<Segment>& s)
{
  os << static_cast<const Segment&>(s);
  if(is_ascii(os))
    os << ' ' << s.surface_index();
  else
    write(os, s.surface_index());
  return os;
}

template <class Segment>
std::istream&
operator>>(std::istream &is, Segment_with_surface_index<Segment>& s)
{
  is >>  static_cast<Segment&>(s);
  int index;
  if(is_ascii(is))
    is >> index;
  else
    read(is, index);
  s.set_surface_index(index);
  return is;
}

} // end namespace CGAL

#endif
