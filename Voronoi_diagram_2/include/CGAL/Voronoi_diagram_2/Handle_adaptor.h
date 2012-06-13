// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_HANDLE_ADAPTOR_H
#define CGAL_VORONOI_DIAGRAM_2_HANDLE_ADAPTOR_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

template<class T>
class Handle_adaptor
{
 private:
  typedef Handle_adaptor<T>  Self;
 public:
  typedef T      value_type;
  typedef T*     pointer;
  typedef T&     reference;
  typedef const T*  const_pointer;
  typedef const T&  const_reference;

 public:
  Handle_adaptor() : t() {}
  Handle_adaptor(const T& t) : t(t) {}

  pointer    operator->() { return &t; }
  reference  operator*() { return t; }

  const_pointer    operator->() const { return &t; }
  const_reference  operator*() const { return t; }

  bool operator==(const Self& other) const {
    return t == other.t;
  }

  bool operator!=(const Self& other) const {
    return t != other.t;
  }

  bool operator<(const Self& other) const {
    return t < other.t;
  }

 private:
  T t;
};

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_HANDLE_ADAPTOR_H
