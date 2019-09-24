// Copyright (c) 2019 GeometryFactory
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERNAL_IO_GENERIC_FACEGRAPH_BUILDER_H
#define CGAL_INTERNAL_IO_GENERIC_FACEGRAPH_BUILDER_H

#include <deque>
#include <iostream>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <iostream>

namespace CGAL{

namespace internal {
namespace IO{
template <class Facegraph, class P, class Derived>
class Generic_facegraph_builder
{
protected:
  typedef P Point_3;
  typedef std::deque<Point_3> Points_3;
  typedef std::vector<std::size_t> Facet;
  typedef std::deque<Facet> Surface;

public:
  std::string name, color;
  Generic_facegraph_builder(std::istream& is_)
    : is(is_), counter(0)
  {}

  void do_construct(){} //specific to Facegraph (declared in BGL)
  void
  read(std::istream&, Points_3&, Surface&) {}

  void operator()( Facegraph& graph)
  {
    static_cast<Derived*>(this)->read(this->is, this->meshPoints, this->mesh);
    static_cast<Derived*>(this)->do_construct(graph);
  }

protected:
  std::istream& is;
  int counter;
  Points_3 meshPoints;
  Surface mesh;
};


} //end IO
}//end internal
}//end CGAL
#endif // CGAL_INTERNAL_IO_GENERIC_FACEGRAPH_BUILDER_H
