// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Möller

#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_2<Kernel> Triangulation;

template<typename T>
void concept_check_triangulation() {
  boost::function_requires< boost::VertexListGraphConcept<T> >();
  boost::function_requires< boost::BidirectionalGraphConcept<T> >();
}

int main()
{
  concept_check_triangulation<Triangulation>();
  return 0;
}
