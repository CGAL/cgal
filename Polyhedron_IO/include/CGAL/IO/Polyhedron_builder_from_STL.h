// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_POLYHEDRON_STL_BUILDER_H
#define CGAL_IO_POLYHEDRON_STL_BUILDER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/STL_reader.h>

#include <iostream>

namespace CGAL{

template <class HDS>
class Polyhedron_builder_from_STL : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Vertex::Point Point_3;
  typedef std::vector<std::array<double, 3> > Points_3;
  typedef std::array<int,3> Facet;
  typedef std::vector<Facet> Surface;

  std::istream& is;
  Points_3 meshPoints;
  Surface mesh;

public:

  Polyhedron_builder_from_STL(std::istream& is_)
    : is(is_)
  {}

  void operator()( HDS& hds) {
    if(!read_STL(is, meshPoints, mesh)) return;

    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
    B.begin_surface( meshPoints.size(), mesh.size());
    typedef typename Points_3::size_type size_type;

    for(size_type i=0; i < meshPoints.size(); i++){
      B.add_vertex(
        Point_3(meshPoints[i][0], meshPoints[i][1], meshPoints[i][2])
      );
    }

    for(size_type i=0; i < mesh.size(); i++){
      B.begin_facet();
      B.add_vertex_to_facet( mesh[i][0]);
      B.add_vertex_to_facet( mesh[i][1]);
      B.add_vertex_to_facet( mesh[i][2]);
      B.end_facet();
    }
    if(B.error())
      {
        std::cerr << "An error occurred while creating a Polyhedron" << std::endl;
        B.rollback();
      }

    B.end_surface();
  }
};

} //end of CGAL namespace

#endif // CGAL_IO_POLYHEDRON_STL_BUILDER_H
