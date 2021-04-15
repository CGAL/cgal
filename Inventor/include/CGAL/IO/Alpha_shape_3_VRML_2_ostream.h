// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_ALPHA_SHAPE_3_VRML_2_OSTREAM_H
#define CGAL_ALPHA_SHAPE_3_VRML_2_OSTREAM_H

#include <CGAL/basic.h>
#include <CGAL/IO/VRML_2_ostream.h>

#ifdef CGAL_ALPHA_SHAPE_3_H
namespace CGAL {

template <class Dt >
VRML_2_ostream&
operator<<(VRML_2_ostream& os,
           Alpha_shape_3<Dt> &as)
{
  // Finite vertices coordinates.
  Alpha_shape_3<Dt>::Alpha_shape_vertices_iterator Vlist_it,
    Vlist_begin = as.alpha_shape_vertices_begin(),
    Vlist_end = as.alpha_shape_vertices_end();

  std::map<Alpha_shape_3<Dt>::Vertex_handle, int> V;
  int number_of_vertex = 0;
  for( Vlist_it = Vlist_begin; Vlist_it != Vlist_end; Vlist_it++) {
    V[*Vlist_it] = number_of_vertex++;
  }

  typename Alpha_shape_3<Dt>::Alpha_shape_facets_iterator Flist_it,
    Flist_begin = as.alpha_shape_facets_begin(),
    Flist_end = as.alpha_shape_facets_end();

  std::map<Alpha_shape_3<Dt>::Facet, int> F;
  int number_of_facets = 0;
  for( Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++) {
    F[*Flist_it] = number_of_facets++;
  }

  const char *Indent = "                                    ";
  os <<      "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance USE A1\n"
             "                    geometry\n"
             "                        IndexedFaceSet {\n"
             "                            coord Coordinate {\n"
             "                                point [ \n"
     <<      Indent << "  ";
  for( Vlist_it = Vlist_begin; Vlist_it != Vlist_end; Vlist_it++) {
    os << CGAL::to_double((*Vlist_it)->point().x()) << " ";
    os << CGAL::to_double((*Vlist_it)->point().y()) << " ";
    os << CGAL::to_double((*Vlist_it)->point().z()) << ",\n" << Indent << "  ";
  }
    os <<    "\n                                ]\n"
             "                            } # coord\n"
             "                            solid   FALSE\n"
     <<      Indent << "coordIndex  [\n";
  // Finite facets indices.
  for( Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++){
    os << Indent << "  ";
      for (int i=0; i<4; i++)
        if (i != (*Flist_it).second){
                os << V[(*Flist_it).first->vertex(i)];
          os << ", ";
        }
    if (Flist_it != Flist_end)
      os << "-1,\n";
    else
      os << "-1 \n";
  }
  os <<      Indent << "]\n";
             "                        } #IndexedFaceSet\n"
             "                } #Shape\n"
             "            ] #children\n"
             "        } #Group\n";

  return os;
}

} //namespace CGAL
#endif // CGAL_ALPHA_SHAPE_3_H

#endif CGAL_ALPHA_SHAPE_3_VRML_2_OSTREAM_H
