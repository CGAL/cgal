// Copyright (c) 2005-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008       GeometryFactory (France)
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_VERTICES_ON_THE_SAME_PSC_ELEMENT_CRITERION_H
#define CGAL_VERTICES_ON_THE_SAME_PSC_ELEMENT_CRITERION_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Surface_mesher/Standard_criteria.h>

namespace CGAL {

namespace Surface_mesher {

template <typename Tr, typename Surface>
class Vertices_on_the_same_psc_element_criterion : 
    public Refine_criterion <Tr> {
public:
  typedef Refine_criterion <Tr> Criterion;
  typedef typename Criterion::Quality Quality;
  
private:
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;

  const Surface& surface;
  
public:
  Vertices_on_the_same_psc_element_criterion(const Surface& surface)
    : surface(surface)
  {
  }

  bool is_bad (const Facet& f, Quality& q) const {
    const Cell_handle& ch = f.first;
    const int i = f.second;
    const Vertex_handle& v1 = ch->vertex((i+1)&3);
    const Vertex_handle& v2 = ch->vertex((i+2)&3);
    const Vertex_handle& v3 = ch->vertex((i+3)&3);

    const bool is_bad = 
      surface.vertices_not_on_same_surface_patch(v1, v2, v3);

    q = (is_bad ? Quality(0) : Quality(1));
    if(is_bad){
#ifdef CGAL_SURFACE_MESHER_VERBOSE
      CGAL_MESHES_OUTPUT_STREAM << "f("
                                << v1->point().element_index() << ","
                                << v2->point().element_index() << ","
                                << v3->point().element_index() << ")";
#endif
    }
    return is_bad;
  }
}; // end Vertices_on_the_same_psc_element_criterion


} // end namespace Surface_mesher

} // end namespace CGAL


#endif // CGAL_VERTICES_ON_THE_SAME_PSC_ELEMENT_CRITERION_H
