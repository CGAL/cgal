// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_FACET_ON_SURFACE_CRITERION_H
#define CGAL_FACET_ON_SURFACE_CRITERION_H

#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>

namespace CGAL {

  namespace Mesh_3 {

  template <class Tr>
  class Facet_on_surface_criterion : 
      public Surface_mesher::Refine_criterion <Tr> {
  public:
    typedef typename Surface_mesher::Refine_criterion <Tr> Criterion;
    typedef typename Criterion::Quality Quality;

  private:
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;

  public:
    bool is_bad (const Facet& f) const {
      Quality q = quality(f);
      if (q==Quality(0))
	{
	  std::cerr << "!";
	  return true;
	}
      else
	{
	  std::cerr << "?";
	  return false;
	}
    }
    
    Quality quality (const Facet& f) const {
      const Cell_handle& ch = f.first;
      const int i = f.second;
      const Vertex_handle& v1 = ch->vertex((i+1)&3);
      const Vertex_handle& v2 = ch->vertex((i+2)&3);
      const Vertex_handle& v3 = ch->vertex((i+3)&3);
      // On regarde alors si la facette a tous ces sommets sur la surface.
      if ( (!v1->info()) || (!v2->info()) || (!v3->info()) )
	return Quality(0);
      else return Quality(1);
    }
  }; // end class Facet_on_surface_criterion

 }  // end namespace Mesh_3

}  // end namespace CGAL
  
#endif // CGAL_FACET_ON_SURFACE_CRITERION_H
