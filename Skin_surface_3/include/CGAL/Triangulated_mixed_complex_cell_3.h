// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef TRIANGULATED_MIXED_COMPLEX_CELL_3_H
#define TRIANGULATED_MIXED_COMPLEX_CELL_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>
#include <CGAL/Skin_surface_quadratic_surface_3.h>

#include <CGAL/Triangulation_simplex_3.h>

CGAL_BEGIN_NAMESPACE

template < class GT, 
           class Polyhedron_3Kernel, 
	   class Cb = Triangulation_cell_base_3<GT> >
class Triangulated_mixed_complex_cell_3 : public Cb
{
public:
  typedef typename Cb::Triangulation_data_structure            Triangulation_data_structure;
  typedef typename Triangulation_data_structure::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure::Cell_handle   Cell_handle;

  typedef Polyhedron_3Kernel                                    Polyhedron_3_kernel;
  typedef Skin_surface_quadratic_surface_3<Polyhedron_3_kernel> QuadrSurface;
	
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Triangulated_mixed_complex_cell_3<GT, Polyhedron_3Kernel, Cb2>
                                                           Other;
  };

  Triangulated_mixed_complex_cell_3() : Cb() {
  }
  Triangulated_mixed_complex_cell_3(Vertex_handle v0, Vertex_handle v1,
				    Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {
  }

  QuadrSurface *surf;
};
	


CGAL_END_NAMESPACE

#endif // TRIANGULATED_MIXED_COMPLEX_CELL_3_H
