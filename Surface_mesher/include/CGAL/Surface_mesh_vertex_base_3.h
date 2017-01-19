// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_VERTEX_BASE_3_H
#define CGAL_SURFACE_MESH_VERTEX_BASE_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>

namespace CGAL {

  template < class GT, class Vb = Triangulation_vertex_base_3 <GT> > 
  class Surface_mesh_vertex_base_3 
    : public Complex_2_in_triangulation_vertex_base_3<GT, Vb> {    
    
  public:
    typedef Surface_mesh_vertex_base_3 <GT, Vb> Self;
    
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS3>::Other  Vb3;
      typedef Surface_mesh_vertex_base_3 <GT, Vb3> Other;
    };
    
  public:  
    Surface_mesh_vertex_base_3()
      : Complex_2_in_triangulation_vertex_base_3<GT, Vb>()
    {}
  };  // end Surface_mesh_vertex_base_3

}  // namespace CGAL

#endif  // CGAL_SURFACE_MESH_CELL_BASE_3_H
