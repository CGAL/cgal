// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Steve OUDOT


#ifndef COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
#define COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/Triangulation_ds_cell_base_3.h>

namespace CGAL {
  namespace Mesh_3 {
  
  template < class GT, class Cb = Triangulation_ds_cell_base_3 <> > 
  class Complex_2_in_triangulation_cell_base_3 : public Cb {    
    
  public:
    typedef Complex_2_in_triangulation_cell_base_3 <GT, Cb> Self;
    
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS3>::Other  Cb3;
      typedef Complex_2_in_triangulation_cell_base_3 <GT, Cb3> Other;
    };
    
    
    typedef typename GT::Point_3 Point;
    
    typedef typename Cb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    
    
  private:
    bool domain;
    
  public:
    
    // Constructors

    Complex_2_in_triangulation_cell_base_3() : Cb(), domain(false)
    {
    }
    
    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3)
      :  Cb (v0, v1, v2, v3), domain(false)
    {
    }
    
    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3,
					    Cell_handle n0,
					    Cell_handle n1,
					    Cell_handle n2,
					    Cell_handle n3)
      : Cb (v0, v1, v2, v3, n0, n1, n2, n3), domain(false) 
    {
    }
    
    bool is_in_domain() const
    {
      return domain;
    }

    void set_in_domain(bool b)
    {
      domain = b;
    }
  };  // end Complex_2_in_triangulation_cell_base_3

  }//namespace Mesh_3
}  // namespace CGAL


#endif  // COMPLEX_2_IN_TRIANGULATION_3_H
