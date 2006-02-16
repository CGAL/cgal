// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
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
// $Id$
// 
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri




#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_vertex_base_3.h>

namespace CGAL {

  template < class GT, class Vb = Triangulation_vertex_base_3 <GT> > 
    class Complex_2_in_triangulation_vertex_base_3 : public Vb {    
    
  public:
    typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb> Self;
    
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS3>::Other  Vb3;
      typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb3> Other;
    };
    
    typedef typename Vb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    typedef typename Tds::Facet Facet;


  private:
    bool visited;

    public: // AF: todo: make private and wrap in functions
      bool regular_is_cached;
      bool regular;

  public:  
    // Constructors

    Complex_2_in_triangulation_vertex_base_3()
      : Vb(), visited(false), regular_is_cached(false), regular(false)
      {}


    bool is_visited() const {
      return visited;
    }
    
    void set_visited(const bool b) {
      visited = b;
    }

  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
