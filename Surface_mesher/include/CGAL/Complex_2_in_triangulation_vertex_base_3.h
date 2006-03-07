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
// $URL$
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
    bool in_complex;
    bool regular;
      
    bool status_cached;

  public:
    // Constructors

    Complex_2_in_triangulation_vertex_base_3()
      : Vb(), in_complex(false), regular(false), status_cached(false)
    {}

    bool is_in_complex() const {
      return in_complex;
    }

    void set_in_complex(const bool b) {
      in_complex = b;
    }

    bool is_regular() const {
      return regular;
    }

    void set_regular(const bool b) {
      regular = b;
    }

    bool is_status_cached() const {
      return status_cached;
    }

    void set_status_cached(const bool b) {
      status_cached = b;
    }
  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
