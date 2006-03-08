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
    bool in_complex_mark_;
    bool regular_or_boundary_mark_;
    bool in_complex_validity_mark_;
    bool regular_or_boundary_validity_mark_;


  public:
    // Constructors

    Complex_2_in_triangulation_vertex_base_3()
      : Vb(), 
	in_complex_mark_(false), 
	regular_or_boundary_mark_(false), 
	in_complex_validity_mark_(false),
	regular_or_boundary_validity_mark_(false)
    {}

    bool in_complex_mark() const {
      return in_complex_mark_;
    }

    void set_in_complex_mark(const bool b) {
      in_complex_mark_ = b;
    }

    bool regular_or_boundary_mark() const {
      return regular_or_boundary_mark_;
    }

    void set_regular_or_boundary_mark(const bool b) {
      regular_or_boundary_mark_ = b;
    }

   
    bool in_complex_validity_mark() const {
      return in_complex_validity_mark_;
    }

    void set_in_complex_validity_mark(const bool b) {
      in_complex_validity_mark_ = b;
    }

    bool regular_or_boundary_mark() const {
      return regular_or_boundary_validity_mark_;
    }

    void set_regular_or_boundary_validity_mark(const bool b) {
      regular_or_boundary_validity_mark_ = b;
    }



  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
