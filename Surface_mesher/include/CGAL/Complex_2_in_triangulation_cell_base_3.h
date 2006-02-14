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
// $Id$ $Date$
// 
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri



#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/Triangulation_ds_cell_base_3.h>

namespace CGAL {

  
  template < class GT, class Cb = Triangulation_ds_cell_base_3 <> > 
  class Complex_2_in_triangulation_cell_base_3 : public Cb {    
    
  public:
    typedef Complex_2_in_triangulation_cell_base_3 <GT, Cb> Self;
    
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS3>::Other  Cb3;
      typedef Complex_2_in_triangulation_cell_base_3 <GT, Cb3> Other;
    };
    
    
    typedef typename Cb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    
    
  private:
    
    // Champs ajoutes a la classe
    
    // Facets
    char c_surface_facets;
    
  public:
    
    // Constructors

    Complex_2_in_triangulation_cell_base_3() 
      : Cb(), c_surface_facets(0) 
    {}
    
    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3) 
      : Cb (v0, v1, v2, v3), c_surface_facets(0) 
    {}
    
    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3,
					    Cell_handle n0,
					    Cell_handle n1,
					    Cell_handle n2,
					    Cell_handle n3) 
      : Cb (v0, v1, v2, v3, n0, n1, n2, n3), c_surface_facets(0) 
    {}
    
    
    
    // Access functions
    
    // Facets
  bool is_facet_on_surface(const int facet) const {
    CGAL_assertion (facet>=0 && facet <4);
    return c_surface_facets & (1 << facet);
  }
  
  // Setting functions
  
  // Facets
  void set_facet_on_surface(const int facet, const bool f) {
    CGAL_assertion (facet>=0 && facet <4);
    if(f){
      c_surface_facets |= (1 << facet);
    }else {
      c_surface_facets &= (15  & ~(1 << facet));
    }
  }
  };  // end Complex_2_in_triangulation_cell_base_3

template < class GT, class Cb >
inline
std::istream&
operator>>(std::istream &is, Complex_2_in_triangulation_cell_base_3<GT, Cb> &c)
{
  bool b;
  is >> static_cast<Cb&>(c);
  if( ! is )
      std:: cerr << "Complex_2_in_triangulation_cell_base_3::operator>>\n"
                 << "bad avant!\n";
  if( is.eof() )
    std::cerr << "Complex_2_in_triangulation_cell_base_3::operator>>\n"
              << "eof!\n";
 for(int i = 0; i < 4; ++i)
  {
    is >> b;
    if( ! is )
      std:: cerr << "Complex_2_in_triangulation_cell_base_3::operator>>\n"
                 << i << " bad\n";
    c.set_facet_on_surface(i, b);
  }
  return is;
}

template < class GT, class Cb >
inline
std::ostream&
operator<<(std::ostream &os,
           const Complex_2_in_triangulation_cell_base_3<GT, Cb> &c)
{
  os << static_cast<const Cb&>(c);
  for(int i = 0; i < 4; ++i)
  {
    if(is_ascii(os))
      os << ' ';
    os << c.is_facet_on_surface(i);
  }
  return os;
}

}  // namespace CGAL



#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H

