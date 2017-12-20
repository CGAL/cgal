// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri



#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/license/Surface_mesher.h>



#include <CGAL/Triangulation_cell_base_3.h>

#include <bitset>
#include <string>

namespace CGAL {


  template < class GT, class Cb = Triangulation_cell_base_3 <GT> >
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
    std::bitset<4> bits;
//     char c_surface_facets;

  public:

    // Constructors

    Complex_2_in_triangulation_cell_base_3()
      : Cb(), bits()// c_surface_facets(0)
    {}

    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3)
      : Cb (v0, v1, v2, v3), bits()// , c_surface_facets(0)
    {}

    Complex_2_in_triangulation_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3,
					    Cell_handle n0,
					    Cell_handle n1,
					    Cell_handle n2,
					    Cell_handle n3)
      : Cb (v0, v1, v2, v3, n0, n1, n2, n3), bits()// c_surface_facets(0)
    {}



    // Access functions

    // Facets
  bool is_facet_on_surface(const int facet) const {
    CGAL_assertion (facet>=0 && facet <4);
    return bits[facet];
//     return c_surface_facets & (1 << facet);
  }

  // Setting functions

  // Facets
  void set_facet_on_surface(const int facet, const bool f) {
    CGAL_assertion (facet>=0 && facet <4);
    bits[facet] = f;
//     if(f){
//       c_surface_facets |= (1 << facet);
//     }else {
//       c_surface_facets &= (15  & ~(1 << facet));
//     }
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Cb>()() + "+4b";
  }
#endif

  };  // end Complex_2_in_triangulation_cell_base_3

template < class GT, class Cb >
inline
std::istream&
operator>>(std::istream &is, Complex_2_in_triangulation_cell_base_3<GT, Cb> &c)
{
  bool b;
  is >> static_cast<Cb&>(c);
  for(int i = 0; i < 4; ++i)
  {
    if(is_ascii(is))
      is >> b;
    else
    {
      int i;
      read(is, i);
      b = static_cast<bool>(i);
    }
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
      os << ' ' << c.is_facet_on_surface(i);
    else
      write(os, static_cast<int>(c.is_facet_on_surface(i)));
  }
  return os;
}

}  // namespace CGAL



#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
