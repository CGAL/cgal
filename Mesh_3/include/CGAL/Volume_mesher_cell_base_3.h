// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU


#ifndef CGAL_VOLUME_MESHER_CELL_BASE_3_H
#define CGAL_VOLUME_MESHER_CELL_BASE_3_H

#include <CGAL/Triangulation_ds_cell_base_3.h>

#include <string>

namespace CGAL {

  template < class GT, class Cb = Triangulation_ds_cell_base_3 <> >
  class Volume_mesher_cell_base_3 : public Cb {

  public:
    typedef Volume_mesher_cell_base_3 <GT, Cb> Self;

    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS3>::Other  Cb3;
      typedef Volume_mesher_cell_base_3 <GT, Cb3> Other;
    };


    typedef typename GT::Point_3 Point;

    typedef typename Cb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;


  private:
    bool domain;

  public:

    // Constructors

    Volume_mesher_cell_base_3() : Cb(), domain(false)
    {
    }

    Volume_mesher_cell_base_3 (Vertex_handle v0,
					    Vertex_handle v1,
					    Vertex_handle v2,
					    Vertex_handle v3)
      :  Cb (v0, v1, v2, v3), domain(false)
    {
    }

    Volume_mesher_cell_base_3 (Vertex_handle v0,
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

#ifdef CGAL_MESH_3_IO_H
    static
    std::string io_signature()
    {
      return Get_io_signature<Cb>()() + "+b";
    }
#endif
  };  // end Volume_mesher_cell_base_3

template < class GT, class Cb >
inline
std::istream&
operator>>(std::istream &is, Volume_mesher_cell_base_3<GT, Cb> &c)
{
  bool b;
  is >> static_cast<Cb&>(c);
  is >> b;
  c.set_in_domain(b);
  return is;
}

template < class GT, class Cb >
inline
std::ostream&
operator<<(std::ostream &os,
           const Volume_mesher_cell_base_3<GT, Cb> &c)
{
  os << static_cast<const Cb&>(c);
  if(is_ascii(os))
    os << ' ';
  return os << c.is_in_domain();
}

}  // namespace CGAL


#endif  // CGAL_VOLUME_MESHER_CELL_BASE_3_H
