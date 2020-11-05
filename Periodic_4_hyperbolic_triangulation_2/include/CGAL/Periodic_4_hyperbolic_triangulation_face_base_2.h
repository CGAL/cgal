// Copyright (c) 1999-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2
#define CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/basic.h>
#include <CGAL/Dummy_tds_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {

template< typename GT,
          typename FB = Triangulation_face_base_2<GT> >
class Periodic_4_hyperbolic_triangulation_face_base_2
    : public FB
{
public:
  typedef typename FB::Vertex_handle                  Vertex_handle;
  typedef typename FB::Face_handle                    Face_handle;
  typedef typename GT::Hyperbolic_translation         Hyperbolic_translation;

  template< typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename FB::template Rebind_TDS<TDS2>::Other              FB2;
    typedef Periodic_4_hyperbolic_triangulation_face_base_2<GT, FB2>   Other;
  };

private:
  Hyperbolic_translation o[3]; // Hyperbolic_translations for vertices

public:
  Periodic_4_hyperbolic_triangulation_face_base_2()
    : FB()
  {
    o[0] = Hyperbolic_translation();
    o[1] = Hyperbolic_translation();
    o[2] = Hyperbolic_translation();
  }

  Periodic_4_hyperbolic_triangulation_face_base_2(const Vertex_handle& v0,
                                                  const Vertex_handle& v1,
                                                  const Vertex_handle& v2)
    : FB(v0, v1, v2)
  {
    o[0] = Hyperbolic_translation();
    o[1] = Hyperbolic_translation();
    o[2] = Hyperbolic_translation();
  }

  Periodic_4_hyperbolic_triangulation_face_base_2(const Vertex_handle& v0, const Vertex_handle& v1,
                                                  const Vertex_handle& v2, const Face_handle& n0,
                                                  const Face_handle& n1, const Face_handle& n2)
    : FB(v0,v1,v2,n0,n1,n2)
  {
    o[0] = Hyperbolic_translation();
    o[1] = Hyperbolic_translation();
    o[2] = Hyperbolic_translation();
  }

  Hyperbolic_translation translation(int i) const
  {
    CGAL_triangulation_precondition(i >= 0 && i <= 2);
    return o[i];
  }

  void set_translation(const int& k, const Hyperbolic_translation& new_o)
  {
    CGAL_triangulation_precondition(k >= 0 && k <= 2);
    o[k] = new_o;
  }

  // CHECKING

  void reorient()
  {
    // swaps vertex 0 with vertex 1, and neighbor 0 with neighbor 1
    FB::reorient();

    // manually swap translation 0 with translation 1
    Hyperbolic_translation tmp = o[0];
    o[0] = o[1];
    o[1] = tmp;
  }
};

template < class TDS >
inline
std::istream&
operator>>(std::istream &is, Periodic_4_hyperbolic_triangulation_face_base_2<TDS> &)
// non combinatorial information. Default = nothing
{
  return is;
}

template < class TDS >
inline
std::ostream&
operator<<(std::ostream &os,
           const Periodic_4_hyperbolic_triangulation_face_base_2<TDS> &)
// non combinatorial information. Default = nothing
{
  return os;
}

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TIANGULATION_FACE_BASE_2
