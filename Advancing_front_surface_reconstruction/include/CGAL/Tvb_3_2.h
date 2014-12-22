// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $Source: /CVSROOT/CGAL/Packages/Triangulation_2/include/CGAL/Triangulation_vertex_base_2.h,v $
// $Revision$ $Date$
// $Name: current_submission $
//
// Author(s)     : Mariette Yvinec


#ifndef CGAL_TVB_3_2_H
#define CGAL_TVB_3_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>

namespace CGAL {

template < typename GT,
           typename Vb = Triangulation_ds_vertex_base_2<> >
class Tvb_3_2 
  : public Vb

{
  typedef typename Vb::Triangulation_data_structure    Tds;
public:
  typedef GT                                    Geom_traits;
  typedef typename GT::Point_3                  Point;
  typedef Tds                                   Triangulation_data_structure;
  typedef typename Tds::Face_handle             Face_handle;
  typedef typename Tds::Vertex_handle           Vertex_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Tvb_3_2<GT, Vb2>           Other;
  };

private:
  Point _p;

public:
  Tvb_3_2 () : Vb() {}
  Tvb_3_2(const Point & p) : Vb(), _p(p) {}
  Tvb_3_2(const Point & p, Face_handle f)
    : Vb(f), _p(p) {}
  Tvb_3_2(Face_handle f) : Vb(f) {} 

  void set_point(const Point & p) { _p = p; }
  const Point&  point() const { return _p; }

  // the non const version of point() is undocument
  // but needed to make the point iterator works
  // using Lutz projection scheme
  Point&        point() { return _p; }

  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
    {return true;}
};

template < class GT, class Vb >
std::istream&
operator>>(std::istream &is, Tvb_3_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return is >> static_cast<Vb&>(v) >> v.point();
}

template < class GT, class Vb >
std::ostream&
operator<<(std::ostream &os, const Tvb_3_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return os << static_cast<const Vb&>(v) << v.point();
}



} //namespace CGAL

#endif //CGAL_TVB_3_2_H
