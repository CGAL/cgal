// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : Triangulation/include/CGAL/Triangulation_vertex_base_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>

CGAL_BEGIN_NAMESPACE

template < typename GT,
           typename Vb = Triangulation_ds_vertex_base_2<> >
class Triangulation_vertex_base_2 
  : public Vb

{
  typedef typename Vb::Triangulation_data_structure    Tds;
public:
  typedef GT                                    Geom_traits;
  typedef typename GT::Point_2                  Point;
  typedef Tds                                   Triangulation_data_structure;
  typedef typename Tds::Face_handle             Face_handle;
  typedef typename Tds::Vertex_handle           Vertex_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Triangulation_vertex_base_2<GT, Vb2>           Other;
  };

private:
  Point _p;

public:
  Triangulation_vertex_base_2 () : Vb() {}
  Triangulation_vertex_base_2(const Point & p) : Vb(), _p(p) {}
  Triangulation_vertex_base_2(const Point & p, Face_handle f)
    : Vb(f), _p(p) {}

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
operator>>(std::istream &is, Triangulation_vertex_base_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return is >> static_cast<Vb&>(v) >> v.point();
}

template < class GT, class Vb >
std::ostream&
operator<<(std::ostream &os, const Triangulation_vertex_base_2<GT, Vb> &v)
  // non combinatorial information. Default = point
{
  return os << static_cast<const Vb&>(v) << v.point();
}



CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_VERTEX_BASE_2_H
