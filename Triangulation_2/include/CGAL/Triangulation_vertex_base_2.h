// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>

namespace CGAL {

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
#ifndef CGAL_NO_DEPRECATED_CODE
  typedef typename Tds::Vertex_circulator       Vertex_circulator; // needed for degree
  typedef typename Tds::Edge_circulator         Edge_circulator;   // needed for degree
  typedef typename Tds::Face_circulator         Face_circulator;   // needed for degree
  typedef typename Tds::size_type               size_type;         // needed for degree
#endif
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
  Triangulation_vertex_base_2(Face_handle f) : Vb(f) {} 

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

#ifndef CGAL_NO_DEPRECATED_CODE
  size_type degree(); //should be const

  Vertex_circulator incident_vertices()     
    {return Vertex_circulator(handle());}
 
  Vertex_circulator incident_vertices( Face_handle f)  
    {return Vertex_circulator(handle(),f);}
  
  Face_circulator incident_faces()  
    { return Face_circulator(handle()) ;}
  
  Face_circulator incident_faces( Face_handle f)    
    { return Face_circulator(handle(), f);}
  
  Edge_circulator incident_edges()   
    { return Edge_circulator(handle());}
  
  Edge_circulator incident_edges( Face_handle f)  
    { return Edge_circulator(handle(), f);}
 
  Vertex_handle handle();

#endif //  CGAL_NO_DEPRECATED_CODE
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

#ifndef CGAL_NO_DEPRECATED_CODE
template <class GT, class Vb>
typename Triangulation_vertex_base_2 <GT, Vb>::size_type
Triangulation_vertex_base_2 <GT, Vb>::
degree() //const
{
  int count = 0;
  Vertex_circulator vc = incident_vertices(), done(vc);
  if ( ! vc.is_empty()) {
    do { 
      count += 1;
    } while (++vc != done);
  }
  return count;
}
template <class GT, class Vb>
typename Triangulation_vertex_base_2<GT, Vb>::Vertex_handle
Triangulation_vertex_base_2 <GT, Vb> ::
handle()
{
  Face_handle fh = this->face();
  for(int i = 0 ; i < 3 ; ++i){
    if ( &*fh->vertex(i) == this) return fh->vertex(i);
  }
  return Vertex_handle();				    
}
#endif
} //namespace CGAL

#endif //CGAL_TRIANGULATION_VERTEX_BASE_2_H
