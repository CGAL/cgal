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
// file          : include/CGAL/Triangulation_vertex_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_2_H
#define CGAL_TRIANGULATION_VERTEX_2_H

#include <utility>
#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_circulators_2.h>

CGAL_BEGIN_NAMESPACE 

template < class Gt, class Tds >
class Triangulation_face_2;

template < class Gt, class Tds >
class Triangulation_vertex_handle_2;

template < class Gt, class Tds >
class Triangulation_face__handle_2;

template<class Gt,  class Tds>
class Triangulation_face_circulator_2;

template<class Gt,  class Tds>
class Triangulation_vertex_circulator_2;

template<class Gt, class Tds>
class Triangulation_edge_circulator_2;


template<class Gt, class Tds >
class Triangulation_vertex_2  : public Tds::Vertex
{
public:
  
  typedef Tds Triangulation_data_structure;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point_2 Point;
  typedef typename Geom_traits::Segment_2 Segment;
  typedef typename Geom_traits::Triangle_2 Triangle;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;

  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Tr_vertex;
  
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef std::pair<Face_handle, int>     Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  
  Triangulation_vertex_2()     : Ve()   {}
  Triangulation_vertex_2(const Point & p)    :  Ve(p)   {}
  Triangulation_vertex_2(const Point & p, const Face_handle& f)
    :  Ve(p, &(*f))
  {}

  void set_face(const Face_handle& f);
  Face_handle face() const;
  Vertex_handle handle() const;
        
  //inherited :
  // cw(i), ccw(i)
  // degree(), point(), is_on_boundary, set_point(p)

  Vertex_circulator incident_vertices() const;
  Vertex_circulator incident_vertices(Face_handle f) const;
  Face_circulator incident_faces() const;
  Face_circulator incident_faces(Face_handle f) const;
  Edge_circulator incident_edges() const;
  Edge_circulator incident_edges(Face_handle f) const;
   
};

template<class Gt, class Tds>
inline void
Triangulation_vertex_2<Gt,Tds>::
set_face(const Face_handle& f)
{
  Ve::set_face(&(*f));
}
    
template<class Gt, class Tds>
inline   
Triangulation_face_handle_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>::
face() const
{
  return  (Face *)Ve::face();
}
        
template<class Gt, class Tds>
inline   
Triangulation_vertex_handle_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>::
handle() const
{
  Tr_vertex * ncthis = (Tr_vertex *) this;
  return Vertex_handle(ncthis);
}
      
template<class Gt, class Tds>
inline   
Triangulation_vertex_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>:: 
incident_vertices() const
{
  return Vertex_circulator(handle(), face());
}

template<class Gt, class Tds>
inline   
Triangulation_vertex_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>::
incident_vertices(Face_handle f) const
{
  return Vertex_circulator(handle(), f);
} 

template<class Gt, class Tds>
inline   
Triangulation_face_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>::
incident_faces() const
{
  return Face_circulator(handle(), face());
}
    
template<class Gt, class Tds>
inline   
Triangulation_face_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>::
incident_faces(Face_handle f) const
{
  return Face_circulator(handle(), f);
}
    
template<class Gt, class Tds>
inline   
Triangulation_edge_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>:: 
incident_edges() const
{
  return Edge_circulator(handle(), face());
}
 
template<class Gt, class Tds>
inline   
Triangulation_edge_circulator_2<Gt,Tds>
Triangulation_vertex_2<Gt,Tds>:: 
incident_edges(Face_handle f) const
{
  return Edge_circulator(handle(), f);
}

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_VERTEX_2_H
