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
// source        : $Source$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_2_H
#define CGAL_TRIANGULATION_VERTEX_2_H

#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_circulators_2.h>


template < class Gt, class Tds >
class CGAL_Triangulation_face_2;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_handle_2;

template < class Gt, class Tds >
class CGAL_Triangulation_face__handle_2;

template<class Gt,  class Tds>
class CGAL_Triangulation_face_circulator_2;

template<class Gt,  class Tds>
class CGAL_Triangulation_vertex_circulator_2;

template<class Gt, class Tds>
class CGAL_Triangulation_edge_circulator_2;


template<class Gt, class Tds >
class CGAL_Triangulation_vertex_2
  : public Tds::Vertex
{
public:
  
  typedef Tds Triangulation_data_structure;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;

  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  
  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  
  CGAL_Triangulation_vertex_2()
     : Ve()
  {}

  CGAL_Triangulation_vertex_2(const Point & p)
    :  Ve(p)
  {}
    
  CGAL_Triangulation_vertex_2(const Point & p, const Face_handle& f)
    :  Ve(p, &(*f))
  {}

  inline void set_face(const Face_handle& f)
  {
    Ve::set_face(&(*f));
  }
    
    
  inline Face_handle face() const
  {
        return  (Face *)Ve::face();
  }
        
  inline Vertex_handle handle() const
  {
    return Vertex_handle(this);
  }
      
  //inherited :
  // cw(i), ccw(i)
  // degree(), point(), is_on_boundary, set_point(p)

  inline Vertex_circulator incident_vertices() const
  {
    return Vertex_circulator(handle(), face());
  }

  inline Vertex_circulator incident_vertices(const Face_handle& f) const
  {
    return Vertex_circulator(handle(), f);
  } 

  inline 
  Face_circulator incident_faces() const
  {
    return Face_circulator(handle(), face());
  }
    
  inline Face_circulator incident_faces(const Face_handle& f) const
  {
    return Face_circulator(handle(), f);
  }
    
  inline Edge_circulator incident_edges() const
  {
    return Edge_circulator(handle(), face());
  }
 
  inline Edge_circulator incident_edges(const Face_handle& f) const
  {
    return Edge_circulator(handle(), f);
  }
    
};
#endif CGAL_TRIANGULATION_2_H
