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
// file          : include/CGAL/Triangulation_face_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_FACE_2_H
#define CGAL_TRIANGULATION_FACE_2_H

#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_vertex_2;

template < class Gt, class Tds >
class Triangulation_vertex_handle_2;

template < class Gt, class Tds >
class Triangulation_face_handle_2;

template < class Gt, class Tds >
class Triangulation_face_2  : public  Tds::Face
{
public:
  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;

  typedef Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef Triangulation_face_2<Gt,Tds> Face;

  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  //  typedef std::pair<Face_handle, int>     Edge;


  Triangulation_face_2()
    : Fa()
  { }

  Triangulation_face_2(const Vertex_handle& v0,
		       const Vertex_handle& v1,
		       const Vertex_handle& v2)
    : Fa(&(*v0), &(*v1), &(*v2))
  {}
        
  Triangulation_face_2(const Vertex_handle& v0,
		       const Vertex_handle& v1,
		       const Vertex_handle& v2,
		       const Face_handle& n0,
		       const Face_handle& n1,
		       const Face_handle& n2)
    : Fa(&(*v0), &(*v1), &(*v2),&(*n0), &(*n1), &(*n2)) 
  {}
        
    

  // Vertex access functions
  inline Vertex_handle vertex(int i) const
  {
    return  ((Vertex *)(Fa::vertex(i)));
  }
    
  inline Vertex_handle opposite_vertex(int i) const
    {
      return ((Vertex *)(Fa::opposite_vertex(i)));
    }
 
  inline bool has_vertex(const Vertex_handle& v) const
  {
        return (Fa::has_vertex( & (*v)) );
  }
    
    
  inline bool has_vertex(const Vertex_handle& v, int& i) const
  {
    return Fa::has_vertex( &(*v), i);
  }

  inline int index(const Vertex_handle& v) const
  {
    return Fa::index( &(*v));
  }
  
  
 

  //ACCESS FUNCTIONS
  inline
  Face_handle neighbor(int i) const
  {
    return (Face *)(Fa::neighbor(i));
  }

   inline int index(const Face_handle& f) const
  {
    return Fa::index( &(*f));
  }
  
  inline bool has_neighbor(const Face_handle& f) const
  {
    return Fa::has_neighbor( &(*f));
  }

  inline bool has_neighbor(const Face_handle& f, int &i) const
  {
    return Fa::has_neighbor( &(*f), i);
  }

  inline int opposite_index(int i) const
    {
      return Fa::opposite_index(i);
    }

  inline Face_handle handle() const
  {
    return Face_handle(this);
  }

 //Setting
  inline
  void set_vertices(const Vertex_handle& v0,
		    const Vertex_handle& v1,
		    const Vertex_handle& v2)
    {
        Fa::set_vertices(&(*v0), &(*v1), &(*v2));
    }
    
  inline
    void set_neighbors(const Face_handle& n0,
                       const Face_handle& n1,
                       const Face_handle& n2)
    {
        Fa::set_neighbors(&(*n0), &(*n1), &(*n2));
    }

  inline  
  void set_vertices() 
  {
    Fa::set_vertices();
  }
   
 inline
  void set_neighbors() 
  {
    Fa::set_neighbors();
  }
    
  inline
    void set_vertex(int i, const Vertex_handle& v)
    {
        Fa::set_vertex(i, &(*v));
    }

    inline
    void set_neighbor(int i, const Face_handle& n)
    {
        Fa::set_neighbor(i, &(*n));
    }

};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_FACE_2_H
