// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Map_overlay_default_dcel.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_MAP_OVERLAY_DEFAULT_DCEL_H
#define CGAL_MAP_OVERLAY_DEFAULT_DCEL_H 

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <fstream.h>
#include <vector>
#include <list>
#include <string>

#ifndef CGAL_PM_DEFAULT_DCEL_H
#include<CGAL/Pm_default_dcel.h>
#endif

CGAL_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////
//               OVERLAY DCEL
///////////////////////////////////////////////////////////////

template <class Face_base>
class Face_overlay : public Face_base {
public:
  typedef Face_base                       face_base;
  typedef Face_overlay                    face_overlay;
  typedef const face_overlay              const_face_overlay;
  typedef const_face_overlay*             const_pointer;
  typedef const_face_overlay&             const_ref;
  typedef void*                           Face_pointer;
  
  enum COLOR {WHITE=0, GRAY=1, BLACK=2};
  
  Face_overlay() : face_base()
  { 
    faces_above[0] = faces_above[1] = 0; 
    color = WHITE;
  }

  Face_overlay(const_pointer f) : face_base(*f) 
  {
    set_first_face_above(f->get_first_face_above());
    set_second_face_above(f->get_second_face_above());
    set_color(f->get_color());
  }

  Face_overlay(const_ref f) : face_base(f) 
  { 
    set_first_face_above(f.get_first_face_above());
    set_second_face_above(f.get_second_face_above());
    set_color(f.get_color());
  } 

  virtual ~Face_overlay() {}
  
  const_ref operator=(const_ref f) 
  {
    if (this == &f)
      return *this;
    
    face_base::assign(f);
    reset();
    
    set_first_face_above(f.get_first_face_above());
    set_second_face_above(f.get_second_face_above());
    set_color(f.get_color());

    return *this;
  }

  void  set_first_face_above(const void* const face)
  {
    // updating only if the value is not NULL and not already in the array.
    if (face && faces_above[0] != face && faces_above[1] != face){
      faces_above[0] = face;
    }
  }

  void  set_second_face_above(const void* const face)
  {
    // updating only if the value is not NULL and not already in the array.
    if (face && faces_above[0] != face && faces_above[1] != face){
      faces_above[1] = face;
    }
  }
  
  const void* const get_first_face_above() const
    {
      return  faces_above[0];
    }

  const void* const get_second_face_above() const 
    {
      return  faces_above[1];
    }

  void  reset()
  {
    faces_above[0] = faces_above[1] = 0; 
  }
  
  void set_color(COLOR c) { color = c;}
  
  COLOR get_color() const { return color;}

  virtual void assign(const_ref f)
  {
    operator=(f);  
  }

protected:
  const void*   faces_above[2];
  COLOR         color;
};

//template <class X_curve>
//class Pm_halfedge_overlay : public Pm_halfedge_base<X_curve>

template <class Halfedge_base>
class Halfedge_overlay : public Halfedge_base
{
public:
  typedef Halfedge_base                   halfedge_base;
  typedef Halfedge_overlay                halfedge_overlay;
  typedef const halfedge_overlay          const_halfedge_overlay;
  typedef const_halfedge_overlay*         const_pointer;
  typedef const_halfedge_overlay&         const_ref;
  typedef void*                           Halfedge_pointer;
  
  Halfedge_overlay() : halfedge_base()
  {
    halfedges_above[0] = halfedges_above[1] = 0;
    faces_above[0] = faces_above[1] = 0;
  }

  Halfedge_overlay(const_pointer e) : halfedge_base(*e) 
  {  
    set_first_halfedge_above(e->get_first_halfedge_above());
    set_second_halfedge_above(e->get_second_halfedge_above());

    set_first_face_above(e->get_first_face_above());
    set_second_face_above(e->get_second_face_above());

    //halfedges_above[0] = halfedges_above[1] = 0;
    //faces_above[0] = faces_above[1] = 0;
  }

  Halfedge_overlay(const_ref e) : halfedge_base(e) 
  { 
    set_first_halfedge_above(e.get_first_halfedge_above());
    set_second_halfedge_above(e.get_second_halfedge_above());
    
    set_first_face_above(e.get_first_face_above());
    set_second_face_above(e.get_second_face_above());

    //faces_above[0] = faces_above[1] = 0;
  } 
 
  virtual ~Halfedge_overlay() {}
  
  const_ref operator=(const_ref e) 
  {
    if (this == &e)
      return *this;

    halfedge_base::assign(e);

    reset();
    
    set_first_halfedge_above(e.get_first_halfedge_above());
    set_second_halfedge_above(e.get_second_halfedge_above());
    
    set_first_face_above(e.get_first_face_above());
    set_second_face_above(e.get_second_face_above());

    return *this;
  }

  void  set_first_halfedge_above(const void* const halfedge)
  {
    // updating only if the value is not 0 and not already in the array.
    if (halfedge &&  halfedges_above[0] != halfedge && 
        halfedges_above[1] != halfedge)
      {
        halfedges_above[0] = halfedge;
      }
  }
  
   void  set_second_halfedge_above(const void* const halfedge)
  {
    // updating only if the value is not 0 and not already in the array.
    if (halfedge &&  halfedges_above[0] != halfedge && 
        halfedges_above[1] != halfedge)
      {
        halfedges_above[1] = halfedge;
      }
  }
  
  const void* const get_first_halfedge_above() const
    {
      return halfedges_above[0];
    }

  const void* const get_second_halfedge_above() const 
    {
      return halfedges_above[1];
    }

  // ----  add faces above halfedges (not to their sides).
  void  set_first_face_above(const void*  face)
  {
    // updating only if the value is not 0 and not already in the array.
    if (face &&  faces_above[0] != face && faces_above[1] != face)
      {
        faces_above[0] = face;
      }
  }
  
  void  set_second_face_above (const void*  face)
  {
    // updating only if the value is not 0 and not already in the array.
    if (face &&  faces_above[0] != face && faces_above[1] != face)
      {
        faces_above[1] = face;
      }
  }

  void reset()
  {
    reset_halfedges_above();
    reset_faces_above();
  }

  const void* const get_first_face_above() const
  {
    return faces_above[0];
  }

  const void* const get_second_face_above() const 
  {
    return faces_above[1];
  }
  
  virtual void assign(const_ref e)
  {
    operator=(e);  
  }

protected:
  
  void reset_halfedges_above()
  { 
    halfedges_above[0] = halfedges_above[1] = 0;
  }
  
  void reset_faces_above()
  { 
    faces_above[0] = faces_above[1] = 0;
  }

  const void*  halfedges_above[2];
  const void*  faces_above[2];
};


template <class Vertex_base>
class Vertex_overlay : public Vertex_base
{
public:
  typedef Vertex_base                   vertex_base;
  typedef Vertex_overlay                vertex_overlay;
  typedef const vertex_overlay          const_vertex_overlay;
  typedef const_vertex_overlay*         const_pointer;
  typedef const_vertex_overlay&         const_ref;
  typedef void*                         Vertex_pointer;

  Vertex_overlay() : vertex_base() {
    vertices_above[0] = vertices_above[1] = 0;
    halfedges_above[0] = halfedges_above[1] = 0;
    faces_above[0] = faces_above[1] = 0;
  }

  Vertex_overlay(const_pointer v) : vertex_base(*v) {
    set_first_vertex_above(v->get_first_vertex_above());
    set_second_vertex_above(v->get_second_vertex_above());
    
    set_first_halfedge_above(v->get_first_halfedge_above());
    set_second_halfedge_above(v->get_second_halfedge_above());

    set_first_face_above(v->get_first_face_above());
    set_second_face_above(v->get_second_face_above());
  }
  
  Vertex_overlay(const_ref v) : vertex_base(v) {  
    set_first_vertex_above(v.get_first_vertex_above());
    set_second_vertex_above(v.get_second_vertex_above());

    set_first_halfedge_above(v.get_first_halfedge_above());
    set_second_halfedge_above(v.get_second_halfedge_above());

    set_first_face_above(v.get_first_face_above());
    set_second_face_above(v.get_second_face_above());
  } 

  virtual ~Vertex_overlay() {}
  
  const_ref operator=(const_ref v) 
  {
    if (this == &v)
      return *this;

    vertex_base::assign(v);

    reset();
    
    set_first_vertex_above(v.get_first_vertex_above());
    set_second_vertex_above(v.get_second_vertex_above());

    set_first_halfedge_above(v.get_first_halfedge_above());
    set_second_halfedge_above(v.get_second_halfedge_above());

    set_first_face_above(v.get_first_face_above());
    set_second_face_above(v.get_second_face_above());

    return *this;
  }

  void  set_first_vertex_above(const void* const vertex) 
  {
    // updating only if the value is not 0 and not already in the array.
    if (vertex && vertices_above[0] != vertex && vertices_above[1] != vertex)
      {
        vertices_above[0] = vertex;
      }
  }

  void  set_second_vertex_above(const void* const vertex) 
  {
    // updating only if the value is not 0 and not already in the array.
    if (vertex && vertices_above[0] != vertex && vertices_above[1] != vertex)
      {
        vertices_above[1] = vertex;
      }
  }
  
  const void* const get_first_vertex_above() const 
  {
    return  vertices_above[0];
  }

  const void* const get_second_vertex_above() const 
  {
    return  vertices_above[1];
  }
  
  void  set_first_halfedge_above(const void* const halfedge) 
  {
    
    if (halfedge && halfedges_above[0] != halfedge && 
        halfedges_above[1] != halfedge){
      halfedges_above[0] = halfedge;
    }
  }
  
  void  set_second_halfedge_above(const void* const halfedge) 
  {
    if (halfedge && halfedges_above[0] != halfedge && 
        halfedges_above[1] != halfedge){
      halfedges_above[1] = halfedge;
    }
  }
  
  const void* const get_first_halfedge_above() const 
  {
    return halfedges_above[0];
  }

  const void* const get_second_halfedge_above() const 
  {
    return halfedges_above[1];
  }
  
  void  set_first_face_above(const void*  face) 
  {
    // updating only if the value is not 0 and not already in the array.
    if (face && faces_above[0] != face && faces_above[1] != face){
      faces_above[0] = face;
    }
  }
  
  void  set_second_face_above (const void*  face) 
  {
    // updating only if the value is not 0 and not already in the array.
    if (face && faces_above[0] != face && faces_above[1] != face){   
      faces_above[1] = face;
    }
  }

  const void* const get_first_face_above() const 
  {
    return faces_above[0];
  }
  
  const void* const get_second_face_above() const 
  {
    return faces_above[1];
  }

  void reset() 
  {
    vertices_above[0] = vertices_above[1] = 0; 
    halfedges_above[0] = halfedges_above[1] = 0; 
    faces_above[0] = faces_above[1] = 0; 
  }

  virtual void assign(const_ref v)
  {
    operator=(v);  
  }
  
protected:
  const void*    vertices_above[2];
  const void*    halfedges_above[2];
  const void*    faces_above[2];
};

template <class Traits, 
  class Vertex_base = Pm_vertex_base<typename Traits::Point>, 
  class Halfedge_base = Pm_halfedge_base<typename Traits::X_curve> , 
  class Face_base =  Pm_face_base> 
class Map_overlay_default_dcel : public Pm_dcel<Vertex_overlay<Vertex_base> , 
                                              Halfedge_overlay<Halfedge_base>, 
                                              Face_overlay<Face_base> > 
{
public:  // CREATION
  Map_overlay_default_dcel() {} 
};

CGAL_END_NAMESPACE

#endif













