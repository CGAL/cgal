// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Eti Ezra          <estere@post.tau.ac.il>
#ifndef CGAL_HOLES_SPLIT_DCEL_H
#define CGAL_SOLES_SPLIT_DCEL_H 

#include <CGAL/Bop_default_dcel.h>

CGAL_BEGIN_NAMESPACE

template <class Face_base>
class Holes_split_face : public  Face_bop<Face_base> {
public:
  typedef Face_bop<Face_base>         face_base;
  typedef Holes_split_face            face_bop;
  typedef const face_bop              const_face_bop;
  typedef const_face_bop*             const_pointer;
  typedef const_face_bop&             const_ref;
  
  Holes_split_face() : face_base(), in_hole_(false) {}

  Holes_split_face(const_pointer f) : face_base(*f), in_hole_(false) {}

  Holes_split_face(const_ref f) : face_base(f), in_hole_(false) {}

  virtual ~Holes_split_face() {}
  
  const_ref operator=(const_ref f) 
  {
    if (this == &f)
      return *this;

    face_base::assign(f);
    in_hole_ = f.in_hole_;
    
    return *this;
  }

  void set_in_hole(bool b) { in_hole_ = b ; }
  
  bool in_hole() const { return in_hole_; }

  virtual void assign(const_ref f)
  {
    operator=(f);  
  }
  
private:
  bool in_hole_;
};

template <class Halfedge_base>
class Holes_split_halfedge : public Halfedge_bop<Halfedge_base>
{
public:
  typedef Halfedge_bop<Halfedge_base>   halfedge_base;
  typedef Holes_split_halfedge          halfedge;
  typedef const halfedge                const_halfedge;
  typedef const_halfedge*               const_pointer;
  typedef const_halfedge&               const_ref;
  
  Holes_split_halfedge() : 
    halfedge_base(), 
    decomposing_(false) {}
  
  Holes_split_halfedge(const_pointer e) : 
    halfedge_base(*e),
    decomposing_(e->decomposing_) {}

  Holes_split_halfedge(const_ref e) : 
    halfedge_base(e) ,
    decomposing_(e.decomposing_) {} 
  
  virtual ~Holes_split_halfedge() {}
  
  const_ref operator=(const_ref e) 
  {
    halfedge_base::assign(e);
    
    decomposing_ = e.decomposing_;
    
    return *this;
  }
  
  virtual void assign(const_ref e)
  {
    operator=(e);  
  }
  
  void set_decomposing(bool b) { decomposing_ = b; }
  
  bool decomposing() const { return decomposing_; } 
  
  
protected:
  bool decomposing_;
};

template <class Traits, 
  class Vertex_base = Pm_vertex_base<typename Traits::Point>, 
  class Halfedge_base = Pm_halfedge_base<typename Traits::X_curve> , 
  class Face_base =  Pm_face_base> 
class Holes_split_dcel : public Pm_dcel<Vertex_bop<Vertex_base> , 
                                        Holes_split_halfedge<Halfedge_base>, 
                                        Holes_split_face<Face_base> > 
{
public:  // CREATION
  Holes_split_dcel() {} 
};

CGAL_END_NAMESPACE

#endif








