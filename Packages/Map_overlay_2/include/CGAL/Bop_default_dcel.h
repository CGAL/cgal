#ifndef CGAL_BOP_DEFAULT_DCEL_H
#define CGAL_BOP_DEFAULT_DCEL_H 

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <fstream.h>
#include <vector>
#include <list>
#include <string>

#ifndef CGAL_MAP_OVERLAY_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

CGAL_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////
//               OVERLAY DCEL
///////////////////////////////////////////////////////////////

template <class Face_base>
class Face_bop : public  Face_overlay<Face_base> {
public:
  typedef Face_overlay<Face_base>     face_overlay;
  typedef Face_bop                    face_bop;
  typedef const face_bop              const_face_bop;
  typedef const_face_bop*             const_pointer;
  typedef const_face_bop&             const_ref;
  
  Face_bop() : face_overlay() { bop_ = true; }

  Face_bop(const_pointer f) : face_overlay(*f) { set_ignore_bop(f->bop()); }

  Face_bop(const_ref f) : face_overlay(f) { set_ignore_bop(f.bop()); } 

  virtual ~Face_bop() {}
  
  const_ref operator=(const_ref f) 
  {
    if (this == &f)
      return *this;

    face_overlay::assign(f);
    bop_ = f.bop();
    
    return *this;
  }

  void set_ignore_bop(bool b) { bop_ = b; }

  bool bop() const { return bop_; }

  virtual void assign(const_ref f)
  {
    operator=(f);  
  }
  
private:
  mutable bool bop_;
};

template <class Halfedge_base>
class Halfedge_bop : public Halfedge_overlay<Halfedge_base>
{
public:
  typedef Halfedge_overlay<Halfedge_base>   halfedge_overlay;
  typedef Halfedge_bop                      halfedge_bop;
  typedef const halfedge_bop                      const_halfedge_bop;
  typedef const_halfedge_bop*                     const_pointer;
  typedef const_halfedge_bop&                     const_ref;
  
  Halfedge_bop() : halfedge_overlay() { bop_ = true; }

  Halfedge_bop(const_pointer e) : halfedge_overlay(*e) { set_ignore_bop(e->bop()); }

  Halfedge_bop(const_ref e) : halfedge_overlay(e) { set_ignore_bop(e.bop()); } 
  
  virtual ~Halfedge_bop() {}
  
  const_ref operator=(const_ref e) 
  {
    if (this == &e)
      return *this;

    halfedge_overlay::assign(e);
    bop_ = e.bop();
    
    return *this;
  }

  void set_ignore_bop(bool b) { bop_ = b; }

  bool bop() const { return bop_; }

  virtual void assign(const_ref e)
  {
    operator=(e);  
  }
  
private:
  mutable bool bop_;
};


template <class Vertex_base>
class Vertex_bop : public Vertex_overlay<Vertex_base>
{
public:
  typedef Vertex_overlay<Vertex_base>  vertex_overlay;
  typedef Vertex_bop                   vertex_bop;
  typedef const vertex_bop                   const_vertex_bop;
  typedef const_vertex_bop*                  const_pointer;
  typedef const_vertex_bop&                  const_ref;

  Vertex_bop() : vertex_overlay (){ bop_ = true; }

  Vertex_bop(const_pointer v) : vertex_overlay(*v) { set_ignore_bop(v->bop()); }

  Vertex_bop(const_ref v) : vertex_overlay(v) { set_ignore_bop(v.bop()); } 
  
  virtual ~Vertex_bop() {}
  
  const_ref operator=(const_ref v) 
  {
    if (this == &v)
      return *this;

    vertex_overlay::assign(v);
    bop_ = v.bop();
    
    return *this;
  }

  void set_ignore_bop(bool b) { bop_ = b; }
  
  bool bop() const { return bop_; }

  virtual void assign(const_ref v)
  {
    operator=(v);  
  }

private:
  mutable bool bop_;
};

template <class Traits, class Vertex_base = Pm_vertex_base<typename Traits::Point>, 
  class Halfedge_base = Pm_halfedge_base<typename Traits::X_curve> , 
  class Face_base =  Pm_face_base> 
class Bop_default_dcel: public Pm_dcel<Vertex_bop<Vertex_base>, Halfedge_bop<Halfedge_base>, Face_bop<Face_base> > 
{
public:  // CREATION
  Bop_default_dcel() {} 
};

CGAL_END_NAMESPACE

#endif









