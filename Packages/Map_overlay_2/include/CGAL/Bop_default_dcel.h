#ifndef CGAL_BOP_DEFAULT_DCEL_H
#define CGAL_BOP_DEFAULT_DCEL_H 

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <fstream.h>
#include <vector>
#include <list>
#include <string>

//#ifndef CGAL_ARRANGEMENT_2_H
//#include <CGAL/Arrangement_2.h>
//#endif

#ifndef CGAL_ARR_2_OVERLAY_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

CGAL_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////
//               OVERLAY DCEL
///////////////////////////////////////////////////////////////

template <class Face_base>
class Arr_2_face_bop : public  Arr_2_face_overlay<Face_base> {
public:
  typedef Arr_2_face_overlay<Face_base>     face_overlay;
  typedef Arr_2_face_bop                    face_bop;
  typedef const face_bop                    const_face_bop;
  typedef const_face_bop*                   const_pointer;
  typedef const_face_bop&                   const_ref;
  
  //enum COLOR {WHITE=0, GRAY=1, BLACK=2};
  
  Arr_2_face_bop() : face_overlay() { bop_ = true; }

  Arr_2_face_bop(const_pointer f) : face_overlay(*f) { set_bop(f->bop()); }

  Arr_2_face_bop(const_ref f) : face_overlay(f) { set_bop(f.bop()); } 
  
  const_ref operator=(const_ref f) 
  {
    face_overlay::assign(f);
    bop_ = f.bop();
    
    return *this;
  }

  void set_bop(bool b) { bop_ = b; }

  bool bop() const { return bop_; }

  virtual void assign(const_ref f)
  {
    operator=(f);  
  }
  
private:
  mutable bool bop_;
};

//template <class Base_node>
template <class Halfedge_base>
class Arr_2_halfedge_bop : public Arr_2_halfedge_overlay<Halfedge_base>
{
public:
  typedef Arr_2_halfedge_overlay<Halfedge_base>   halfedge_overlay;
  typedef Arr_2_halfedge_bop                      halfedge_bop;
  typedef const halfedge_bop                      const_halfedge_bop;
  typedef const_halfedge_bop*                     const_pointer;
  typedef const_halfedge_bop&                     const_ref;
  
  Arr_2_halfedge_bop() : halfedge_overlay() { bop_ = true; }

  Arr_2_halfedge_bop(const_pointer e) : halfedge_overlay(*e) { set_bop(e->bop()); }

  Arr_2_halfedge_bop(const_ref e) : halfedge_overlay(e) { set_bop(e.bop()); } 
  
  const_ref operator=(const_ref e) 
  {
    halfedge_overlay::assign(e);
    bop_ = e.bop();
    
    return *this;
  }

  void set_bop(bool b) { bop_ = b; }

  bool bop() const { return bop_; }

  virtual void assign(const_ref e)
  {
    operator=(e);  
  }
  
private:
  mutable bool bop_;
};


//template <class Point>
template <class Vertex_base>
class Arr_2_vertex_bop : public Arr_2_vertex_overlay<Vertex_base>
{
public:
  typedef Arr_2_vertex_overlay<Vertex_base>  vertex_overlay;
  typedef Arr_2_vertex_bop                   vertex_bop;
  typedef const vertex_bop                   const_vertex_bop;
  typedef const_vertex_bop*                  const_pointer;
  typedef const_vertex_bop&                  const_ref;

  Arr_2_vertex_bop() : vertex_overlay (){ bop_ = true; }

  Arr_2_vertex_bop(const_pointer v) : vertex_overlay(*v) { set_bop(v->bop()); }

  Arr_2_vertex_bop(const_ref v) : vertex_overlay(v) { set_bop(v.bop()); } 
  
  const_ref operator=(const_ref v) 
  {
    vertex_overlay::assign(v);
    bop_ = v.bop();
    
    return *this;
  }

  void set_bop(bool b) { bop_ = b; }
  
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
class Bop_default_dcel: public Pm_dcel<Arr_2_vertex_bop<Vertex_base>, Arr_2_halfedge_bop<Halfedge_base>, Arr_2_face_bop<Face_base> > 
{
public:  // CREATION
  Bop_default_dcel() {} 
};

CGAL_END_NAMESPACE

#endif









