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
// Author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_PM_DEFAULT_DCEL_H
#define CGAL_PM_DEFAULT_DCEL_H

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <list>
#include <map>
#include <CGAL/N_step_adaptor.h>
#include <CGAL/In_place_list.h>
#include <CGAL/HalfedgeDS_iterator.h>

CGAL_BEGIN_NAMESPACE

/*! Base vertex class */
template <class Pt>
class Pm_vertex_base {
protected:
  /*! An incident halfedge pointing at the vertex */
  void * hdg;
  Pt pt;

public:
  typedef Pt           Point;
  
  Pm_vertex_base() {}
  Pm_vertex_base(const Pt & p) : pt(p) {}
  
  virtual ~Pm_vertex_base() {}

  void * halfedge() { return hdg; }
  const void * halfedge() const { return hdg; }
  void set_halfedge(void * h) { hdg = h; }
  
  Point & point() { return pt; }
  const Point & point() const { return pt; }
  
  void set_point(const Point & p) { pt = p; }

  // assign function for non-connectivity data
  virtual void assign(const Pm_vertex_base<Pt> & v)
  {
    pt = v.pt;
  }
};

/*! Base halfedge class */
template <class X_curve>
class Pm_halfedge_base {
public:
  typedef X_curve Curve;

  /*! Parameterless Constructor */
  Pm_halfedge_base() {}

  /*! Constructor */
  Pm_halfedge_base(const X_curve & c) : cv(c) {}

  /*! Destructor */
  virtual ~Pm_halfedge_base() {}

  /*! \brief obtains the opposite halfedge */
  void * opposite() { return opp; }
  const void * opposite() const { return opp; }

  /*! \brief obtains the next halfedge along the face */
  void * next() { return nxt; }
  const void * next() const { return nxt; }
  
  /*! \brief sets the opposite halfedge */
  void set_opposite(void * h) { opp = h; }

  /*! \brief sets the next halfedge along the face */
  void set_next(void * h) { nxt = h; }
  
  /*! \brief obtains the incident target vertex */
  void * vertex() { return v; }
  const void * vertex() const { return v; }
  
  /*! \brief obtains the face to the left */
  void * face() { return f; }
  const void * face() const { return f; }

  /*! \brief sets the incident target vertex */
  void set_vertex(void * _v) { v = _v; }

  /*! \brief sets the face to the left */
  void set_face(void * _f) { f = _f; }

  /*! \brief obtains the geometric x-monotone curve */
  Curve & curve() { return cv; }
  const Curve & curve() const { return cv; }

  /*! \brief sets the geometric x-monotone curve */
  void set_curve(const X_curve& c) { cv = c; }

  /*! assign function for non-connectivity data */
  virtual void assign(const Pm_halfedge_base<X_curve> &e)
  {
    cv = e.cv;
  }

protected:
  void * opp;
  void * nxt;
  
  void * v;     // target
  void * f;     // face
  
  X_curve cv;
};

/*! Base face class */
class Pm_face_base {
public:
  typedef std::list<void*> Holes_container; 
  typedef Holes_container::iterator Holes_iterator; 
  typedef Holes_container::const_iterator Holes_const_iterator;

  /*! Parameterless Constructor */
  Pm_face_base() : holes() {};

  /*! Destructor */
  virtual ~Pm_face_base() {}

  void * halfedge() { return hdg;}
  const void * halfedge() const { return hdg;}

  void set_halfedge(void * h) {hdg = h;}

  Holes_iterator  holes_begin() {return holes.begin();}
  Holes_iterator  holes_end() {return holes.end();}

  Holes_const_iterator  holes_begin() const {return holes.begin();}
  Holes_const_iterator  holes_end() const {return holes.end();}

  void add_hole(void * halfedge_ptr) { holes.push_back(halfedge_ptr); }

  void erase_hole(Holes_iterator hit) { holes.erase(hit); }
    
  void erase_holes(Holes_iterator first, Holes_iterator last) {
    holes.erase(first,last);
  }

  // assign function for non-connectivity data
  virtual void assign(const Pm_face_base & f)
  { 
    // The reason we do not assign anything here is because
    // we can't copy pointers of a face from another Pm_dcel.
    // The assign function of Pm_dcel does all the necessary updates.

    (void) f; // We avoid the `unused parameter' warning.
  }
  
private:
  void * hdg;
  Holes_container holes;
};

// Forward declarations:
template <class V, class H, class F> class _Pm_Vertex;
template <class V, class H, class F> class _Pm_Halfedge;
template <class V, class H, class F> class _Pm_Face;

/*! The vertex class */
template <class V, class H, class F>
class _Pm_Vertex : public V,
                   public In_place_list_base< _Pm_Vertex<V,H,F> >
{
public:
  typedef V                     Base;
  typedef _Pm_Vertex<V,H,F>     Vertex;
  typedef _Pm_Halfedge<V,H,F>   Halfedge;
  typedef _Pm_Face<V,H,F>       Face;

  _Pm_Vertex() {}
  // _Pm_Vertex(const Point & p) : V(p) {}

  /*! \brief obtains an incident halfedge */
  Halfedge * halfedge() { return (Halfedge*)(V::halfedge()); }

  /*! \brief obtains an incident halfedge */
  const Halfedge * halfedge() const
  {
    return (const Halfedge*)(V::halfedge());
  }

  /*! \brief sets an incident halfedge */
  void set_halfedge(Halfedge * h) { V::set_halfedge(h); }

#if 0
  // Apparently the implementation of the In_place_list requires a copy
  // constructor of the value type!
protected:
  // forbid copy constructor and assignment (only derived classes can use them)

  /*! Copy Constructor */
  _Pm_Vertex(const _Pm_Vertex & v) {}

  /*! Assignment operator */
  _Pm_Vertex & operator=(const _Pm_Vertex &) { return *this; }
#endif
};

/*! The halfdege class */
template <class V, class H, class F>
class _Pm_Halfedge : public  H,
                     public  In_place_list_base< _Pm_Halfedge<V,H,F> >
{
public:
  typedef H                     Base;
  typedef _Pm_Vertex<V,H,F>     Vertex;
  typedef _Pm_Halfedge<V,H,F>   Halfedge;
  typedef _Pm_Face<V,H,F>       Face;

  _Pm_Halfedge() : H() {}  
  //_Pm_Halfedge( const Curve& c) : H(c) {}

  /*! \brief obtains the opposite halfedge */
  Halfedge * opposite() { return (Halfedge*)(H::opposite()); }
  
  //in the future will probably be implemented in a max base
  //    const Halfedge* prev()     const {return (const Halfedge*)(H::prev());}

  /*! \brief obtains the base halfedge */
  Vertex * vertex() { return (Vertex*)(H::vertex()); }

  /*! \brief obtains the base halfedge */
  const Vertex * vertex() const {return (const Vertex*)(H::vertex()); }

  /*! \brief obtains the incident base face */
  Face * face() { return (Face*)(H::face()); }
  
  /*! \brief obtains the incident base face */
  const Face * face() const {return (const Face*)(H::face()); }
  
  /*! \brief obtains the next halfedge along the face */
  const Halfedge * next() const { return (const Halfedge*)(H::next()); }

  /*! \brief obtains the next halfedge along the face. */
  Halfedge * next() { return (Halfedge*)(H::next()); }
  
  //in the future will probably be implemented in a max base
  // const Halfedge* prev() const {return (const Halfedge*)(H::prev()); }

  /*! \brief obtains the opposite halfedge */
  const Halfedge * opposite() const
  {
    return (const Halfedge*)(H::opposite());
  }

  /*! \brief sets the base vertex */
  void set_vertex(Vertex * ve) { H::set_vertex(ve); }

  /*! \brief sets the next halfedge along the face */
  void set_next(Halfedge * h) { H::set_next(h); }

  /*! \brief sets the incident base face */
  void set_face(Face * face) { H::set_face(face); }
  
//private:
  /*! \brief sets the opposite halfedge */
  void set_opposite(void * h) { H::set_opposite(h); }

#if 0
  // Apparently the implementation of the In_place_list requires a copy
  // constructor of the value type!
protected:
  //forbid copy constructor and assignment (only derived classes can use them)

  /*! Copy Constructor */
  _Pm_Halfedge(const _Pm_Halfedge &) {}

  /*! Assignment Operator */
  _Pm_Halfedge & operator=(const _Pm_Halfedge &) { return *this; }
#endif
};

/*! The face class */
template <class V, class H, class F>
class _Pm_Face : public F,
                 public In_place_list_base< _Pm_Face<V,H,F> >
{
public:
  typedef F                     Base;
  typedef _Pm_Vertex<V,H,F>     Vertex;
  typedef _Pm_Halfedge<V,H,F>   Halfedge;
  typedef _Pm_Face<V,H,F>       Face;

  /*! Parameterless Constructor */
  _Pm_Face() {}  
  
  /*! \brief obtains an incident halfedge */
  Halfedge * halfedge() { return (Halfedge*)(F::halfedge()); }

  /*! \brief obtains an incident halfedge */
  const Halfedge * halfedge() const
  {
    return (const Halfedge*)(F::halfedge());
  }

  /*! \brief sets an incident halfedge */
  void set_halfedge(Halfedge * h) { F::set_halfedge(h); }

  typedef I_HalfedgeDS_iterator< typename F::Holes_iterator, 
    Halfedge*, 
    typename F::Holes_iterator::difference_type,
    typename F::Holes_iterator::iterator_category>       Holes_iterator;

  typedef I_HalfedgeDS_const_iterator<
    typename F::Holes_const_iterator, 
    typename F::Holes_iterator,
    const Halfedge*,
    typename F::Holes_const_iterator::difference_type,
    typename F::Holes_const_iterator::iterator_category> Holes_const_iterator;

  /*! \brief adds a hole */
  void add_hole(Halfedge * h) { F::add_hole(h); }

  /*! \brief erases a hole */
  void erase_hole(Holes_iterator hit)
  {
    F::erase_hole(hit.current_iterator());
  }

  /*! \brief erases a range of holes */
  void erase_holes(Holes_iterator first, Holes_iterator last) 
  {
    F::erase_holes(first.current_iterator(), last.current_iterator());
  }

  Holes_iterator holes_begin() {return F::holes_begin(); }
  Holes_iterator holes_end() {return F::holes_end(); }

  Holes_const_iterator holes_begin() const {return F::holes_begin(); }
  Holes_const_iterator holes_end() const {return F::holes_end(); }
  
#if 0
  // Apparently the implementation of the In_place_list requires a copy
  // constructor of the value type!
protected:
  //forbid copy constructor and assignment (only derived classes can use them)

  /*! Copy Constructor */
  _Pm_Face(const _Pm_Face &) {}

  /*! Assignment Operator */
  _Pm_Face & operator=(const _Pm_Face &) { return * this; }
#endif
};

/*! A Dcel Class Using Lists */
template < class V, class H, class F ,
           class Allocator = CGAL_ALLOCATOR(int) >
class Pm_dcel {
public:
  typedef Pm_dcel<V,H,F>        Self;
  typedef _Pm_Vertex<V,H,F>     Vertex;
  typedef _Pm_Halfedge<V,H,F>   Halfedge;
  typedef _Pm_Face<V,H,F>       Face;
  
protected:
  // Three managed in-place lists for the elements. 
  typedef In_place_list<Vertex,false>   Vertex_list;
  typedef In_place_list<Halfedge,false> Halfedge_list;
  typedef In_place_list<Face,false>     Face_list;

  // Vertex allocator
  typedef typename Allocator::template rebind<Vertex>    Vertex_alloc_rebind;
  typedef typename Vertex_alloc_rebind::other            Vertex_allocator;

  // Halfedge allocator
  typedef typename Allocator::template rebind<Halfedge>  Halfedge_alloc_rebind;
  typedef typename Halfedge_alloc_rebind::other          Halfedge_allocator;

  // Face allocator
  typedef typename Allocator::template rebind<Face>      Face_alloc_rebind;
  typedef typename Face_alloc_rebind::other              Face_allocator;




public:
  typedef typename Halfedge_list::size_type             Size;
  typedef typename Halfedge_list::size_type             size_type;
  typedef typename Halfedge_list::difference_type       difference_type;
  typedef typename Halfedge_list::difference_type       Difference;
  typedef std::bidirectional_iterator_tag               iterator_category;

protected:
  Vertex_list vertices;
  Halfedge_list halfedges;
  Face_list faces;
  
  // three allocators (Vertex , Halfedge , Face)
  Vertex_allocator   m_vertex_allocator;
  Halfedge_allocator m_halfedge_allocator;
  Face_allocator     m_face_allocator;
 


public:
  typedef typename Vertex_list::iterator                Vertex_iterator;
  typedef typename Halfedge_list::iterator              Halfedge_iterator;
  typedef typename Face_list::iterator                  Face_iterator;
  typedef CGAL::N_step_adaptor<Halfedge_iterator, 2>    Edge_iterator;
  
  typedef typename Vertex_list::const_iterator          Vertex_const_iterator;
  typedef typename Halfedge_list::const_iterator
    Halfedge_const_iterator;
  typedef typename Face_list::const_iterator            Face_const_iterator;
  typedef CGAL::N_step_adaptor<Halfedge_const_iterator,2>
    Edge_const_iterator;

  /*! Parameterless Constructor */
  Pm_dcel() {}

private:
  // Forbid copy constructor and assignment (will be implemented later).

  /*! Copy Constructor */
  Pm_dcel(const Self & dcel)
  {
    // vertices = dcel.vertices;
    // halfedges = dcel.halfedges;
    // faces = dcel.faces;
  }

  /*! Assignment operator */
  Self & operator=(const Self & dcel)
  {
    // vertices = dcel.vertices;
    // halfedges = dcel.halfedges;
    // faces = dcel.faces;
    return *this;
  }
  
public:
  /*! Destructor */
  ~Pm_dcel() { delete_all(); }

public:
  Size size_of_vertices() const  { return vertices.size(); }
  Size size_of_halfedges() const { return halfedges.size(); }
  Size size_of_faces() const     { return faces.size(); }

  Vertex_iterator   vertices_begin()  { return vertices.begin(); }
  Vertex_iterator   vertices_end()    { return vertices.end(); }
  Halfedge_iterator halfedges_begin() { return halfedges.begin();}
  Halfedge_iterator halfedges_end()   { return halfedges.end(); }
  Face_iterator     faces_begin()     { return faces.begin(); }
  Face_iterator     faces_end()       { return faces.end(); }
  Edge_iterator     edges_begin()     { return halfedges.begin(); }  
  Edge_iterator     edges_end()       { return halfedges.end(); }  

  // The constant iterators and circulators.

  Vertex_const_iterator   vertices_begin() const { return vertices.begin(); }
  Vertex_const_iterator   vertices_end() const { return vertices.end(); }
  Halfedge_const_iterator halfedges_begin() const { return halfedges.begin(); }
  Halfedge_const_iterator halfedges_end() const { return halfedges.end(); }
  Face_const_iterator     faces_begin() const { return faces.begin(); }
  Face_const_iterator     faces_end() const { return faces.end(); }
  Edge_const_iterator     edges_begin() const { return halfedges.begin(); }  
  Edge_const_iterator     edges_end() const { return halfedges.end(); }  
  
  // Insertion
  //
  // The following operations just allocate a new element of that type.
  // Halfedges are always allocated in pairs of opposite halfedges. The
  // opposite pointers are automatically set.
  
  Vertex* new_vertex() {
    //Vertex* v = new Vertex;
    Vertex *v = m_vertex_allocator.allocate(1);
    m_vertex_allocator.construct(v, Vertex());
    vertices.push_back( *v);
    return v;
  }
  
  Vertex* new_vertex( const Vertex* w) {
    //Vertex* v = new Vertex(*w);
    Vertex *v = m_vertex_allocator.allocate(1);
    m_vertex_allocator.construct(v, *w);
    vertices.push_back( *v);
    return v;
  }
  
  /*
  Vertex* new_vertex( const Point& p) {
    Vertex* v = new Vertex(p);
    vertices.push_back( *v);
    return v;
  }
  */

  Halfedge * new_halfedge() {
    //Halfedge * h = new Halfedge;
    Halfedge * h = m_halfedge_allocator.allocate(1);
    m_halfedge_allocator.construct(h, Halfedge());
    halfedges.push_back(*h);
    return h;
  }

  Halfedge * new_halfedge(const Halfedge * he) {
    //Halfedge * h = new Halfedge( *he);
    Halfedge * h = m_halfedge_allocator.allocate(1);
    m_halfedge_allocator.construct(h, *he);
    halfedges.push_back(*h);
    return h;
  }

  Halfedge * new_edge() {
    // creates a new pair of opposite halfedges.
    Halfedge * h = new_halfedge();
    Halfedge * g = new_halfedge();
    h->H::set_opposite(g);
    g->H::set_opposite(h);
    return h;
  }

  Halfedge * new_edge( const Halfedge * he) {
    Halfedge * h = new_halfedge(he);
    Halfedge * g = new_halfedge(he->opposite());
    h->H::set_opposite(g);
    g->H::set_opposite(h);
    return h;
  }
  
  /*
    Halfedge* new_edge(const Curve& c) {
    Halfedge* h = new Halfedge(c);
    Halfedge* g = new Halfedge(c);  //maybe change to flip??
    //    h->H::set_twin(g);
    //g->H::set_twin(h);
    h->H::set_opposite(g);
    g->H::set_opposite(h);

    halfedges.push_back( *h);
    halfedges.push_back( *g);
    return h;
  }*/

  Face * new_face() {
    //Face * f = new Face;
    Face * f = m_face_allocator.allocate(1);
    m_face_allocator.construct(f, Face());
    faces.push_back(*f);
    return f;
  }
  
  Face * new_face(const Face * g) {
    //Face * f = new Face(*g);
    Face * f = m_face_allocator.allocate(1);
    m_face_allocator.construct(f, *g);
    faces.push_back(*f);
    return f;
  }
  
  // Removal
  //
  // The following operations erase an element referenced by a pointer.
  // Halfedges are always deallocated in pairs of opposite halfedges. Erase
  // of single elements is optional. The deletion of all is mandatory.

  void delete_vertex(Vertex * v) {
    vertices.erase(v);
    m_vertex_allocator.destroy(v);
    m_vertex_allocator.deallocate(v,1);
    //delete v;
  }
  
  void delete_halfedge(Halfedge * h) {
    halfedges.erase(h);
    m_halfedge_allocator.destroy(h);
    m_halfedge_allocator.deallocate(h,1);
    //delete h;
  }
    
  void delete_edge(Halfedge * h) {
    // deletes the pair of opposite halfedges h.
    Halfedge * g = h->opposite();
    delete_halfedge(h);
    delete_halfedge(g);
  }

  void delete_face(Face * f) {
    faces.erase(f);
    m_face_allocator.destroy(f);
    m_face_allocator.deallocate(f,1);
    //delete f;
  }
  
  void delete_all() {
    for (Vertex_iterator vi = vertices.begin(); vi != vertices.end();) {
      Vertex_iterator curr = vi;
      ++vi;
      delete_vertex(&(*curr));
    }

    for (Halfedge_iterator hi = halfedges.begin(); hi != halfedges.end();) {
      Halfedge_iterator curr = hi;
      ++hi;
      delete_halfedge(&(*curr));
    }

    for (Face_iterator fi = faces.begin(); fi != faces.end();) {
      Face_iterator curr = fi;
      ++fi;
      delete_face(&(*curr));
    }

    // vertices.destroy();
    // halfedges.destroy();
    // faces.destroy();
  }
  
  /*! returns the unbounded face in the assigned map */
  void * assign(const Self & d, void * u_face)
  {
    //typedef std::map<Vertex_list::iterator, Vertex_list::iterator> VertexMap;
    //typedef std::map<Halfedge_list::iterator, 
    //                 Halfedge_list::iterator> HalfedgeMap;
    //typedef std::map<Face_list::iterator, Face_list::iterator> FaceMap;
    typedef std::map<void*, void*> ConnectMap;

    delete_all();

    ConnectMap vm, hm, fm;
    //VertexMap vm;
    //HalfedgeMap hm;
    //FaceMap fm;

    Vertex_const_iterator vit;
    Halfedge_const_iterator hit;
    Face_const_iterator fit;

    for (vit = d.vertices_begin(); vit != d.vertices_end(); vit++)
    {
      Vertex * nv = new_vertex();
      nv->assign(*vit);
      vm.insert(ConnectMap::value_type((void*)&(*vit), (void*)nv));
    }

    for (hit = d.halfedges_begin(); hit != d.halfedges_end(); hit++)
    {
      Halfedge * nh = new_halfedge();
      nh->assign(*hit);
      hm.insert(ConnectMap::value_type((void*)(&(*hit)), (void*)nh));
    }

    for (fit = d.faces_begin(); fit != d.faces_end(); fit++)
    {
      Face * nf = new_face();
      nf->assign(*fit);
      fm.insert(ConnectMap::value_type((void*)&(*fit), (void*)nf));
    }

    // update pointers
    for (vit = d.vertices_begin(); vit != d.vertices_end(); vit++)
    {
      void *he, *nhe, *nv, *v;
      v = (void*)(&(*vit));
      nv = (void*)(vm.find(v)->second);
      he = (void*)vit->halfedge();
      nhe = (void*)(hm.find(he)->second);
      ((Vertex*)nv)->set_halfedge((Halfedge*)nhe);
    }

    for (hit = d.halfedges_begin(); hit != d.halfedges_end(); hit++)
    {
      void *he, *nhe, *v, *nv, *f, *nf, *op, *nop, *xt, *nxt;
      he = (void*)(&(*hit));
      nhe = hm.find(he)->second;
      v = (void*)hit->vertex();
      f = (void*)hit->face();
      op = (void*)hit->opposite();
      xt = (void*)hit->next();

      nv = vm.find(v)->second;
      nf = fm.find(f)->second;
      nop = hm.find(op)->second;
      nxt = hm.find(xt)->second;

      ((Halfedge*)nhe)->set_vertex((Vertex*)nv);
      ((Halfedge*)nhe)->set_face((Face*)nf);
      ((Halfedge*)nhe)->set_opposite((Halfedge*)nop);
      ((Halfedge*)nhe)->set_next((Halfedge*)nxt);
    }

    for (fit = d.faces_begin(); fit != d.faces_end(); fit++)
    {
      void *f, *nf, *he, *nhe, *h, *nh;
      typename Face::Holes_const_iterator holes;
      f = (void*)(&(*fit));
      nf = fm.find(f)->second;
      he = (void*)fit->halfedge();
      if (he != NULL)
        nhe = hm.find(he)->second;
      else
        nhe = NULL;
      ((Face*)nf)->set_halfedge((Halfedge*)nhe);

      for (holes = fit->holes_begin(); holes != fit->holes_end(); holes++)
      {
        h = (void*)(*holes);
        nh = hm.find(h)->second;
        ((Face*)nf)->add_hole((Halfedge*)nh);
      }
    }
    return fm.find(u_face)->second;
  }
};

/*!
 * DEFAULT DCEL
 */
template <class Traits>
class Pm_default_dcel :
  public Pm_dcel<Pm_vertex_base<typename Traits::Point>,
                 Pm_halfedge_base<typename Traits::X_curve>,
                 Pm_face_base> 
{
public:
  Pm_default_dcel() {}
};

CGAL_END_NAMESPACE

#endif 
