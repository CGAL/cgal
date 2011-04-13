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
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Pm_default_dcel.h
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_DEFAULT_DCEL_H
#define CGAL_PM_DEFAULT_DCEL_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#include <list>
#include <map>
#ifndef CGAL_IN_PLACE_LIST_H
#include <CGAL/In_place_list.h>
#endif
//needed for holes iterator
#ifndef CGAL_PLANAR_MAP_MISC_H
#include <CGAL/Planar_map_2/Planar_map_misc.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Pt>
class Pm_vertex_base {
protected:
  void* hdg;
  Pt    pt;
public:
  typedef Pt           Point;
  
  Pm_vertex_base() {}
  Pm_vertex_base( const Pt& p) : pt(p) {}
  
  virtual ~Pm_vertex_base() {}

  void*       halfedge()               {return hdg;}
  const void* halfedge() const         {return hdg;}
  void        set_halfedge( void* h)   { hdg = h;}
  // an incident halfedge pointing at `v'.
  
  Point&       point()       { return pt;}
  const Point& point() const { return pt;}
  
  void set_point(const Point& p) {pt=p;}

  // assign function for non-connectivity data
  virtual void assign(const Pm_vertex_base<Pt> &v)
  {
	pt = v.pt;
  }
};



template <class X_curve>
class Pm_halfedge_base {
public:
  typedef X_curve Curve;

  Pm_halfedge_base() {}

  Pm_halfedge_base(const X_curve& c) : cv(c) {}

  virtual ~Pm_halfedge_base() {}

  void*       opposite()       { return opp;}
  const void* opposite() const { return opp;}

  void*       next()           { return nxt;}
  const void* next() const     { return nxt;}
  // the next halfedge along the face.
  
  void  set_opposite( void* h)  { opp = h;}

  void  set_next( void* h)      { nxt = h;}
  
  void*       vertex()       { return v;}
  const void* vertex() const { return v;}
  
  void*       face()       { return f;}
  const void* face() const { return f;}
  // the face to the left.

  void  set_vertex( void* _v)     { v = _v;}

  void  set_face( void* _f)      { f = _f;}


  Curve&       curve()       { return cv;}
  const Curve& curve() const { return cv;}

  void set_curve(const X_curve& c) {cv=c;}

  // assign function for non-connectivity data
  virtual void assign(const Pm_halfedge_base<X_curve> &e)
  {
	cv = e.cv;
  }

protected:
  void* opp;
  void* nxt;
  
  void* v; //target
  void* f; //face
  
  X_curve cv;

};


class Pm_face_base {
public:
  typedef std::list<void*> Holes_container; 
  typedef Holes_container::iterator Holes_iterator; 
  typedef Holes_container::const_iterator Holes_const_iterator;


  Pm_face_base() : holes() {};
  virtual ~Pm_face_base() {}

  void* halfedge() { return hdg;}
  const void* halfedge() const { return hdg;}

  void set_halfedge(void* h) {hdg=h;}

  Holes_iterator  holes_begin() {return holes.begin();}
  Holes_iterator  holes_end() {return holes.end();}

  Holes_const_iterator  holes_begin() const {return holes.begin();}
  Holes_const_iterator  holes_end() const {return holes.end();}


  void add_hole(void* halfedge_ptr)
  {

    holes.push_back(halfedge_ptr);

  }


  void erase_hole(Holes_iterator hit) {
    holes.erase(hit);
  }
  void erase_holes(Holes_iterator first, Holes_iterator last) {
    holes.erase(first,last);
  }

  // assign function for non-connectivity data
  virtual void assign(const Pm_face_base &f)
  {
  }
  

private:
  void* hdg;
  Holes_container holes ;

};





template <class V, class H, class F> class _Pm_Vertex;
template <class V, class H, class F> class _Pm_Halfedge;
template <class V, class H, class F> class _Pm_Face;

template <class V, class H, class F>
class _Pm_Vertex
    : public  V,
      public  In_place_list_base< _Pm_Vertex<V,H,F> >
{
public:
  typedef V                              Base;
  typedef _Pm_Vertex<V,H,F>    Vertex;
  typedef _Pm_Halfedge<V,H,F>  Halfedge;
  typedef _Pm_Face<V,H,F>     Face;

  _Pm_Vertex() {}
  //    _Pm_Vertex( const Point& p) : V(p) {}

  Halfedge*       halfedge()       {return (Halfedge*)(V::halfedge());}
  const Halfedge* halfedge() const {return (const Halfedge*)(V::halfedge());}
  void            set_halfedge( Halfedge* h) { V::set_halfedge(h);}

  /* irrelevant
#if _MSC_VER>=1100
public:
#else
protected:
#endif
  //forbid copy constructor and assignment (only derived classes can use them)
  _Pm_Vertex( const _Pm_Vertex&);
  _Pm_Vertex& operator=(const _Pm_Vertex&);
  */

};

template <class V, class H, class F>
class _Pm_Halfedge
  : public  H,
    public  In_place_list_base< _Pm_Halfedge<V,H,F> >
{
public:
  typedef H                              Base;
  typedef _Pm_Vertex<V,H,F>    Vertex;
  typedef _Pm_Halfedge<V,H,F>  Halfedge;
  typedef _Pm_Face<V,H,F>     Face;

  _Pm_Halfedge() : H() {}  
  //_Pm_Halfedge( const Curve& c) : H(c) {}

  Halfedge*       opposite()       {return (Halfedge*)(H::opposite());}
  
  Halfedge*       next()           {return (Halfedge*)(H::next());}
  //in the future will probably be implemented in a max base
  //    const Halfedge* prev()     const {return (const Halfedge*)(H::prev());}

  Vertex*         vertex()         {return (Vertex*)(H::vertex());}

  Face*          face()          {return (Face*)(H::face());}
  
  const Halfedge* opposite() const {return (const Halfedge*)(H::opposite());}

  const Halfedge* next()     const {return (const Halfedge*)(H::next());}
  //in the future will probably be implemented in a max base
  //    const Halfedge* prev()     const {return (const Halfedge*)(H::prev());}

  const Vertex*   vertex()   const {return (const Vertex*)(H::vertex());}

  const Face*    face()    const {return (const Face*)(H::face());}
  
  void  set_next( Halfedge* h)     { H::set_next(h);}

  void  set_vertex( Vertex* ve)    { H::set_vertex(ve);}

  void  set_face( Face* face)   { H::set_face(face);}
  
//private:
  void  set_opposite( void* h)     { H::set_opposite(h);}

  /*
#if _MSC_VER>=1100
public:
#else
protected:
#endif
  //forbid copy constructor and assignment (only derived classes can use them)
  
  _Pm_Halfedge( const _Pm_Halfedge&);
  _Pm_Halfedge& operator=(const _Pm_Halfedge&);
  */  

};


template <class V, class H, class F>
class _Pm_Face
    : public  F,
      public  In_place_list_base< _Pm_Face<V,H,F> >
{
public:
  typedef F                              Base;
  typedef _Pm_Vertex<V,H,F>    Vertex;
  typedef _Pm_Halfedge<V,H,F>  Halfedge;
  typedef _Pm_Face<V,H,F>     Face;

  _Pm_Face() {}  
  
  Halfedge*       halfedge()       {return (Halfedge*)(F::halfedge());}
  const Halfedge* halfedge() const {return (const Halfedge*)(F::halfedge());}
  void            set_halfedge( Halfedge* h) { F::set_halfedge(h);}


  typedef _Polyhedron_iterator< typename F::Holes_iterator, 
    Halfedge*, 
    typename F::Holes_iterator::difference_type,
    typename F::Holes_iterator::iterator_category>       Holes_iterator;

  typedef _Polyhedron_const_iterator<
  typename F::Holes_const_iterator, 
  typename F::Holes_iterator,
  const Halfedge*,
  typename F::Holes_const_iterator::difference_type,
  typename F::Holes_const_iterator::iterator_category>   Holes_const_iterator;

  void add_hole(Halfedge* h) { F::add_hole(h); }
  void erase_hole(Holes_iterator hit) {F::erase_hole(hit.current_iterator());}
  void erase_holes(Holes_iterator first, Holes_iterator last) 
  {F::erase_holes(first.current_iterator(), last.current_iterator());}


  Holes_iterator holes_begin() {return F::holes_begin();}
  Holes_iterator holes_end() {return F::holes_end();}

  Holes_const_iterator holes_begin() const {return F::holes_begin();}
  Holes_const_iterator holes_end() const {return F::holes_end();}
  
  /*
#if _MSC_VER>=1100
public:
#else
protected:
#endif
  //forbid copy constructor and assignment (only derived classes can use them)
  _Pm_Face( const _Pm_Face&);
  _Pm_Face& operator=(const _Pm_Face&);
  */

};


// A Dcel Class Using Lists
// ----------------------------------
//

template < class V, class H, class F>
class Pm_dcel {
public:
  typedef Pm_dcel<V,H,F>   Self;
  typedef _Pm_Vertex<V,H,F>         Vertex;
  typedef _Pm_Halfedge<V,H,F>       Halfedge;
  typedef _Pm_Face<V,H,F>          Face;
  
protected:
  // Three managed in-place lists for the elements. 
  typedef In_place_list<Vertex,true>      Vertex_list;
  typedef In_place_list<Halfedge,true>    Halfedge_list;
  typedef In_place_list<Face,true>       Face_list;
public:
  typedef typename Halfedge_list::size_type   Size;
  typedef typename Halfedge_list::size_type   size_type;
  typedef typename Halfedge_list::difference_type difference_type;
  typedef typename Halfedge_list::difference_type Difference;
  typedef std::bidirectional_iterator_tag          iterator_category;

protected:
    Vertex_list    vertices;
    Halfedge_list  halfedges;
    Face_list     faces;

public:


  typedef typename Vertex_list::iterator      Vertex_iterator;
  typedef typename Halfedge_list::iterator    Halfedge_iterator;
  typedef typename Face_list::iterator       Face_iterator;
  
  typedef typename Vertex_list::const_iterator    Vertex_const_iterator;
  typedef typename Halfedge_list::const_iterator  Halfedge_const_iterator;
  typedef typename Face_list::const_iterator     Face_const_iterator;
  
  // CREATION

  Pm_dcel() {}

#if _MSC_VER>=1100
public:
#else
private:
#endif
  // Forbid copy constructor and assignment (will be implemented later).
  Pm_dcel( const Self&) {}
  Self& operator=( const Self&)            { return *this;}
  
public:
  // Access Member Functions
  
  Size size_of_vertices() const  { return vertices.size();}
  Size size_of_halfedges() const { return halfedges.size();}
  Size size_of_faces() const    { return faces.size();}

  Vertex_iterator   vertices_begin()   { return vertices.begin();}
  Vertex_iterator   vertices_end()     { return vertices.end();}
  Halfedge_iterator halfedges_begin()  { return halfedges.begin();}
  Halfedge_iterator halfedges_end()    { return halfedges.end();}
  Face_iterator    faces_begin()     { return faces.begin();}
  Face_iterator    faces_end()       { return faces.end();}
  
  // The constant iterators and circulators.

  Vertex_const_iterator   vertices_begin()  const{ return vertices.begin();}
  Vertex_const_iterator   vertices_end()    const{ return vertices.end();}
  Halfedge_const_iterator halfedges_begin() const{ return halfedges.begin();}
  Halfedge_const_iterator halfedges_end()   const{ return halfedges.end();}
  Face_const_iterator    faces_begin()    const{ return faces.begin();}
  Face_const_iterator    faces_end()      const{ return faces.end();}
  
  // Insertion
  //
  // The following operations just allocate a new element of that type.
  // Halfedges are always allocated in pairs of opposite halfedges. The
  // opposite pointers are automatically set.
  
  Vertex* new_vertex() {
    Vertex* v = new Vertex;
    vertices.push_back( *v);
    return v;
  }
  
  Vertex* new_vertex( const Vertex* w) {
    Vertex* v = new Vertex(*w);
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

  Halfedge* new_edge() {
    // creates a new pair of opposite halfedges.
    Halfedge* h = new Halfedge;
    Halfedge* g = new Halfedge;
    h->H::set_opposite(g);
    g->H::set_opposite(h);

    halfedges.push_back( *h);
    halfedges.push_back( *g);
    return h;
  }
  
  Halfedge* new_edge( const Halfedge* he) {
    Halfedge* h = new Halfedge( *he);
    Halfedge* g = new Halfedge(* (he->opposite()));
    h->H::set_opposite(g);
    g->H::set_opposite(h);

    halfedges.push_back( *h);
    halfedges.push_back( *g);
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
  
  Face* new_face(){
    Face* f = new Face;
    faces.push_back( *f);
    return f;
  }
  
  Face* new_face( const Face* g) {
    Face* f = new Face(*g);
    faces.push_back( *f);
    return f;
  }
  
  // Removal
  //
  // The following operations erase an element referenced by a pointer.
  // Halfedges are always deallocated in pairs of opposite halfedges. Erase
  // of single elements is optional. The deletion of all is mandatory.

  void delete_vertex( Vertex* v) { vertices.erase(v);}
  
  void delete_edge( Halfedge* h) {
    // deletes the pair of opposite halfedges h.
    Halfedge* g = h->opposite();
    halfedges.erase(h);
    halfedges.erase(g);
  }

  void delete_face( Face* f) { faces.erase(f);}
  
  void delete_all() {
    vertices.destroy();
    halfedges.destroy();
    faces.destroy();
  }
  
  // returns the unbounded face in the assigned map
  void *assign(const Self &d, void *u_face)
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
	Vertex* nv = new Vertex;
	nv->assign(*vit);
	vertices.push_back(*nv);
	vm.insert(ConnectMap::value_type((void*)&(*vit), (void*)nv));
      }

    for (hit = d.halfedges_begin(); hit != d.halfedges_end(); hit++)
      {
	Halfedge* nh = new Halfedge;
	nh->assign(*hit);
	halfedges.push_back(*nh);
	hm.insert(ConnectMap::value_type((void*)(&(*hit)), (void*)nh));
      }

    for (fit = d.faces_begin(); fit != d.faces_end(); fit++)
      {
	Face* nf = new Face;
	nf->assign(*fit);
	faces.push_back(*nf);
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



///////////////////////////////////////////////////////////////
//               DEFAULT DCEL
///////////////////////////////////////////////////////////////

template <class Traits>
class Pm_default_dcel
  : public Pm_dcel<
Pm_vertex_base<typename Traits::Point>,
Pm_halfedge_base<typename Traits::X_curve>,
Pm_face_base
> 
{
public:  // CREATION
  
  Pm_default_dcel() {}
  
};


CGAL_END_NAMESPACE

#endif 
// EOF //









