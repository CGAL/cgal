// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Arr_2_bases.h
// package       : arr (1.03)
// author(s)     : Iddo Hanniel 
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_2_BASES_H
#define CGAL_ARR_2_BASES_H

#include <list>

CGAL_BEGIN_NAMESPACE

template <class Pt>
class Arr_2_vertex_base {
protected:
  void* hdg;
  Pt    pt;
public:
  typedef Pt           Point;

  Arr_2_vertex_base() {}

  virtual ~Arr_2_vertex_base() {}
  
  void*       halfedge()               { return hdg; }
  const void* halfedge() const         { return hdg; }
  void        set_halfedge(void* h)    { hdg = h; }
  // an incident halfedge pointing at `v'.

  //    Point&       point()       { return pt;}
  const Point& point() const { return pt;}
  void set_point(const Point& p) { pt = p; }

  // assign function for non-connectivity data
  virtual void assign(const Arr_2_vertex_base<Pt> & v) { pt = v.pt; }
};

template <class Base_node>
class Arr_2_halfedge_base {
public:
  typedef typename Base_node::X_monotone_curve_2 X_monotone_curve_2;

  Arr_2_halfedge_base() : bn(0) {}
  
  virtual ~Arr_2_halfedge_base() {}
  
  void*       opposite()       { return opp; }
  const void* opposite() const { return opp; }

  void*       next()           { return nxt; }
  const void* next() const     { return nxt; }
  // the next halfedge along the face.
  
  void  set_opposite(void* h)  { opp = h; }

  void  set_next(void* h)      { nxt = h; }
  

  //    void*       prev()       { return prv;}
  //  const void* prev() const { return prv;}


  void*       vertex()       { return v; }
  const void* vertex() const { return v; }
  
  void*       face()       { return f; }
  const void* face() const { return f; }
  // the face to the left.

  void  set_vertex( void* _v)     { v = _v; }

  void  set_face( void* _f)      { f = _f; }


  //for debug only !!
  //const Curve&       curve() const  { return cv;}
  //void set_curve(const Curve& c) {cv=c;}

  //WATCH OUT:
  //we make the curve and set_curve empty so the pm can find them but not
  //use them , the curves are set in the halfedge via the edge_node!!

  
  const X_monotone_curve_2 & curve() const 
  { 
    return bn->x_curve(); 
  }
//  void set_curve(const Curve& cv) {bn->set_curve(cv);}
//the setting of the curve is done only in the arrangement level
  void set_curve(const X_monotone_curve_2 &) {}
    
  Base_node* edge_node() {return bn;} //will become private in the arrangement
  const Base_node* edge_node() const {return bn;} 
  void set_edge_node(Base_node* b) {bn=b;}

  // assign function for non-connectivity data
  virtual void assign(const Arr_2_halfedge_base<Base_node> &)
  {
    //bn = new Base_node(*e.bn);
  }

protected:

  void* opp;
  void* nxt;
  
  void* v; 
  void* f; //face
  
  Base_node* bn;

  //debug
  //Curve cv;

};


class Arr_2_face_base {
public:
  typedef std::list<void*> Holes_container; 
  
  typedef Holes_container::iterator Holes_iterator; 
  typedef Holes_container::const_iterator Holes_const_iterator;


  Arr_2_face_base() : holes() {};
  
  virtual ~Arr_2_face_base() {}
  
  void* halfedge() { return hdg;}
  const void* halfedge() const { return hdg;}

  void set_halfedge(void* h) {hdg=h;}

  //mine

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

  
  //this is not documented but needed for a private project
  Holes_container::size_type number_of_holes() const { return holes.size(); }

  // assign function for non-connectivity data
  virtual void assign(const  Arr_2_face_base &) { }

protected:
  void* hdg;
  Holes_container holes;
};



template <class _Curve_2, class _X_monotone_curve_2>
class Arr_base_node {
public:
  typedef _Curve_2 Curve_2;
  typedef _X_monotone_curve_2 X_monotone_curve_2;

  class Curve_wrap{
  public:
    Curve_2* cv;
    X_monotone_curve_2* x_cv;

    Curve_wrap ():cv(NULL), x_cv(NULL)
    {}

    ~Curve_wrap()
    {
      if(cv) delete cv;
      if(x_cv) delete x_cv;
    }

    Curve_wrap(const Curve_wrap& cv_wrap)
    {
      if(cv_wrap.cv!=NULL)
      {
	cv = new Curve_2;
	*cv = *(cv_wrap.cv);
      }
      else
	cv=NULL;

      if(cv_wrap.x_cv!=NULL)
      {
	x_cv = new X_monotone_curve_2;
	*x_cv = *(cv_wrap.x_cv);
      }
      else
	x_cv=NULL;
    }

    Curve_wrap& operator= (const Curve_wrap& cv_wrap)
    {
      if (this == &cv_wrap)
	return (*this);

      if(cv_wrap.cv!=NULL)
      {
	cv = new Curve_2;
	*cv = *(cv_wrap.cv);
      }
      else
	cv=NULL;

      if(cv_wrap.x_cv!=NULL)
      {
	x_cv = new X_monotone_curve_2;
	*x_cv = *(cv_wrap.x_cv);
      }
      else
	x_cv=NULL;

      return *this;
    }
  };

  Arr_base_node() {}
  virtual ~Arr_base_node() {}

  const Curve_2& curve() const {return *(cv_wrap.cv);}
  void set_curve(const Curve_2& c) 
  {
    *cv_wrap.cv = c;
  }

  const X_monotone_curve_2& x_curve() const 
  {
    return *(cv_wrap.x_cv);
  }
  void set_x_monotone_curve(const X_monotone_curve_2& c) {*cv_wrap.x_cv = c;}

  // assign function for non-connectivity data
  virtual void assign(const Arr_base_node<_Curve_2, _X_monotone_curve_2> &bn)
  {
    cv_wrap = bn.cv_wrap;
  }

protected: 
//  Curve cv;
  Curve_wrap cv_wrap;
};

CGAL_END_NAMESPACE

#endif
