// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Mariette Yvinec


#ifndef CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H
#define CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H

// #include <CGAL/circulator.h>
// #include <CGAL/Triangulation_utils_2.h>
// #include <CGAL/triangulation_assertions.h>
// #include <CGAL/Triangulation_face_2.h>
// #include <CGAL/Triangulation_vertex_2.h>
// #include <CGAL/Triangulation_handles_2.h>

namespace CGAL {


template <class Triangulation_ > //  < class Gt, class Tds >
class Triangulation_line_face_circulator_2
  :   public Bidirectional_circulator_base< typename Triangulation_::Triangulation_data_structure::Face,
	                                    std::ptrdiff_t,
                                            std::size_t>,
      public Triangulation_cw_ccw_2
{
public:

  typedef Triangulation_line_face_circulator_2<Triangulation_> Line_face_circulator;
  typedef Triangulation_                      Triangulation;
  typedef typename Triangulation::Geom_traits  Gt;
  typedef typename Triangulation_::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex                 Vertex;
  typedef typename Tds::Face                   Face;
  typedef typename Tds::Edge                   Edge;
  typedef typename Tds::Vertex_handle          Vertex_handle;
  typedef typename Tds::Face_handle            Face_handle;
  typedef typename Tds::Face_circulator        Face_circulator;
  
  typedef typename Gt::Point_2 Point;
  typedef typename Triangulation::Locate_type Locate_type;

   enum State {undefined = -1,
	       vertex_vertex,
	       vertex_edge,
	       edge_vertex,
	       edge_edge};
            
private:
  Face_handle pos;
  const Triangulation* _tr;
  State s;
  int i;
  Point p, q;
            
public:
  Triangulation_line_face_circulator_2()
    : pos(), _tr(NULL), s(undefined), i(-1)
    {}
            
  Triangulation_line_face_circulator_2(Vertex_handle v,
				       const Triangulation* tr,
				       const Point& dir);

  Triangulation_line_face_circulator_2(const Point& pp,
				       const Point& qq,
				       const Triangulation * t);
           
  Triangulation_line_face_circulator_2(const Point& pp,
				       const Point& qq,
				       const Face_handle& ff,
				       const Triangulation* t); 

  Line_face_circulator&   operator++() ;
  Line_face_circulator&   operator--() ;
  Line_face_circulator    operator++(int);
  Line_face_circulator    operator--(int);
  Face*                   operator->() {return &*pos;}
  Face&                   operator*() { return *pos;}
  Face_handle             handle() {return pos;}
  operator const Face_handle() const {return pos;}
  bool  operator==(const Line_face_circulator& lfc) const;
  bool  operator!=(const Line_face_circulator& lfc) const;

  bool  operator==(const Face_handle& fh) const { return fh == pos; }
  bool  operator!=(const Face_handle& fh) const { return fh != pos; }

  bool  operator==(Nullptr_t  CGAL_triangulation_assertion_code(n)) const;
  bool  operator!=(Nullptr_t n) const;
  bool  is_empty() const;
  bool  collinear_outside() const;
  bool locate(const Point& t, Locate_type &lt,  int &li);

  //private:
  Triangulation_line_face_circulator_2(const Face_handle& face,
				       int index,
				       State state,
				       const Triangulation * t,
				       const Point& pp,
				       const Point& qq);
private:
  void increment();
  void decrement();
};

template < class Triangulation >
inline
bool
operator==(typename Triangulation::Triangulation_data_structure::Face_handle fh, 
	   Triangulation_line_face_circulator_2<Triangulation> fc)
{
  return (fc==fh);
}

template < class Triangulation >
inline
bool
operator!=(typename Triangulation::Triangulation_data_structure::Face_handle fh, 
	   Triangulation_line_face_circulator_2<Triangulation> fc)
{
  return (fc!=fh);
}

template < class Triangulation >
Triangulation_line_face_circulator_2<Triangulation>::
Triangulation_line_face_circulator_2(const Face_handle& face,
				     int index,
				     State state,
				     const Triangulation * t,
				     const Point& pp,
				     const Point& qq)
    : pos(face), _tr(t), s(state), i(index),  
      p(pp), q(qq)            
{
  CGAL_triangulation_precondition(! t->xy_equal(p, q));
}


template < class Triangulation >
Triangulation_line_face_circulator_2<Triangulation>::
Triangulation_line_face_circulator_2(Vertex_handle v,
				     const Triangulation* tr,
				     const Point& dir)
  :pos(), _tr(tr), s(undefined)
  // begin at the face incident to v, traversed by the ray from v to
  // dir 
  // or null iterator
{
  CGAL_triangulation_precondition((!_tr->is_infinite(v)) &&
			       (_tr->dimension() == 2)  &&
			       (! _tr->xy_equal(v->point(),dir)));
  p=v->point();
  q=dir;

  // find a finite vertex to the left of pq
  // if there is no, the line_face_circulator is null
  Face_circulator fc = _tr->incident_faces(v);
  Face_circulator done(fc);
  int ic = fc->index(v);
  Vertex_handle  vt= fc->vertex(cw(ic));
   while( _tr->is_infinite(vt) || 
	  _tr->orientation(p, q, vt->point()) != LEFT_TURN) {
    ++fc;
    ic = fc->index(v);
    vt= fc->vertex(cw(ic));
    if (fc == done) { *this = Line_face_circulator(); return;}
  }
  
  // now vt is finite and to the left of pq
  Vertex_handle vr = fc-> vertex(ccw(ic));
  Orientation pqr = RIGHT_TURN; // warning "pqr might be used uninitialized"
  while ( (!_tr->is_infinite(vr)) &&
	  (pqr = _tr->orientation(p, q, vr->point()))== LEFT_TURN ) {
    --fc;
    ic = fc->index(v);
    vr = fc-> vertex(ccw(ic));
  }

  // vr can be infinite or finite. 
  // If finite [pqr] is COLLINEAR or RIGHT_TURN
  // reset vt and conclude.  vt is still finite and [pqt] still LEFT_TURN
  ic = fc->index(v);
  vt= fc->vertex(cw(ic));
  CGAL_triangulation_assertion (_tr->orientation(p,q, vt->point())==
				LEFT_TURN );
  if (_tr->is_infinite(vr)) {
    --fc;
    ic = fc->index(v);
    vr = fc->vertex(ccw(ic));
    pqr = _tr->orientation(p, q, vr->point());
    switch(pqr){
    case RIGHT_TURN:
    case COLLINEAR:
      ++fc;
      ic = fc->index(_tr->infinite_vertex());
      pos = fc;
      s = vertex_vertex;
      i = ic;
      break;
    case LEFT_TURN:
     *this = Line_face_circulator(); 
     break;
    }
  }
  else if (pqr == COLLINEAR) {
    pos = fc;
    s = vertex_vertex;
    i = ccw(ic);
   }
  else { // pqr==RIGHT_TURN 
    pos = fc;
    s = vertex_edge;
    i = ic ;
  }
  return;
}
  



template < class Triangulation >
Triangulation_line_face_circulator_2<Triangulation>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Triangulation * t)
     : pos(), _tr(t), s(undefined), p(pp), q(qq)
  //begins at the  first finite face traversed be the oriented line pq
{
  Vertex_handle inf = _tr->infinite_vertex();
  Face_circulator fc = _tr->incident_faces(inf),
    done(fc);

  i = fc->index(inf);
  Point l = fc->vertex(cw(i))->point();
  Point r = fc->vertex(ccw(i))->point();
  Orientation pql = _tr->orientation(p, q, l);
  Orientation pqr = _tr->orientation(p, q, r);
            
   do{
    if( (pql == LEFT_TURN) && (pqr == RIGHT_TURN) ){
      *this = ++Line_face_circulator( fc, i, vertex_edge, t, p, q);
           return;
    } 
    else if ( (pql == LEFT_TURN) && (pqr == COLLINEAR) ){
      --fc;
      i = fc->index(inf);
      Point  ss = fc->vertex(ccw(i))->point();
      Orientation pqs  = _tr->orientation(p, q, ss);
      Face_handle fn;
      int in;
      switch(pqs) {
      case LEFT_TURN:
	*this = Line_face_circulator();
	return;
      case COLLINEAR:
	fn = fc->neighbor(i);
	in = fn->index(fc);
	*this = Line_face_circulator( fn, cw(in),vertex_vertex,t,p,q);
	return;
      case RIGHT_TURN:
	fn = fc->neighbor(i);
	Vertex_handle vr = fc->vertex(cw(i)); // vertex corresponding to r
	in = fn->index(vr);
	ss = fn->vertex(cw(in))->point();
	pqs = _tr->orientation(p, q, ss);
	Orientation pqss = RIGHT_TURN;
	while ( pqs != LEFT_TURN) {
	  pqss = pqs;
	  fn = fn->neighbor(ccw(in));
	  in = fn->index(vr);
	  ss = fn->vertex(cw(in))->point();
	  pqs = _tr->orientation(p, q, ss);
	}
	if (pqss == RIGHT_TURN)
	  *this = Line_face_circulator( fn, in ,vertex_edge,t,p,q);
	else // pqss = COLLINEAR
	  *this = Line_face_circulator(fn,ccw(in),vertex_vertex,t,p,q);
	return;
      }
    }

    // going CCW around convex hull is CW around infinite vertex
    --fc; 
    l = r;
    pql = pqr;
    i = fc->index(inf);
    r = fc->vertex(ccw(i))->point();
    pqr = _tr->orientation(p, q, r);
  }while(fc != done);

   // if line (p,q) does not intersect the convex hull in an edge
   // the circulator has a singular value
   *this=Line_face_circulator();
   return;
}


template < class Triangulation >
Triangulation_line_face_circulator_2<Triangulation>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Face_handle& ff,
				     const Triangulation* t)
  : pos(ff), _tr(t), s(undefined), p(pp), q(qq)
  // precondition : face ff contain p
  // the walk  begins at face ff if ff is a finite face traversed by the
  // circulator
  // if ff is finite but not traversed by the circulator
  // (this happens when p is a vertex of ff or on an edge) :
  // the circulator may be empty, or the walk begins at a finite face
  // incident to p 
  // if ff is infinite, the walk begin at the first finite face traversed
{
  CGAL_triangulation_precondition(_tr->is_infinite(ff) ||
			      _tr->oriented_side(ff,p) != ON_NEGATIVE_SIDE);
  int j;
  if(_tr->is_infinite(pos)){
    *this  = Line_face_circulator(p, q, t);
    return;
  }

  // Test whether p lies on a vertex
  for(j = 0; j < 3; j++){
    if(_tr->xy_equal(pos->vertex(j)->point(), p)){
      *this = Line_face_circulator( pos->vertex(j), t, q);
      if( (!is_empty()) && _tr->is_infinite(pos )) --(*this);
      return;
    }
  }
            
  // Test whether p lies on an edge
  for(j = 0; j < 3; j++) {
    if(_tr->orientation(pos->vertex(ccw(j))->point(),
			pos->vertex(cw(j))->point(),
			p) == COLLINEAR){
      Orientation pqj =
	_tr->orientation(p, q, pos->vertex(j)->point());
      Orientation pqcwj =
	_tr->orientation(p, q, pos->vertex(cw(j))->point());
      switch(pqcwj) {
      case COLLINEAR :
	if(pqj == LEFT_TURN){
	  s = vertex_vertex;
	  i = cw(j);
	  return;
	} 
	else if(! _tr->is_infinite(pos->neighbor(j))){
	  Face_handle n = pos->neighbor(j);
	  i = cw(n->index(pos)); 
	  pos = n;
	  s = vertex_vertex;
	  return;
	} else {
           // singular value
	  *this = Line_face_circulator();
	  return;
	}
      case LEFT_TURN :
	i = j;
	s = (pqj == COLLINEAR) ? vertex_edge :  
	  edge_edge;
	break;
      case RIGHT_TURN :
	switch(pqj){
	case COLLINEAR:
	  s = edge_vertex;
	  i = j;
	  return;
	case LEFT_TURN:
	  s = edge_edge;
	  i = ccw(j);
	  return;
	case RIGHT_TURN:
	  s = edge_edge;
	  i = cw(j);
	  return;
	}
      }
    }
  }

  // p lies in the interior of the face
  Orientation orient[3];
  for(j=0; j<3; j++) {
    orient[j] =
      _tr->orientation(p,q,pos->vertex(j)->point());
  }
  for(j=0; j<3; j++) {
    if(orient[j] == COLLINEAR) {
      i = j;
      s = (orient[ccw(j)] == LEFT_TURN) ? edge_vertex : 
	vertex_edge;
      return;
    }
  }
  s = edge_edge;
  for(j=0; j<3; j++){
    if(orient[j] == RIGHT_TURN){
      i = (orient[ccw(j)] == RIGHT_TURN) ? j : cw(j);
      return;
    }
  }
}


template < class Triangulation >
inline
void
Triangulation_line_face_circulator_2<Triangulation>::
increment()
{
  CGAL_triangulation_precondition(pos != Face_handle());
  if(s == vertex_vertex || s == edge_vertex) {
    Orientation o;
    do{
      Face_handle n = pos->neighbor(cw(i));
      i = n->index(pos);
      pos = n;
      if (pos->vertex(i) == _tr->infinite_vertex()){
	o = COLLINEAR;
	i = cw(i);
	break;
      }
      o = _tr->orientation(p, q, pos->vertex(i)->point());
      i = cw(i);
    }while(o == LEFT_TURN);
            
    if(o == COLLINEAR) {
      s = vertex_vertex;
      i = ccw(i);
    } 
    else {
      s = vertex_edge;
    }
  } 
  else {
    Face_handle n = pos->neighbor(i);
    int ni = n->index(pos);
    pos = n ;
    Orientation o = _tr->is_infinite(pos->vertex(ni)) ?
      COLLINEAR :
      _tr->orientation(p,q,pos->vertex(ni)->point());
            
    switch(o){
    case LEFT_TURN:
      s = edge_edge;
      i = ccw(ni);
      break;
    case RIGHT_TURN:
      s = edge_edge;
      i = cw(ni);
      break;
    default:
      s = edge_vertex;
      i = ni;
    }
  }
} 
            
template < class Triangulation >
void
Triangulation_line_face_circulator_2<Triangulation>::             
decrement()
{
  CGAL_triangulation_precondition(pos != Face_handle());
  if(s == vertex_vertex || s == vertex_edge) {
    if(s == vertex_vertex){
      i = cw(i);
    }
    Orientation o;
    do{
      Face_handle n = pos->neighbor(ccw(i));
      i = n->index(pos);
      pos = n;
      if (pos->vertex(i) == _tr->infinite_vertex()){
	o = COLLINEAR;
	i = ccw(i);
	break;
      }
      o = _tr->orientation(p, q, pos->vertex(i)->point());
      i = ccw(i);
    }while(o == LEFT_TURN);
            
    s = (o == COLLINEAR) ? vertex_vertex : edge_vertex;
            
  } 
  else { // s == edge_edge  ||  s == edge_vertex
    // the following is not nice. A better solution is to say
    // that index i is at the vertex that is alone on one side of l(p,q)
    if(s == edge_edge){
      i = (_tr->orientation
	   (p, q,
	    pos->vertex(i)->point()) == 
	   LEFT_TURN)
	? cw(i) : ccw(i);
    }
    Face_handle n = pos->neighbor(i);
    i = n->index(pos);
    pos = n;
    Orientation o = _tr->is_infinite(pos->vertex(i)) ?
      COLLINEAR :
      _tr->orientation(p, q, pos->vertex(i)->point());
            
    s = (o == COLLINEAR) ? vertex_edge : edge_edge;
  }
}

template < class Triangulation >
bool
Triangulation_line_face_circulator_2<Triangulation>::
locate(const Point& t, Locate_type &lt,  int &li)
{
  switch(s){            
  case edge_edge:
  case vertex_edge:
    {
      Orientation o =
	_tr->orientation(pos->vertex(ccw(i))->point(),
			 pos->vertex(cw(i))->point(),
			 t);
      if(o == RIGHT_TURN)      return false;
      if(o == COLLINEAR){
	lt = Triangulation::EDGE;
	li = i;
	return true;
      }
      lt = Triangulation::FACE;
      li = 4;//li unused in this case
      return true;
    }
  case vertex_vertex:
    {
      if(_tr->is_infinite(pos->vertex(i))){
	CGAL_triangulation_assertion(
	       _tr->orientation( pos->vertex(cw(i))->point(),
				 pos->vertex(ccw(i))->point(),
				 t) != LEFT_TURN);
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
      const Point &u = pos->vertex(cw(i))->point();
      const Point &v = pos->vertex(i)->point();
      // u == t  was detected earlier
      if(_tr->compare_x(v,t)==EQUAL && 
	 _tr->compare_y(v,t)==EQUAL){
	lt = Triangulation::VERTEX;
	li = i;
	return true;
      }
      if(_tr->collinear_between(u, t, v)) {
	lt = Triangulation::EDGE;
	li = ccw(i);
	return true;
      }
      return false;
    }
  default: // edge_vertex
    {
      if(_tr->is_infinite(pos->vertex(i))){
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
      if(_tr->xy_equal(t,pos->vertex(i) ->point()) ){
	li = i;
	lt = Triangulation::VERTEX;
	return true;
      }
      if(_tr->collinear_between(p, t, pos->vertex(i)->point())) {
	lt = Triangulation::FACE;
	return true;
      }
      return false;
    }
  }
}
           
template < class Triangulation >
inline
Triangulation_line_face_circulator_2<Triangulation>&
Triangulation_line_face_circulator_2<Triangulation>::
operator++()
{
  CGAL_triangulation_precondition( pos != Face_handle()) ;
  increment();
  return *this;
}
            
template < class Triangulation >
inline
Triangulation_line_face_circulator_2<Triangulation>&
Triangulation_line_face_circulator_2<Triangulation>::            
operator--()
{
  CGAL_triangulation_precondition(pos != Face_handle()) ;
  decrement();
  return *this;
}
            
template < class Triangulation >
inline
Triangulation_line_face_circulator_2<Triangulation>
Triangulation_line_face_circulator_2<Triangulation>::             
operator++(int)
{
  Line_face_circulator tmp(*this);
  ++(*this);
  return tmp;
}
            
template < class Triangulation >
inline
Triangulation_line_face_circulator_2<Triangulation>
Triangulation_line_face_circulator_2<Triangulation>::             
operator--(int)
{
  Line_face_circulator tmp(*this);
  --(*this);
  return tmp;
}

template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>::    
operator==(const Line_face_circulator& lfc) const
{
  CGAL_triangulation_precondition( pos != Face_handle() &&
			       lfc.pos != Face_handle());
  return ( pos == lfc.pos &&  _tr == lfc._tr &&
            s== lfc.s && p==lfc.p && q==lfc.q);
}

template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>:: 
operator!=(const Line_face_circulator& lfc) const
{
  return !(*this == lfc);
}
            
template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>::   
is_empty() const
{
  return pos == Face_handle();
}

template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>::            
operator==(Nullptr_t CGAL_triangulation_assertion_code(n)) const
{
  CGAL_triangulation_assertion( n == NULL);
  return pos == Face_handle();
}
            
template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>::            
operator!=(Nullptr_t n) const
{
  CGAL_triangulation_assertion( n == NULL);
  return !(*this == n);
}
            
template < class Triangulation >
inline bool
Triangulation_line_face_circulator_2<Triangulation>:: 
collinear_outside() const
{
//   return (_tr->is_infinite(*this))
//     && (s == vertex_vertex)
//     && (! _tr->is_infinite((*this)->vertex(i)));
  // return true if  the circulator is non null
  // the line is collinear with a convex hull edge
  // and the convex hull is on the left of the line
  Face_handle fh = pos;
  return ( s == vertex_vertex  &&
	  (! _tr->is_infinite(fh)) &&
           _tr->is_infinite(fh->neighbor(ccw(i))));
}

} //namespace CGAL
#endif //CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H
