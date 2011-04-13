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
// file          : include/CGAL/Triangulation_line_face_circulator_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H
#define CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H

// #include <CGAL/circulator.h>
// #include <CGAL/Triangulation_utils_2.h>
// #include <CGAL/triangulation_assertions.h>
// #include <CGAL/Triangulation_face_2.h>
// #include <CGAL/Triangulation_vertex_2.h>
// #include <CGAL/Triangulation_handles_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_2;

template < class Gt, class Tds >
class Triangulation_line_face_circulator_2
  :   public Bidirectional_circulator_base<Triangulation_face_2<Gt,Tds>,
	                                      std::ptrdiff_t,std::size_t>,
      public Triangulation_cw_ccw_2,
      public  Triangulation_face_handle_2<Gt,Tds>
{
public:
  typedef Triangulation_line_face_circulator_2<Gt,Tds> Line_face_circulator;
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  
  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef std::pair<Face_handle, int>                Edge;
 
  typedef typename Gt::Point_2 Point;
  typedef typename Triangulation::Locate_type Locate_type;

   enum State {undefined = -1,
	       vertex_vertex,
	       vertex_edge,
	       edge_vertex,
	       edge_edge};
            
private:
  const Triangulation_2<Gt, Tds>* _tr;
  State s;
  int i;
  Point p, q;
            
public:
  Triangulation_line_face_circulator_2()
    : Face_handle(), _tr(NULL), s(undefined), i(-1)
    {}
            
  Triangulation_line_face_circulator_2(const Line_face_circulator& lfc)
    : Face_handle(& (*lfc)), _tr(lfc._tr), s(lfc.s), 
					i(lfc.i),  p(lfc.p), q(lfc.q)
    {}
            
  Triangulation_line_face_circulator_2(Vertex_handle v,
				       const Triangulation_2<Gt,Tds>* tr,
				       const Point& dir);

  Triangulation_line_face_circulator_2(const Point& pp,
				       const Point& qq,
				       const Triangulation_2<Gt,Tds> * t);
           
  Triangulation_line_face_circulator_2(const Point& pp,
				       const Point& qq,
				       const Face_handle& ff,
				       const Triangulation_2<Gt,Tds>* t); 

  Triangulation_line_face_circulator_2&
  operator=(const Line_face_circulator& lfc) ;
  

  Line_face_circulator&   operator++() ;
  Line_face_circulator&   operator--() ;
  Line_face_circulator    operator++(int);
  Line_face_circulator    operator--(int);
  bool  operator==(const Line_face_circulator& lfc) const;
  bool  operator!=(const Line_face_circulator& lfc) const;
  bool  is_empty() const;
  bool  operator==(CGAL_NULL_TYPE n) const;
  bool  operator!=(CGAL_NULL_TYPE n) const;
  bool  collinear_outside() const;
  bool locate(const Point& t, Locate_type &lt,  int &li);

  //private:
  Triangulation_line_face_circulator_2(const Face_handle& face,
				       int index,
				       State state,
				       const Triangulation_2<Gt,Tds> * t,
				       const Point& pp,
				       const Point& qq);
private:
  void increment();
  void decrement();
};


template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(const Face_handle& face,
				     int index,
				     State state,
				     const Triangulation_2<Gt,Tds> * t,
				     const Point& pp,
				     const Point& qq)
    : Face_handle(face), _tr(t), s(state), i(index),  
      p(pp), q(qq)            
{
  CGAL_triangulation_precondition(! t->xy_equal(p, q));
}

template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>&
Triangulation_line_face_circulator_2<Gt,Tds>::
operator=(const Line_face_circulator& lfc)
{
  Face_handle::operator=(lfc);
  i = lfc.i;
  s = lfc.s;
   _tr = lfc._tr;
  p = lfc.p;
  q = lfc.q;
  return *this;
}

template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(Vertex_handle v,
				     const Triangulation_2<Gt,Tds>* tr,
				     const Point& dir)
  :Face_handle(NULL), _tr(tr), s(undefined)
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
  Face_circulator fc = v->incident_faces();
  Face_circulator done(fc);
  int ic = fc->index(v);
  Vertex_handle  vt= fc->vertex(cw(ic));
   while( _tr->is_infinite(vt) || 
	  _tr->orientation(p, q, vt->point()) != LEFTTURN) {
    ++fc;
    ic = fc->index(v);
    vt= fc->vertex(cw(ic));
    if (fc == done) { *this = Line_face_circulator(); return;}
  }
  
  // now vt is finite and to the left of pq
  Vertex_handle vr = fc-> vertex(ccw(ic));
  Orientation pqr = RIGHTTURN;
  // this initialisation means nothing; just there to avoid the
  // warning "pqr might be used uninitialized"
  while ( (!_tr->is_infinite(vr)) &&
	  (pqr = _tr->orientation(p, q, vr->point()))== LEFTTURN ) {
    --fc;
    ic = fc->index(v);
    vr = fc-> vertex(ccw(ic));
  }

  // vr can be infinite or finite. 
  // If finite [pqr] is COLLINEAR or RIGHTTURN
  // reset vt and conclude.  vt is still finite and [pqt] still LEFTTURN
  ic = fc->index(v);
  vt= fc->vertex(cw(ic));
  CGAL_triangulation_assertion (_tr->orientation(p,q, vt->point())==
				LEFTTURN );
  if (_tr->is_infinite(vr)) {
    --fc;
    ic = fc->index(v);
    vr = fc->vertex(ccw(ic));
    pqr = _tr->orientation(p, q, vr->point());
    switch(pqr){
    case RIGHTTURN:
    case COLLINEAR:
      ++fc;
      ic = fc->index(_tr->infinite_vertex());
      Face_handle::operator=(&(*fc));
      s = vertex_vertex;
      i = ic;
      break;
    case LEFTTURN:
     *this = Line_face_circulator(); 
     break;
    }
  }
  else if (pqr == COLLINEAR) {
    Face_handle::operator=(&(*fc));
    s = vertex_vertex;
    i = ccw(ic);
   }
  else { // pqr==RIGHTTURN 
    Face_handle::operator=(&(*fc));
    s = vertex_edge;
    i = ic ;
  }
  return;
}
  



template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Triangulation_2<Gt,Tds> * t)
     : Face_handle(NULL), _tr(t), s(undefined), p(pp), q(qq)
  //begins at the  first finite face traversed be the oriented line pq
{
  Vertex_handle inf = _tr->infinite_vertex();
  Face_circulator fc = inf->incident_faces(),
    done(fc);

  i = fc->index(inf);
  Point l = fc->vertex(cw(i))->point();
  Point      r = fc->vertex(ccw(i))->point();
  Orientation pql = _tr->orientation(p, q, l);
  Orientation pqr = _tr->orientation(p, q, r);
            
   do{
    if( (pql == LEFTTURN) && (pqr == RIGHTTURN) ){
      *this = ++Line_face_circulator( fc, i, vertex_edge, t, p, q);
           return;
    } 
    else if ( (pql == LEFTTURN) && (pqr == COLLINEAR) ){
      --fc;
      i = fc->index(inf);
      Point s = fc->vertex(ccw(i))->point();
      Orientation pqs  = _tr->orientation(p, q, s);
      Face_handle fn;
      int in;
      switch(pqs) {
      case LEFTTURN:
	*this = Line_face_circulator();
	return;
      case COLLINEAR:
	fn = fc->neighbor(i);
	in = fn->index(fc);
	*this = Line_face_circulator( fn, cw(in),vertex_vertex,t,p,q);
	return;
      case RIGHTTURN:
	fn = fc->neighbor(i);
	Vertex_handle vr = fc->vertex(cw(i)); // vertex corresponding to r
	in = fn->index(vr);
	s = fn->vertex(cw(in))->point();
	pqs = _tr->orientation(p, q, s);
	Orientation pqss = RIGHTTURN;
	while ( pqs != LEFTTURN) {
	  pqss = pqs;
	  fn = fn->neighbor(ccw(in));
	  in = fn->index(vr);
	  s = fn->vertex(cw(in))->point();
	  pqs = _tr->orientation(p, q, s);
	}
	if (pqss == RIGHTTURN)
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


template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Face_handle& ff,
				     const Triangulation_2<Gt,Tds>* t)
  : Face_handle(ff), _tr(t), s(undefined), p(pp), q(qq)
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
  if(_tr->is_infinite(*this)){
    *this  = Line_face_circulator(p, q, t);
    return;
  }

  // Test whether p lies on a vertex
  for(j = 0; j < 3; j++){
    if((*this)->vertex(j)->point() == p){
      *this = Line_face_circulator( (*this)->vertex(j), t, q);
      if( (!is_empty()) && _tr->is_infinite(*this )) --(*this);
      return;
    }
  }
            
  // Test whether p lies on an edge
  for(j = 0; j < 3; j++) {
    if(_tr->orientation((*this)->vertex(ccw(j))->point(),
			(*this)->vertex(cw(j))->point(),
			p) == COLLINEAR){
      Orientation pqj =
	_tr->orientation(p, q, (*this)->vertex(j)->point());
      Orientation pqcwj =
	_tr->orientation(p, q, (*this)->vertex(cw(j))->point());
      switch(pqcwj) {
      case COLLINEAR :
	if(pqj == LEFTTURN){
	  s = vertex_vertex;
	  i = cw(j);
	  return;
	} 
	else if(! _tr->is_infinite((*this)->neighbor(j))){
	  Face_handle n = (*this)->neighbor(j);
	  i = cw(n->index(*this));
	  Face_handle::operator=(n);
	  s = vertex_vertex;
	  return;
	} else {
           // singular value
	  *this = Line_face_circulator();
	  return;
	}
      case LEFTTURN :
	i = j;
	s = (pqj == COLLINEAR) ? vertex_edge :  
	  edge_edge;
	break;
      case RIGHTTURN :
	switch(pqj){
	case COLLINEAR:
	  s = edge_vertex;
	  i = j;
	  return;
	case LEFTTURN:
	  s = edge_edge;
	  i = cw(j);
	  return;
	case RIGHTTURN:
	  s = edge_edge;
	  i = ccw(j);
	  return;
	}
      }
    }
  }

//  for(j = 0; j < 3; j++) {
//     if(_tr->orientation((*this)->vertex(j)->point(),
// 			(*this)->vertex(ccw(j))->point(),
// 			p) == COLLINEAR){
//       Orientation jpq =
// 	_tr->orientation((*this)->vertex(j)->point(), p, q);
//       Orientation p_cwj_q =
// 	_tr->orientation(p, (*this)->vertex(cw(j))->point(), q);
//       switch(jpq) {
//       case COLLINEAR:
// 	if(p_cwj_q == RIGHTTURN){
// 	  s = vertex_vertex;
// 	  i = ccw(j);
// 	  return;
// 	} 
// 	else if(! _tr->is_infinite((*this)->neighbor(cw(j)))){
// 	  Face_handle n = (*this)->neighbor(cw(j));
// 	  i = cw(n->index(*this));
// 	  Face_handle::operator=(n);
// 	  s = vertex_vertex;
// 	  return;
// 	} else {
//            // singular value
// 	  *this = Line_face_circulator();
// 	  return;
// 	}
//       case RIGHTTURN:
// 	i = cw(j);
// 	s = (p_cwj_q == COLLINEAR) ? vertex_edge :  
// 	  edge_edge;
// 	break;
//       default: //  LEFTTURN
// 	switch(p_cwj_q){
// 	case COLLINEAR:
// 	  s = edge_vertex;
// 	  i = cw(j);
// 	  return;
// 	case RIGHTTURN:
// 	  s = edge_edge;
// 	  i = j;
// 	  return;
// 	default:
// 	  s = edge_edge;
// 	  i = ccw(j);
// 	  return;
// 	}
//       }
//     }
//   }
            
  // p lies in the interior of the face
  Orientation orient[3];
  for(j=0; j<3; j++) {
    orient[j] =
      _tr->orientation(p,q,(*this)->vertex(j)->point());
  }
  for(j=0; j<3; j++) {
    if(orient[j] == COLLINEAR) {
      i = j;
      s = (orient[ccw(j)] == LEFTTURN) ? edge_vertex : 
	vertex_edge;
      return;
    }
  }
  s = edge_edge;
  for(j=0; j<3; j++){
    if(orient[j] == RIGHTTURN){
      i = (orient[ccw(j)] == RIGHTTURN) ? j : cw(j);
      return;
    }
  }
}


template < class Gt, class Tds >
void
Triangulation_line_face_circulator_2<Gt,Tds>::
increment()
{
  CGAL_triangulation_precondition(*this != CGAL_CIRC_NULL);
  if(s == vertex_vertex || s == edge_vertex) {
    Orientation o;
    Point r;
    do{
      Face_handle n = (*this)->neighbor(cw(i));
      i = n->index(*this);
      Face_handle::operator=(n);
      if (n->vertex(i) == _tr->infinite_vertex()){
	o = COLLINEAR;
	i = cw(i);
	break;
      }
      r = n->vertex(i)->point();
      i = cw(i);
    }while((o = _tr->orientation(p, q, r)) ==  LEFTTURN);
            
    if(o == COLLINEAR) {
      s = vertex_vertex;
      i = ccw(i);
    } 
    else {
      s = vertex_edge;
    }
  } 
  else {
    Face_handle n = (*this)->neighbor(i);
    int ni = n->index(*this);
    Face_handle::operator=(n);
    Orientation o = _tr->is_infinite((*this)->vertex(ni)) ?
      COLLINEAR :
      _tr->orientation(p,q,(*this)->vertex(ni)->point());
            
    switch(o){
    case LEFTTURN:
      s = edge_edge;
      i = ccw(ni);
      break;
    case RIGHTTURN:
      s = edge_edge;
      i = cw(ni);
      break;
    default:
      s = edge_vertex;
      i = ni;
    }
  }
} 
            
template < class Gt, class Tds >
void
Triangulation_line_face_circulator_2<Gt,Tds>::             
decrement()
{
  CGAL_triangulation_precondition(*this != CGAL_CIRC_NULL);
  if(s == vertex_vertex || s == vertex_edge) {
    if(s == vertex_vertex){
      i = cw(i);
    }
    Orientation o;
    Point r;
    do{
      Face_handle n = (*this)->neighbor(ccw(i));
      i = n->index(*this);
      Face_handle::operator=(n);
      if (n->vertex(i) == _tr->infinite_vertex()){
	o = COLLINEAR;
	i = ccw(i);
	break;
      }
      r = n->vertex(i)->point();
      i = ccw(i);
    }while((o = _tr->orientation(p, q, r)) == 
	   LEFTTURN);
            
    s = (o == COLLINEAR) ? vertex_vertex : edge_vertex;
            
  } 
  else { // s == edge_edge  ||  s == edge_vertex
    // the following is not nice. A better solution is to say
    // that index i is at the vertex that is alone on one side of l(p,q)
    if(s == edge_edge){
      i = (_tr->orientation
	   (p, q,
	    (*this)->vertex(i)->point()) == 
	   LEFTTURN)
	? cw(i) : ccw(i);
    }
    Face_handle n = (*this)->neighbor(i);
    i = n->index(*this);
    Face_handle::operator=(n);
    Orientation o = _tr->is_infinite((*this)->vertex(i)) ?
      COLLINEAR :
      _tr->orientation(p, q, (*this)->vertex(i)->point());
            
    s = (o == COLLINEAR) ? vertex_edge : edge_edge;
  }
}

template < class Gt, class Tds >
bool
Triangulation_line_face_circulator_2<Gt,Tds>::
locate(const Point& t, Locate_type &lt,  int &li)
{
  switch(s){            
  case edge_edge:
  case vertex_edge:
    {
      Orientation o =
	_tr->orientation((*this)->vertex(ccw(i))->point(),
			 (*this)->vertex(cw(i))->point(),
			 t);
      if(o == RIGHTTURN)      return false;
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
      if(_tr->is_infinite((*this)->vertex(i))){
	CGAL_triangulation_assertion(
	       _tr->orientation( (*this)->vertex(cw(i))->point(),
				 (*this)->vertex(ccw(i))->point(),
				 t) != LEFTTURN);
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
      Point u = (*this)->vertex(cw(i))->point();
      Point v = (*this)->vertex(i)->point();
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
      if(_tr->is_infinite((*this)->vertex(i))){
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
      if(_tr->xy_equal(t,(*this)->vertex(i) ->point()) ){
	li = i;
	lt = Triangulation::VERTEX;
	return true;
      }
      if(_tr->collinear_between(p, t, (*this)->vertex(i)->point())) {
	lt = Triangulation::FACE;
	return true;
      }
      return false;
    }
  }
}
           
template < class Gt, class Tds >
inline
Triangulation_line_face_circulator_2<Gt,Tds>&
Triangulation_line_face_circulator_2<Gt,Tds>::
operator++()
{
  CGAL_triangulation_precondition( *this != CGAL_CIRC_NULL) ;
  increment();
  return *this;
}
            
template < class Gt, class Tds >
inline
Triangulation_line_face_circulator_2<Gt,Tds>&
Triangulation_line_face_circulator_2<Gt,Tds>::            
operator--()
{
  CGAL_triangulation_precondition(*this != CGAL_CIRC_NULL) ;
  decrement();
  return *this;
}
            
template < class Gt, class Tds >
inline
Triangulation_line_face_circulator_2<Gt,Tds>
Triangulation_line_face_circulator_2<Gt,Tds>::             
operator++(int)
{
  Line_face_circulator tmp(*this);
  ++(*this);
  return tmp;
}
            
template < class Gt, class Tds >
inline
Triangulation_line_face_circulator_2<Gt,Tds>
Triangulation_line_face_circulator_2<Gt,Tds>::             
operator--(int)
{
  Line_face_circulator tmp(*this);
  --(*this);
  return tmp;
}

template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>::    
operator==(const Line_face_circulator& lfc) const
{
  CGAL_triangulation_precondition( *this != CGAL_CIRC_NULL &&
				   lfc != CGAL_CIRC_NULL);
  return ( Face_handle::operator==(lfc) &&  _tr == lfc._tr &&
            s== lfc.s && p==lfc.p && q==lfc.q);
}

template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>:: 
operator!=(const Line_face_circulator& lfc) const
{
  return !(*this == lfc);
}
            
template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>::   
is_empty() const
{
  return Face_handle::operator==(CGAL_CIRC_NULL);
}

template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>::            
operator==(CGAL_NULL_TYPE n) const
{
  CGAL_triangulation_assertion( n == NULL);
  return Face_handle::operator==(CGAL_CIRC_NULL);
}
            
template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>:: 
operator!=(CGAL_NULL_TYPE n) const
{
  return !(*this == n);
}
            
template < class Gt, class Tds >
inline bool
Triangulation_line_face_circulator_2<Gt,Tds>:: 
collinear_outside() const
{
//   return (_tr->is_infinite(*this))
//     && (s == vertex_vertex)
//     && (! _tr->is_infinite((*this)->vertex(i)));
  // return true if  the circulator is non null
  // the line is collinear with a convex hull edge
  // and the convex hull is on the left of the line
  Face_handle fh = *this;
  return ( s == vertex_vertex  &&
	  (! _tr->is_infinite(fh)) &&
           _tr->is_infinite(fh->neighbor(ccw(i))));
}

CGAL_END_NAMESPACE
#endif //CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H

