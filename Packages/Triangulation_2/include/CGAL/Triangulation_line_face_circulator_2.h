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

  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef std::pair<Face_handle, int>                Edge;
 
  typedef typename Gt::Point Point;
  typedef typename Triangulation::Locate_type Locate_type;

   enum State {undefined = -1,
	       vertex_vertex,
	       vertex_edge,
	       edge_vertex,
	       edge_edge};
            
public:
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
            
  ~Triangulation_line_face_circulator_2()
    {}
            
            
  Triangulation_line_face_circulator_2(const Face_handle& face,
				       int index,
				       State state,
				       const Triangulation_2<Gt,Tds> * t,
				       const Point& pp,
				       const Point& qq)
    : Face_handle(face), _tr(t), s(state), i(index),  
      p(pp), q(qq)            {
    CGAL_triangulation_precondition(! t->geom_traits().compare(p, q));
  }
            
            
  Triangulation_line_face_circulator_2&
  operator=(const Line_face_circulator& lfc)
    {
      ptr() = lfc.ptr();
      i = lfc.i;
      s = lfc.s;
      _tr = lfc._tr;
      p = lfc.p;
      q = lfc.q;
      return *this;
    }

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
  
void increment();
void decrement();
bool locate(const Point& t, 
	    Locate_type &lt, 
	    int &li);



  Line_face_circulator&
  operator++()
    {
      if (ptr()==NULL) {
	return *this; // circulator has singular value
      }
      //xfc while (_tr->is_infinite(ptr()->handle()))
      //    increment();
      //we are looking for a finite face but stop 
      //if back to the origin
      // strange behavoiur
      //why not simply increment and step through infinite faces
      //as other circulators do !!!
      //Face_handle origin = ptr()->handle();
      //do
      increment();
      //while ((_tr->is_infinite(ptr()->handle())) 
      // && (ptr()->handle() != origin));
            
      return *this;
    }
            
            
  Line_face_circulator&
  operator--()
    {
      if (ptr()==NULL) {
	return *this; // circulator has singular value
      }
      //  while (_tr->is_infinite(ptr()->handle()))
      decrement();
      return *this;
    }
            
            
  Line_face_circulator
  operator++(int)
    {
      Line_face_circulator tmp(*this);
      ++(*this);
      return tmp;
    }
            
            
  Line_face_circulator
  operator--(int)
    {
      Line_face_circulator tmp(*this);
      --(*this);
      return tmp;
    }

  bool
  operator==(const Line_face_circulator& lfc) const
    {
      CGAL_triangulation_precondition
	( ptr() != NULL  &&  lfc.ptr() != NULL );
      return ptr() == lfc.ptr();
    }
            
  bool
  operator!=(const Line_face_circulator& lfc) const
    {
      CGAL_triangulation_precondition
	( ptr() != NULL  &&  lfc.ptr() != NULL );
      return ptr() != lfc.ptr();
    }
            
  inline bool
  is_empty()
    {
      return s == undefined;
    }

            
  inline bool
  operator==(CGAL_NULL_TYPE n) const
    {
      CGAL_triangulation_assertion( n == NULL);
      return s == undefined;
    }
            
  inline bool
  operator!=(CGAL_NULL_TYPE n) const
    {
      return !(*this == n);
    }
            
  bool
  collinear_outside() const
    {
            
      return (_tr->is_infinite(ptr()->handle()))
	&& (s == vertex_vertex)
	&& (! _tr->is_infinite(ptr()->vertex(i)));
    }
            
};


template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(Vertex_handle v,
				     const Triangulation_2<Gt,Tds>* tr,
				     const Point& dir)
  : _tr(tr)
{
  CGAL_triangulation_precondition(
				  (! _tr->is_infinite(v)) &&
				  (_tr->dimension() == 2) &&
				  (! _tr->geom_traits().compare(v->point(),dir)));

  p=v->point();
  q=dir;

  //cerr << " p " << p << " q " << q << endl;
		
  Face_circulator fc = v->incident_faces();
  Face_circulator done = fc;

		//cerr  << "(" << fc->vertex(0)->point() << ", "
		//      <<fc->vertex(1)->point() << ", "
		//      << fc->vertex(2)->point() << ")" << endl ;

  int ic = fc->index(v);
  Vertex_handle  vt= fc->vertex(ccw(ic));
  Orientation ptq = RIGHTTURN; 
  //this	initialisation means nothing;
  // just there to avoid a warning
  if (! _tr->is_infinite(vt)) 
    ptq = _tr->geom_traits().orientation(p, vt->point(), q);

		
  while( _tr->is_infinite(vt) || ptq == RIGHTTURN) {
    ++fc;
    if (fc == done) {
      // no edge on the left of pq , pq is a supporting line
      // set ptr() to the right infinite face
      while ( ! _tr->is_infinite(fc)) 
	{  ++fc;}
      ic = fc->index(_tr->infinite_vertex());
      if( _tr->geom_traits().orientation(
					 fc->vertex( cw(ic))->point(),
					 fc->vertex( ccw(ic))->point(),
					 q) != RIGHTTURN) {  ++fc;}
      ptr() = &(*fc);
      i = fc->index(_tr->infinite_vertex());
      s = vertex_vertex;
      return;
    }

    ic = fc->index(v);
    vt= fc->vertex(ccw(ic));
    if (! _tr->is_infinite(vt)) 
      ptq = _tr->geom_traits().orientation(p,  vt->point(), q);
  }
		
	
  // now vt is a finite vertex and ptq is COLLINEAR or LEFTTURN
  Vertex_handle vr = fc-> vertex(cw(ic));
  Orientation prq= RIGHTTURN; 
  //this	initialisation means nothing;
  // just there to avoid a warning
  if (! _tr->is_infinite(vr))
    prq = _tr->geom_traits().orientation(p, vr->point(), q);

  while ( (!_tr->is_infinite(vr)) && (!(prq == RIGHTTURN ))){
    ++fc;
    ic = fc->index(v);
    vr = fc-> vertex(cw(ic));
    if (! _tr->is_infinite(vr))
      prq = _tr->geom_traits().orientation(p, vr->point(), q);
  }

	
  ptr() = &(*fc);
  // reset vt, vt is finite and ptq is still COLLINEAR or LEFTTURN
  ic = fc->index(v);
  vt= fc->vertex(ccw(ic));
  ptq = _tr->geom_traits().orientation(p,  vt->point(), q);

		

  if (_tr->is_infinite(vr)) {		  
    s = vertex_vertex;
    if (ptq == LEFTTURN) {
      i = fc->index(vr); 
    }
    else {// ptq == COLLINEAR
      i= fc->index(vt);
    }
  }
  else{ // vr is a finite vertex}
    if (ptq == LEFTTURN) {
      s = vertex_edge;
      i = ic;
    }
    else { // ptq == COLLINEAR
      s = vertex_vertex;
      i = fc->index(vt);
    }
  }
	
}


template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Triangulation_2<Gt,Tds> * t)
     : _tr(t), s(undefined), p(pp), q(qq)
{
  Vertex_handle inf = _tr->infinite_vertex();
  Face_circulator fc = inf->incident_faces(),
    done(fc);
  i = fc->index(inf);
            
  Point l = fc->vertex(cw(i))->point(),
    r = fc->vertex(ccw(i))->point();
            
  Orientation pql = _tr->geom_traits().orientation(p, q, l),
    pqr = _tr->geom_traits().orientation(p, q, r);
            
            
  do{
    if( (pql == LEFTTURN) && (pqr == RIGHTTURN) ){
      *this = ++Line_face_circulator( fc->handle() ,
				      i,
				      vertex_edge,
				      t,
				      p,
				      q);
      return;
    } else if ( (pql == LEFTTURN) && 
		(pqr == COLLINEAR) ){
      *this = ++Line_face_circulator( fc->handle() ,
				      ccw(i),
				      vertex_vertex,
				      t,
				      p,
				      q);
      return;
    } else if( (pql == COLLINEAR) && 
	       (pqr == COLLINEAR) ){
      Face_handle n = fc->neighbor(i);
      int ni  = n->index( fc->handle() );
      Vertex_handle vn = n->vertex(ni);
      if(_tr->geom_traits().orientation(p, q, vn->point()) ==
	 LEFTTURN){
	// the entire triangulation is to the left of line (p,q).
	// There might be further collinear edges, so we might have
	// to walk back on the hull.
	while(1){
	  ++fc;
	  i = fc->index(inf);
	  l = fc->vertex(cw(i))->point();
	  if(_tr->geom_traits().orientation(p, q, l) == 
	     COLLINEAR){
	    continue;
	  } else {
	    // we went one step to far back
	    --fc;
	    i = fc->index(inf);
	    ptr() = &(*fc->neighbor(i));
	    i = cw(ptr()->index( fc->handle() ));
	    s = vertex_vertex;
	    return;
	  }
	}
      } else {
	// the entire triangulation is to the right of line (p,q).
	// here are no faces to traverse, so we give the circulator
	// a singular value
	return;
      }
    } else {
      --fc;
      l = r;
      pql = pqr;
      i = fc->index(inf);
      r = fc->vertex(ccw(i))->point();
      pqr = _tr->geom_traits().orientation(p, q, r);
    }
  }while(fc != done);
  // if line (p,q) does not intersect the convex hull in an edge
  // the circulator has a singular value
}


template < class Gt, class Tds >
Triangulation_line_face_circulator_2<Gt,Tds>::
Triangulation_line_face_circulator_2(const Point& pp,
				     const Point& qq,
				     const Face_handle& ff,
				     const Triangulation_2<Gt,Tds>* t)
  : Face::Face_handle(ff), _tr(t), s(undefined), p(pp), q(qq)
{
  CGAL_triangulation_precondition(_tr->is_infinite(ff) ||
				  _tr->oriented_side(ff,p) != ON_NEGATIVE_SIDE);
            
  int j;
  if(_tr->is_infinite(ptr()->handle())){
    *this  = Line_face_circulator(p, q, t);
    return;
  }
            
		
  // Test whether p lies on a vertex
  for(j = 0; j < 3; j++){
    if(ptr()->vertex(j)->point() == p){
      *this = Line_face_circulator( ptr()->vertex(j), t, q);
      return;
    }
  }
            
  // Test whether p lies on an edge
  for(j = 0; j < 3; j++){
    if(_tr->geom_traits().orientation(ptr()->vertex(j)->point(),
				      ptr()->vertex(ccw(j))->point(),
				      p) == COLLINEAR){
      Orientation jpq =
	_tr->geom_traits().orientation(ptr()->vertex(j)->point(),
				       p,
				       q);
      Orientation p_cwj_q =
	_tr->geom_traits().orientation(p,
				       ptr()->vertex(cw(j))->point(),
				       q);
      switch(jpq){
      case COLLINEAR:
	if(p_cwj_q == RIGHTTURN){
	  s = vertex_vertex;
	  i = ccw(j);
	  return;
	} 
	else if(! _tr->is_infinite(ptr()->neighbor(cw(j)))){
	  Face_handle n = ptr()->neighbor(cw(j));
	  i = cw(n->index(ptr()->handle()));
	  ptr() = &(*n);
	  s = vertex_vertex;
	  return;
	} else {
                                // singular value
	  return;
	}
      case RIGHTTURN:
	i = cw(j);
	s = (p_cwj_q == COLLINEAR) ? vertex_edge :  
	  edge_edge;
	break;
      default: //  LEFTTURN
	switch(p_cwj_q){
	case COLLINEAR:
	  s = edge_vertex;
	  i = cw(j);
	  return;
	case RIGHTTURN:
	  s = edge_edge;
	  i = j;
	  return;
	default:
	  s = edge_edge;
	  i = ccw(j);
	  return;
	}
      }
    }
  }
            
  // p lies in the interior of the face
  Orientation or[3];
  for(j=0; j<3; j++){
    or[j] =
      _tr->geom_traits().orientation(p,q,ptr()->vertex(j)->point());
  }
  for(j=0; j<3; j++){
    if(or[j] == COLLINEAR){
      i = j;
      s = (or[ccw(j)] == LEFTTURN) ? edge_vertex : 
	vertex_edge;
      return;
    }
  }
  s = edge_edge;
  for(j=0; j<3; j++){
    if(or[j] == RIGHTTURN){
      i = (or[ccw(j)] == RIGHTTURN) ? j : cw(j);
      return;
    }
  }
}


template < class Gt, class Tds >
void
Triangulation_line_face_circulator_2<Gt,Tds>::
increment()
{
  CGAL_triangulation_precondition(s != undefined);
  if(s == vertex_vertex || s == edge_vertex){
    Orientation o;
    Point r;
    do{
      Face_handle n = ptr()->neighbor(cw(i));
      i = n->index(ptr()->handle());
      ptr() = &(*n);
      if (n->vertex(i) == _tr->infinite_vertex()){
	o = COLLINEAR;
	i = cw(i);
	break;
      }
      r = n->vertex(i)->point();
      i = cw(i);
    }while((o = _tr->geom_traits().orientation(p, q, r)) == 
	   LEFTTURN);
            
    if(o == COLLINEAR){
      s = vertex_vertex;
      i = ccw(i);
    } else {
      s = vertex_edge;
    }
  } else {
    Face_handle n = ptr()->neighbor(i);
    int ni = n->index(ptr()->handle());
    ptr() = &(*n);
    Orientation o = _tr->is_infinite(ptr()->vertex(ni)) ?
      COLLINEAR :
      _tr->geom_traits().orientation(p,q,ptr()->vertex(ni)->point());
            
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
                CGAL_triangulation_precondition(s != undefined);
                if(s == vertex_vertex || s == vertex_edge){
                    if(s == vertex_vertex){
                        i = cw(i);
                    }
                    Orientation o;
                    Point r;
                    do{
                        Face_handle n = ptr()->neighbor(ccw(i));
                        i = n->index(ptr()->handle());
                        ptr() = &(*n);
                        if (n->vertex(i) == _tr->infinite_vertex()){
                            o = COLLINEAR;
                            i = ccw(i);
                            break;
                        }
                        r = n->vertex(i)->point();
                        i = ccw(i);
                    }while((o = _tr->geom_traits().orientation(p, q, r)) == 
                       LEFTTURN);
            
                    s = (o == COLLINEAR) ? vertex_vertex : edge_vertex;
            
                } else { // s == edge_edge  ||  s == edge_vertex
             // the following is not nice. A better solution is to say
          // that index i is at the vertex that is alone on one side of l(p,q)
                    if(s == edge_edge){
                        i = (_tr->geom_traits().orientation
                                            (p, q,
                                            ptr()->vertex(i)->point()) == 
						LEFTTURN)
                            ? cw(i) : ccw(i);
                    }
                    Face_handle n = ptr()->neighbor(i);
                    i = n->index(ptr()->handle());
                    ptr() = &(*n);
                    Orientation o = _tr->is_infinite(ptr()->vertex(i)) ?
                        COLLINEAR :
               _tr->geom_traits().orientation(p, q, ptr()->vertex(i)->point());
            
                    s = (o == COLLINEAR) ? vertex_edge : edge_edge;
                }
}

template < class Gt, class Tds >
bool
Triangulation_line_face_circulator_2<Gt,Tds>::
locate(const Point& t,
       Locate_type &lt,
       int &li)
{
  switch(s){
            
  case edge_edge:
  case vertex_edge:
    {
      Orientation o =
	_tr->geom_traits().orientation(ptr()->vertex(ccw(i))->point(),
				       ptr()->vertex(cw(i))->point(),
				       t);
      if(o == RIGHTTURN){
	return false;
      }
      if(o == COLLINEAR){
	lt = Triangulation::EDGE;
	li = i;
	return true;
      }
      lt = Triangulation::FACE;
      return true;
    }
  case vertex_vertex:
    {
      if(_tr->is_infinite(ptr()->vertex(i))){
	CGAL_triangulation_assertion(
	   _tr->geom_traits().orientation( ptr()->vertex(cw(i))->point(),
					   ptr()->vertex(ccw(i))->point(),
					   t) != LEFTTURN);
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
            
      Point u = ptr()->vertex(cw(i))->point();
      Point v = ptr()->vertex(i)->point();
      // u == t  was detected earlier
      if(_tr->geom_traits().compare_x(v,t)==EQUAL && 
	 _tr->geom_traits().compare_y(v,t)==EQUAL){
	lt = Triangulation::VERTEX;
	li = i;
	return true;
      }
      if(_tr->collinear_between(u, t, v)){
	lt = Triangulation::EDGE;
	li = ccw(i);
	return true;
      }
      return false;
    }
  default: // edge_vertex
    {
      if(_tr->is_infinite(ptr()->vertex(i))){
	lt = Triangulation::OUTSIDE_CONVEX_HULL;
	li = i;
	return true;
      }
      if(_tr->geom_traits().compare_x(t,ptr()->vertex(i)
				      ->point())==EQUAL &&
	 _tr->geom_traits().compare_y(t,ptr()->vertex(i)
				      ->point())==EQUAL ){
	li = i;
	lt = Triangulation::VERTEX;
	return true;
      }
      if(_tr->collinear_between(p, t, ptr()->vertex(i)
				->point())){
	lt = Triangulation::FACE;
	return true;
      }
      return false;
    }
  }
}
           


CGAL_END_NAMESPACE
#endif //CGAL_TRIANGULATION_LINE_FACE_CIRCULATOR_2_H

