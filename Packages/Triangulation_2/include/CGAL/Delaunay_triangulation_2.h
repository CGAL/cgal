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
// release       : $CGAL_Revision: CGAL-2.0-I-12 $
// release_date  : $CGAL_Date: 1999/04/28 $
//
// file          : include/CGAL/Delaunay_triangulation_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_DELAUNAY_TRIANGULATION_2_H
#define CGAL_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Triangulation_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Delaunay_triangulation_2 : public Triangulation_2<Gt,Tds>
{
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Distance Distance;
  typedef typename Geom_traits::Ray Ray;
  typedef typename Geom_traits::Line Line;
  typedef typename Geom_traits::Direction Direction;

  typedef Triangulation_2<Gt,Tds>              Triangulation;
  typedef typename Triangulation::Locate_type           Locate_type;
  typedef typename Triangulation::Face_handle           Face_handle;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;
  typedef typename Triangulation::Edge                  Edge;
  typedef typename Triangulation::Edge_circulator       Edge_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  
  Delaunay_triangulation_2(const Gt& gt = Gt())
  : Triangulation_2<Gt,Tds>(gt) {}
  
  Delaunay_triangulation_2(
	       const Delaunay_triangulation_2<Gt,Tds> &tr)
      : Triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition( is_valid() );  }
  
// CHECK -QUERY
  bool is_valid(bool verbose = false, int level = 0) const;

  Oriented_side
  side_of_oriented_circle(Face_handle f, const Point & p) const;
  
  Vertex_handle
  nearest_vertex(const Point& p, Face_handle f= Face_handle()) const;
   
  // DUAL
  Point dual (Face_handle f) const;
  Object dual(const Edge &e) const ;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;
 
  //INSERTION-REMOVAL
  Vertex_handle  insert(const Point  &p, Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  Vertex_handle push_back(const Point &p);

  void  remove(Vertex_handle v );
  

  
private:
  void restore_Delaunay(Vertex_handle v);
  void propagating_flip(Face_handle& f,int i);
  void remove_2D(Vertex_handle v );
  void fill_hole( Vertex_handle v, std::list<Edge> & hole);

  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const;
  Vertex_handle nearest_vertex_1D(const Point& p) const;

  void  look_nearest_neighbor(const Point& p,
			      Face_handle f,
			      int i,
			      int& min,
			      Vertex_handle& nn,
			      Distance& closer) const;

public:
  template < class Stream>
  Stream& draw_dual(Stream & ps)
    {
      Finite_edges_iterator eit= finite_edges_begin();
      for (; eit != finite_edges_end(); ++eit) {
	Object o = dual(eit);
	Ray r;
	Segment s;
	Line l;
	if (CGAL::assign(s,o)) ps << s;
	if (CGAL::assign(r,o)) ps << r;
	if (CGAL::assign(l,o)) ps << l;
      }
      return ps;
    }

  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
    {
      int n = number_of_vertices();
      while(first != last){
	insert(*first);
	++first;
      }
      return number_of_vertices() - n;
    }

};


template < class Gt, class Tds >
bool
Delaunay_triangulation_2<Gt,Tds>::
is_valid(bool verbose = false, int level = 0) const
{
  bool result = Triangulation_2<Gt,Tds>::is_valid(verbose, level);

  for( Finite_faces_iterator it = finite_faces_begin(); 
       it != finite_faces_end() ; it++) {
    for(int i=0; i<3; i++) {
      if ( ! is_infinite( it->mirror_vertex(i))) {
	result = result &&  ON_POSITIVE_SIDE != 
	  side_of_oriented_circle( it, it->mirror_vertex(i)->point());
      }
      CGAL_triangulation_assertion( result );
    }
  }
  return result;
}

template < class Gt, class Tds >
Oriented_side
Delaunay_triangulation_2<Gt,Tds>::
side_of_oriented_circle(Face_handle f, const Point & p) const
{
  if ( ! is_infinite(f) ) {
    return geom_traits().side_of_oriented_circle(f->vertex(0)->point(),
						 f->vertex(1)->point(),
						 f->vertex(2)->point(),p);
  }

  int i = f->index(infinite_vertex());
  Orientation o =
    geom_traits().orientation(f->vertex(ccw(i))->point(),
			      f->vertex(cw(i))->point(),
			      p);
						     
  return (o == NEGATIVE) ? ON_NEGATIVE_SIDE :
                (o == POSITIVE) ? ON_POSITIVE_SIDE :
                      ON_ORIENTED_BOUNDARY;
}

template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex(const Point  &p, Face_handle f) const
{
  switch (dimension()) {
  case 0:
    if (number_of_vertices() == 0) return NULL;
    if (number_of_vertices() == 1) return finite_vertex();
    break;
  case 1:
    return nearest_vertex_1D(p);
    break;      
  case 2:
    return nearest_vertex_2D(p,f);
    break;
  }
  return NULL;
}
  
template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_2D(const Point& p, Face_handle f) const
{
  CGAL_triangulation_precondition(dimension() == 2);
  Vertex_handle nn;
  if (f== Face_handle()) f = locate(p);
  else CGAL_triangulation_precondition(oriented_side(f,p) != ON_NEGATIVE_SIDE);
  Distance closer(p,&geom_traits());
  int min;
  int i;
  
  i = ( ! is_infinite(f->vertex(0)) ) ? 0 : 1;
  closer.set_point(1,f->vertex(i)->point());
  min = 1;
  nn = f->vertex(i);
  if ( ! is_infinite(f->vertex(ccw(i)))){
    closer.set_point( 3-min, f->vertex(ccw(i))->point() );
    if (  ( (min==1)? LARGER : SMALLER )
	  == closer.compare() ) {
      min = 3-min;
      nn=f->vertex(ccw(i));
    }
  }
  if ( ! is_infinite(f->vertex(cw(i)))){
    closer.set_point( 3-min, f->vertex(cw(i))->point() );
    if (  ( (min==1)? LARGER : SMALLER )
	  == closer.compare() ) {
      min = 3-min;
      nn=f->vertex(cw(i));
    }
  }
  look_nearest_neighbor(p,f,0,min,nn,closer);
  look_nearest_neighbor(p,f,1,min,nn,closer);
  look_nearest_neighbor(p,f,2,min,nn,closer);
  return nn;
}

template < class Gt, class Tds >
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_1D(const Point& p) const
{
  Vertex_handle nn;
  Distance closer(p,&geom_traits());
  int min;

  Finite_vertices_iterator vit=finite_vertices_begin();
  closer.set_point(1,vit->point());
  min = 1;
  nn = vit->handle();
  do {
    closer.set_point( 3-min, (++vit)->point());
     if (  ( (min==1)? LARGER : SMALLER )
	  == closer.compare() ) {
      min = 3-min;
      nn=vit->handle();
    }
  }while( vit != finite_vertices_end());
  return nn;
}
  
template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
look_nearest_neighbor(const Point& p,
                      Face_handle f,
		      int i,
		      int& min,
		      Vertex_handle& nn,
		      Distance& closer) const
{
  Face_handle  ni=f->neighbor(i);
  if ( ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p) ) {
    return;
  }
  i = ni->index(f);
  if ( ! is_infinite(ni->vertex(i))){
    closer.set_point( 3-min, ni->vertex(i)->point() );
    if (  ( (min==1)? LARGER : SMALLER )
	  == closer.compare() ) {
      min = 3-min;
      nn=ni->vertex(i);
    }
  }
  // recursive exploration of triangles whose circumcircle contains p
  look_nearest_neighbor(p, ni, ccw(i), min, nn, closer);
  look_nearest_neighbor(p, ni, cw(i),  min, nn, closer);
} 

//DUALITY
template<class Gt, class Tds>
typename Gt::Point
Delaunay_triangulation_2<Gt,Tds>::
dual (Face_handle f) const
{
  CGAL_triangulation_precondition (dimension()==2);
  return geom_traits().circumcenter(f->vertex(0)->point(),
				    f->vertex(1)->point(),
				    f->vertex(2)->point());
}
  
template < class Gt, class Tds >
Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Edge &e) const
{
  CGAL_triangulation_precondition (!is_infinite(e));
  if( dimension()== 1 ){
    Line l = geom_traits().bisector(segment(e)).opposite();
    return Object(new Wrapper< Line >(l));
  }
  
  // dimension==2
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Segment s(dual(e.first),dual(e.first->neighbor(e.second)));
    return Object(new Wrapper< Segment >(s));
  } 
  // one of the adjacent face is infinite
  Face_handle f; int i;
  if (is_infinite(e.first)) {
    f=e.first->neighbor(e.second); f->has_neighbor(e.first,i);
  } 
  else {
    f=e.first; i=e.second;
  }
  Line l = geom_traits().bisector(segment(f,i)).opposite();
  Ray r(dual(f),l.direction());
  return Object(new Wrapper< Ray >(r));
}
  
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::  
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}
  
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Finite_edges_iterator& ei) const
{
  return dual(*ei);
}
  
template < class Gt, class Tds >
inline
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p,  Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate (p, lt, li, start);
  return insert(p, lt, loc, li);
}
  
template < class Gt, class Tds >
inline
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  return insert(p);
}
  
template < class Gt, class Tds >
inline
Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v = Triangulation_2<Gt,Tds>::insert(p,lt,loc,li);
  restore_Delaunay(v);
  return(v);
}


template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
restore_Delaunay(Vertex_handle v)
{
  if(dimension() <= 1) return;

  Face_handle f=v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do {
    i = f->index(v);
    next = f->neighbor(ccw(i));  // turn ccw around v
    propagating_flip(f,i);
    f=next;
  } while(next != start);
  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove(Vertex_handle v )
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));
        
  if ( dimension() <= 1) Triangulation::remove(v);
  else  remove_2D(v);
        
  return;
}



template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
propagating_flip(Face_handle& f,int i)
{
  Face_handle n = f->neighbor(i);
      
  if ( ON_POSITIVE_SIDE != 
       side_of_oriented_circle(n,  f->vertex(i)->point()) ) {          
    return;
  }
  flip(f, i);
  propagating_flip(f,i);
  i = n->index(f->vertex(i));
  propagating_flip(n,i);
}

template <class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) {  _tds.remove_dim_down(&(*v));  }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole(v, hole);
    delete &(*v);
    set_number_of_vertices(number_of_vertices()-1);
  }
  return;       
}


template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
fill_hole(Vertex_handle v, std::list<Edge> & first_hole)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;
  
  Face_handle  f, ff, fn;
  int i =0,ii =0, in =0;
  Hole_list hole_list;
  Hole hole;
        
  hole_list.push_front(first_hole);
  
  while( ! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      Hole::iterator hit = hole.begin();
  
      // if the hole has only three edges, create the triangle
      if (hole.size() == 3) {
	hit = hole.begin();
	f = (*hit).first;        i = (*hit).second;
	ff = (* ++hit).first;    ii = (*hit).second;
	fn = (* ++hit).first;    in = (*hit).second;
	add_face(f,i,ff,ii,fn,in);
	continue;
      }
  
      // else find an edge with two finite vertices
      // on the hole boundary
      // and the new triangle adjacent to that edge
      //  cut the hole and push it back
  
      // first, ensure that a neighboring face
      // whose vertices on the hole boundary are finite
      // is the first of the hole
      bool finite= false;
      while (!finite){
	ff = (hole.front()).first;
	ii = (hole.front()).second;
	if ( is_infinite(ff->vertex(cw(ii))) ||
	     is_infinite(ff->vertex(ccw(ii)))) {
          hole.push_back(hole.front());
          hole.pop_front();
	}
	else finite=true;
      }
  
      // take the first neighboring face and pop it;
      ff = (hole.front()).first;
      ii =(hole.front()).second;
      hole.pop_front();
  
  
      Vertex_handle v0 = ff->vertex(ff->cw(ii)); Point p0 =v0->point();
      Vertex_handle v1 = ff->vertex(ff->ccw(ii)); Point p1 =v1->point();
      Vertex_handle v2 = infinite_vertex(); Point p2;
      Vertex_handle vv; Point p;
  
      Hole::iterator hdone = hole.end();
      hit =  hole.begin();
      Hole::iterator cut_after(hit);
  
      // if tested vertex is c with respect to the vertex opposite
      // to NULL neighbor,
      // stop at the before last face;
      hdone--;
      while( hit != hdone) {
	fn = (*hit).first;
	in = (*hit).second;
	vv = fn->vertex(ccw(in));
	if (is_infinite(vv)) {
	  if(is_infinite(v2)) cut_after = hit;
	}
	else {     // vv is a finite vertex
	  p = vv->point();
	  if (geom_traits().orientation(p0,p1,p) == COUNTERCLOCKWISE) {
	    if(is_infinite(v2)) { v2=vv; p2=p; cut_after=hit;}
	    else{
	      if( geom_traits().side_of_oriented_circle (p0,p1,p2,p) ==
		  ON_POSITIVE_SIDE){
		v2=vv; p2=p; cut_after=hit;}
	    }
	  }
	}
	++hit;
      }
  
  
      // create new triangle and update adjacency relations
      // Face_handle  newf = new Face(v0,v1,v2);
//       newf->set_neighbor(2,ff);
//       ff->set_neighbor(ii, newf);
      Face_handle newf;
  
  
      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes
  
      fn = (hole.front()).first;
      in = (hole.front()).second;
      if (fn->has_vertex(v2, i) && i == fn->ccw(in)) {
	//newf->set_neighbor(0,fn);
	//fn->set_neighbor(in,newf);
	newf = add_face(ff,ii,fn,in);
	hole.pop_front();
	hole.push_front(Edge( &(*newf),1));
	hole_list.push_front(hole);
      }
      else{
	fn = (hole.back()).first;
	in = (hole.back()).second;
	if (fn->has_vertex(v2, i) && i== fn->cw(in)) {
	  //newf->set_neighbor(1,fn);
	  //fn->set_neighbor(in,newf);
	  newf = add_face(fn,in,ff,ii);
	  hole.pop_back();
	  //hole.push_back(Edge(&(*newf),0));
	  hole.push_back(Edge(&(*newf),1));
	  hole_list.push_front(hole);
	}
	else{
	  // split the hole in two holes
	  newf = add_face(ff,ii,v2);
	  Hole new_hole;
	  ++cut_after;
	  while( hole.begin() != cut_after )
            {
              new_hole.push_back(hole.front());
              hole.pop_front();
            }
  
	  hole.push_front(Edge( &(*newf),1));
	  new_hole.push_front(Edge( &(*newf),0));
	  hole_list.push_front(hole);
	  hole_list.push_front(new_hole);
	}
      }
    }
}

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H
