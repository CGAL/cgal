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
// file          : include/CGAL/Triangulation_2.h
// package       : Triangulation (3.7)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_TRIANGULATION_2_H
#define CGAL_TRIANGULATION_2_H

#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>
#include <CGAL/Pointer.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_face_2.h>
#include <CGAL/Triangulation_vertex_2.h>
#include <CGAL/Triangulation_handles_2.h>
#include <CGAL/Triangulation_iterators_2.h>
#include <CGAL/Triangulation_circulators_2.h>
#include <CGAL/Triangulation_line_face_circulator_2.h>


CGAL_BEGIN_NAMESPACE


template < class Gt, class Tds>
class Triangulation_all_faces_iterator_2;

template < class Gt, class Tds>
class Triangulation_all_vertices__iterator_2;

template < class Gt, class Tds>
class Triangulation_all_edges_iterator_2;


template < class Gt, class Tds >
class Triangulation_2
  : public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
                (std::istream& is, Triangulation_2<Gt,Tds> &tr);
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
                (std::ostream& os, const Triangulation_2<Gt,Tds> &tr);
  friend Triangulation_face_iterator_2<Gt,Tds>;
  friend Triangulation_edge_iterator_2<Gt,Tds>;
  friend Triangulation_vertex_iterator_2<Gt,Tds>;

public:
  typedef Tds Triangulation_data_structure;
  typedef Triangulation_2<Gt,Tds> Triangulation;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;
 
  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;

  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef std::pair<Face_handle, int>                Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  //TODO for compatibility with previous version
  //add old names typedef iterators and functions

  typedef Triangulation_all_faces_iterator_2<Gt,Tds>    All_faces_iterator;
  typedef Triangulation_all_edges_iterator_2<Gt,Tds>    All_edges_iterator;
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> All_vertices_iterator;
 
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds>    
                                                    Finite_faces_iterator;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds>    
                                                    Finite_edges_iterator;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> 
                                                    Finite_vertices_iterator;

  typedef Triangulation_line_face_circulator_2<Gt,Tds>  Line_face_circulator;

  typedef Point value_type; // to have a back_inserter

  enum Locate_type {VERTEX=0, EDGE, FACE, OUTSIDE, COLLINEAR_OUTSIDE};


protected:
  Vertex_handle _infinite_vertex;
  Tds _tds;
  Gt  _gt;

public:
  // CONSTRUCTORS
  Triangulation_2(const Geom_traits& geom_traits=Geom_traits()); 
  Triangulation_2(const Triangulation_2<Gt,Tds> &tr);
 
  //Assignement
  Triangulation_2 &operator=(const Triangulation_2 &tr);

  //Helping
  void copy_triangulation(const Triangulation_2 &tr);
  void swap(Triangulation_2 &tr);
  void clear();


  //ACCESS FUNCTIONs
  int dimension() const { return _tds.dimension();}
  int number_of_vertices() const {return _tds.number_of_vertices() - 1;}
  const Geom_traits& geom_traits() const { return _gt;}
  int number_of_faces() const;
  Vertex_handle infinite_vertex() const;
  Vertex_handle finite_vertex() const;
  Face_handle infinite_face() const;
  
  //SETTING
  void set_number_of_vertices(int n) { _tds.set_number_of_vertices(n+1);}

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
 

  // TEST INFINITE FEATURES
  bool is_infinite(const Face_handle& f) const;
  bool is_infinite(const Vertex_handle& v) const; 
  bool is_infinite(const Face_handle& f, int i) const;
  bool is_infinite(const Edge& e) const;
  bool is_infinite(const Edge_circulator& ec) const;
  bool is_infinite(const Edge_iterator& ei) const;
 

 // GEOMETRIC FEATURES
  Triangle triangle(const Face_handle& f) const;
  Segment segment(const Face_handle& f, int i) const;
  Segment segment(const Edge& e) const;
  Segment segment(const Edge_circulator& ec) const;
  Segment segment(const Edge_iterator& ei) const;


  //INSERTION - DELETION - Flip
public:
  void   flip(Face_handle& f, int i);
  
  Vertex_handle insert_first(const Point& p);
  Vertex_handle insert_second(const Point& p);
  Vertex_handle insert_in_edge(const Point& p, Face_handle f,int i);
  Vertex_handle insert_outside_convex_hull(const Point& p, Face_handle f);
  Vertex_handle insert_outside_affine_hull(const Point& p);
  Vertex_handle insert(const Point& p,
		       Locate_type& lt,
		       Face_handle f = Face_handle() );
  Vertex_handle insert(const Point &p,
		       Face_handle f = Face_handle() );
  Vertex_handle push_back(const Point &p);
 
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
    template < class InputIterator >
    int insert(InputIterator first, InputIterator last)
    {
        int n = number_of_vertices();
        while(first != last){
            insert(*first);
            ++first;
        }
        return number_of_vertices() - n;
    }
#else
#if defined(LIST_H) || defined(__SGI_STL_LIST_H)
    int insert(std::list<Point>::const_iterator first,
               std::list<Point>::const_iterator last)
    {
        int n = number_of_vertices();
        while(first != last){
            insert(*first);
            ++first;
        }
        return number_of_vertices() - n;
    }
#endif // LIST_H
#if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
    int insert(std::vector<Point>::const_iterator first,
               std::vector<Point>::const_iterator last)
    {
        int n = number_of_vertices();
        while(first != last){
            insert(*first);
            ++first;
        }
        return number_of_vertices() - n;
    }
#endif // VECTOR_H
  //#ifdef ITERATOR_H
    int insert(std::istream_iterator<Point, ptrdiff_t> first,
               std::istream_iterator<Point, ptrdiff_t> last)
    {
        int n = number_of_vertices();
        while(first != last){
            insert(*first);
            ++first;
        }
        return number_of_vertices() - n;
    }
  //#endif // ITERATOR_H
    int insert(Point* first,
               Point* last)
    {
        int n = number_of_vertices();
        while(first != last){
            insert(*first);
            ++first;
        }
        return number_of_vertices() - n;
    }
#endif // TEMPLATE_MEMBER_FUNCTIONS
    

  void remove_degree_3(Vertex_handle  v, Face_handle f = Face_handle());
  void remove_first(Vertex_handle  v);
  void remove_second(Vertex_handle v);
  void remove(Vertex_handle  v);

  // POINT LOCATION
  Face_handle
  march_locate_1D(const Point& t, Locate_type& lt, int& li) const ;
  Face_handle
  march_locate_2D(const Face_handle& start,
                    const Point& t,
                    Locate_type& lt,
                    int& li) const;
  Face_handle
  locate(const Point& p,
           Locate_type& lt,
           int& li,
           Face_handle start = Face_handle()) const;

  Face_handle
  locate(const Point &p,
	 Face_handle start = Face_handle()) const;
    

  
  //TRAVERSING : ITERATORS AND CIRCULATORS
  Finite_faces_iterator finite_faces_begin() const;
  Finite_faces_iterator finite_faces_end() const;
  Finite_vertices_iterator finite_vertices_begin() const;
  Finite_vertices_iterator finite_vertices_end() const;
  Finite_edges_iterator finite_edges_begin() const;
  Finite_edges_iterator finite_edges_end() const; 

  All_faces_iterator all_faces_begin() const;
  All_faces_iterator all_faces_end() const;
  All_vertices_iterator all_vertices_begin() const;
  All_vertices_iterator all_vertices_end() const;
  All_edges_iterator all_edges_begin() const;
  All_edges_iterator all_edges_end() const; 

  Face_circulator incident_faces( const Vertex_handle& v) const;
  Face_circulator incident_faces( const Vertex_handle& v, 
				  const Face_handle& f) const;
  Vertex_circulator incident_vertices(const Vertex_handle& v) const;
  Vertex_circulator incident_vertices(const Vertex_handle& v,
				      const Face_handle& f) const;
  Edge_circulator incident_edges(const Vertex_handle& v) const;
  Edge_circulator incident_edges(const Vertex_handle& v,
				 const Face_handle& f) const;
 
  Line_face_circulator    line_walk(const Point& p,
				    const Point& q,
				    Face_handle f = Face_handle()) const;

// not documented

 // TO DEBUG
 void show_all();
 void show_face( typename Tds::Face_iterator fi);

 Oriented_side
 oriented_side(const Point &p0, const Point &p1,
	      const Point &p2, const Point &p) const;
    
 Bounded_side
 bounded_side(const Point &p0, const Point &p1,
	     const Point &p2, const Point &p) const;
    
 Oriented_side
 oriented_side(const Face_handle& f, const Point &p) const;

 bool 
 collinear_between(const Point& p, const Point& q, const Point& r) const;

protected:
  void remove_1D(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  
private:
Vertex_handle insert_outside_convex_hull_1(const Point& p, Face_handle f);
Vertex_handle insert_outside_convex_hull_2(const Point& p, Face_handle f);
void   make_hole ( Vertex_handle v, list<Edge> & hole);
void   fill_hole ( Vertex_handle v, list< Edge > & hole );
    

};

// CONSTRUCTORS
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::
Triangulation_2(const Geom_traits& geom_traits=Geom_traits()) 
     : _gt(geom_traits), _tds()
{
  _infinite_vertex = (Vertex *)_tds.insert_first();
}
  
// copy constructor duplicates vertices and faces
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::
Triangulation_2(const Triangulation_2<Gt,Tds> &tr)
    : _gt(tr._gt)
{
  _infinite_vertex = 
    (Vertex *) _tds.copy_tds(tr._tds, &(*tr.infinite_vertex()));
} 
 

  //Assignement
template <class Gt, class Tds >
Triangulation_2<Gt, Tds> &
Triangulation_2<Gt, Tds>::
operator=(const Triangulation_2 &tr)
{
  copy_triangulation(tr);
  return *this;
}

  // Helping functions

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::   
copy_triangulation(const Triangulation_2 &tr)
{
  clear();
  _gt = tr._gt;
  _infinite_vertex = 
  (Vertex *) _tds.copy_tds(tr._tds, &(*tr._infinite_vertex));
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>:: 
swap(Triangulation_2 &tr)
{
  Vertex_handle v= _infinite_vertex;
  _infinite_vertex = tr._infinite_vertex;
  tr._infinite_vertex = v;
  
  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t; 
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>:: 
clear()
{
  _tds.clear(); //detruit tous les sommets et toutes les faces
  _infinite_vertex = (Vertex *)_tds.insert_first();
}

template <class Gt, class Tds >
int
Triangulation_2<Gt, Tds>::
number_of_faces() const
{
  int count = _tds.number_of_faces();
  Face_circulator fc= infinite_vertex()->incident_faces(),
    done(fc);
  if ( ! fc.is_empty() ) {
    do { 
      --count; ++fc;
    }  while (fc != done);
  }
  return count;
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt,Tds>::Vertex_handle 
Triangulation_2<Gt,Tds>::
infinite_vertex() const
{
  return  _infinite_vertex;
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt,Tds>::Vertex_handle 
Triangulation_2<Gt,Tds>::
finite_vertex() const
{
  CGAL_triangulation_precondition (number_of_vertices() >= 1);
  return (finite_vertices_begin());
}
   
template <class Gt, class Tds >
inline
Triangulation_2<Gt,Tds>::Face_handle 
Triangulation_2<Gt,Tds>:
infinite_face() const
{
  return infinite_vertex()->face();
}


template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
is_valid(bool verbose = false, int level = 0) const
{
  bool result = _tds.is_valid(verbose, level);
  switch(dimension()) {
  case -1:
  case 0 : 
    break;
  case 1 :
    if (number_of_vertices() == 2) break;
    Finite_vertices_iterator it1 = finite_vertices_begin();
    Finite_vertices_iterator it2(it1); ++it2;
    Finite_vertices_iterator it3(it1); ++it3;
    while( it3 != finite_vertices_end()) {
     Orientation s = geom_traits().orientation(it1->vertex(0)->point(),
					       it2->vertex(1)->point(),
					       it3->vertex(2)->point()); 
     result = result && s == COLLINEAR ;
     CGAL_triangulation_assertion(result);
     result = result && collinear_between(it1->vertex(0)->point(),
					  it2->vertex(1)->point(),
					  it3->vertex(2)->point());
     CGAL_triangulation_assertion(result);
     ++it1 ; ++it2; ++it3:
    }    
    break;
  case 2 :
    for(Finite_faces_iterator it=finite_faces_begin(); 
	it!=faces_end(); it++){
      CGAL_triangulation_assertion( !is_infinite(it));
      Orientation s = geom_traits().orientation(it->vertex(0)->point(),
						it->vertex(1)->point(),
						it->vertex(2)->point());
      CGAL_triangulation_assertion( s == LEFTTURN );
      result = result && ( s == LEFTTURN );
    }

    Vertex_circulator start = infinite_vertex()->incident_vertices();
    Vertex_circulator pc(start);
    Vertex_circulator qc(start); ++qc;
    Vertex_circulator rc(start); ++rc; ++rc;
    do{
      Orientation s = geom_traits().orientation(pc->point(),
						qc->point(),
						rc->point());
      CGAL_triangulation_assertion( s != LEFTTURN );
      result = result && ( s != LEFTTURN );
      ++pc ; ++qc ; ++rc;
    }while(pc != start);
    break;
  default: result = false;
  }
  return result;
}


template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Face_handle& f) const 
{
  return f->has_vertex(infinite_vertex());
}


template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Vertex_handle& v) const 
{
  return v == infinite_vertex();
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Face_handle& f, int i) const 
{
  return is_infinite(f->vertex(ccw(i))) ||
         is_infinite(f->vertex(cw(i)));
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Edge& e) const 
{
  return is_infinite(e.first,e.second);
}
   
template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Edge_circulator& ec) const 
{
  return is_infinite(*ec);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Edge_iterator& ei) const 
{
  return is_infinite(*ei);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Triangle
Triangulation_2<Gt, Tds>::
triangle(const Face_handle& f) const
{
  CGAL_triangulation_precondition( ! is_infinite(f) );
  return Triangle(f->vertex(0)->point(),
		  f->vertex(1)->point(),
		  f->vertex(2)->point());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Face_handle& f, int i) const
{
  CGAL_triangulation_precondition( ! is_infinite(f,i));
  return Segment(f->vertex(ccw(i))->point(),
		 f->vertex(cw(i))->point());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge& e) const
{
  CGAL_triangulation_precondition(! is_infinite(e));
  return Segment(e.first->vertex(ccw(e.second))->point(),
		 e.first->vertex( cw(e.second))->point());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge_circulator& ec) const
{
  return segment(*ec);
}
    
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge_iterator& ei) const
{
  return segment(*ei);
}



template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
flip(Face_handle& f, int i)
{
  CGAL_triangulation_precondition ( ! f.is_null() );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2); 
    
  CGAL_triangulation_precondition( !is_infinite(f) && 
				   !is_infinite(f->neighbor(i)) );
  CGAL_triangulation_precondition( 
    geom_traits().orientation(f->vertex(i)->point(),
			      f->vertex(cw(i))->point(),
			      f->opposite_vertex(i)->point()) == RIGHTTURN &&
    geom_traits().orientation(f->vertex(i)->point(),
			      f->vertex(ccw(i))->point(),
			      f->opposite_vertex(i)->point()) == LEFTTURN); 
  _tds.flip( &(*f), i);
  return;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_first(const Point& p)
{
  CGAL_triangulation_precondition(number_of_vertices() == 0);
  Vertex_handle v = static_cast<Vertex*>(_tds.insert_second());
  v->set_point(p);
  return  v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_second(const Point& p)
{
  CGAL_triangulation_precondition(number_of_vertices() == 1);
   Vertex_handle v =  static_cast<Vertex*>
    (_tds.insert_dim_up(&(*infinite_vertex()),true));
   v->set_point(p);
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_in_edge(const Point& p, Face_handle f,int i)
{
  CGAL_triangulation_precondition(
      geom_traits().orientation(f->vertex(cw(i))->point(), p,
				f->vertex(ccw(i))->point() ) == COLLINEAR &&
      collinear_between(f->vertex(cw(i))->point(), p,
			f->vertex(ccw(i))->point()) );
  Vertex_handle v =   static_cast<Vertex*>
                     (_tds.insert_in_edge(&(*v), &(*f), i));
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition(is_infinite(f) && dimension() >= 1);
  Vertex_handle v;
  if (dimension() == 1)  v=insert_outside_convex_hull_1(p, f);
  else   v=insert_outside_convex_hull_2(v, f);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull_1(const Point& p, Face_handle f)
{
  int i = f->index(infinite_vertex());
  Face_handle  n = f->neighbor(i);
  int in = n->index(f);
  CGAL_triangulation_precondition( ! is_infinite(n));
  CGAL_triangulation_precondition(
	 geom_traits().orientation( n->vertex(in)->point(),
				    n->vertex(1-in)->point(),
				    v->point() ) == CGAL_COLLINEAR &&
	 collinear_between( n->vertex(in)->point(),
			    n->vertex(1-in)->point(),
			    v->point()) );
   Vertex_handle v=static_cast<Vertex*>(_tds.insert_in_edge(&(*v), &(*f), 2));
   v->set_point(p);
   return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull_2(const Point& p, Face_handle f)
{ 
  CGAL_triangulation_precondition(is_infinite(f));
  
  int li = f->index(infinite_vertex());
  Point q,r;
  q = f->vertex(ccw(li))->point();
  r = f->vertex(cw(li))->point();
  CGAL_triangulation_precondition( 
	    geom_traits().orientation(p,q,r) == CGAL_LEFTTURN);

  list<Face_handle> ccwlist;
  list<Face_handle> cwlist;
    
  Face_circulator fc = infinite_vertex()->incident_faces(f);
  bool done = false;
  while(! done) {
    fc--;
    li = fc->index(infinite_vertex());
    q = fc->vertex(ccw(li))->point();
    r = fc->vertex(cw(li))->point();
    if(geom_traits().orientation(p,q,r) == CGAL_LEFTTURN ) {
      ccwlist.push_back(&(*fc));
    }
    else {done=true;}
  }

  fc= infinite_vertex()->incident_faces(f);
  done = false;
  while(! done){
    fc++;
    li = fc->index(infinite_vertex());
    q = fc->vertex(ccw(li))->point();
    r = fc->vertex(cw(li))->point();
    if(geom_traits().orientation(p,q,r) == CGAL_LEFTTURN ) {
      cwlist.push_back(&(*fc));
    }
    else {done=true;}
  }

  Vertex_handle v =static_cast<Vertex*>(_tds.insert_in_face( &(*v), &(*f)));
  v->set_point(p);

  Face_handle fh;
  while ( ! ccwlist.empty()) {
    fh = ccwlist.front();
    li = ccw(fh->index(infinite_vertex()));
    _tds.flip( &(*fh) , li);
    ccwlist.pop_front();
  }

  while ( ! cwlist.empty()) {
    fh = cwlist.front();
    li = cw(fh->index(infinite_vertex()));
    _tds.flip( &(*fh) , li);
    cwlist.pop_front();
  }

  //reset infinite_vertex()->face();
  fc = v->incident_faces();
  while( ! is_infinite(&(*fc))) {
    fc++;}
  infinite_vertex()->set_face(&(*fc));

  return v;
} 

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_affine_hull(const Point& p)
{
  CGAL_triangulation_precondition(dimension() == 1);
  Face_handle f = (*edges_begin()).first;
  CGAL_Orientation or = geom_traits().orientation( f->vertex(0)->point(),
						   f->vertex(1)->point(),
						   p);
  CGAL_triangulation_precondition(or != CGAL_COLLINEAR);
  bool conform = ( or == CGAL_COUNTERCLOCKWISE);

  Vertex_handle v = static_cast<Vertex*>
    (_tds.insert_dim_up( &(*infinite_vertex()), conform));
  v->set_point(p);
  return v;
}
 
template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point& p, Locate_type& lt, Face_handle f)
{
  if(number_of_vertices() == 0) {
    lt = OUTSIDE_AFFINE_HULL;
    return(insert_first(p));
  }
  if(number_of_vertices() == 1) {
    if (geom_traits().compare(p,finite_vertex()->point()) ) {
      lt = VERTEX;
      return finite_vertex();
    }
    lt = OUTSIDE_AFFINE_HULL;
    return(insert_second(p));
  }

  int li;
  Face_handle loc = locate(p, lt, li, f);
  switch(lt){
  case FACE:
    return insert_in_face(p,loc);
    break;
  case EDGE:
    return insert_in_edge(p,loc,li);
    break;
  case OUTSIDE_CONVEX_HULL:
    return  insert_outside_convex_hull(p,loc);
    break;
  case OUTSIDE_AFFINE_HULL:
   return  insert_outside_affine_hull(p,v);
   break; 
  case VERTEX:
    return loc->vertex(li);
  default:
    CGAL_triangulation_assertion(false);  // locate step failed
  }
  return Vertex_handle();
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point &p, Face_handle f = Face_handle() );
{
  Locate_type lt;
  return insert(p, lt, f);
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  Locate_type lt;
  return insert(p, lt, NULL);
}

template <class Gt, class Tds >
void inline
Triangulation_2<Gt,Tds>::
remove_degree_3(Vertex_handle  v, Face_handle f = Face_handle())
{
  if (f == Face_handle()) f=v->face();
  _tds.remove_degree_3(&(*v), &(*f));
  return;
}

template <class Gt, class Tds >
void inline
Triangulation_2<Gt,Tds>::
remove_first(Vertex_handle  v)
{
  _tds.remove_second(&(*v));
  return;
}

template <class Gt, class Tds >
void inline
Triangulation_2<Gt,Tds>::
remove_second(Vertex_handle v)
{
  _tds.remove_dim_down(&(*v));
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::      
remove(Vertex_handle  v)
{
  CGAL_triangulation_precondition( ! v.is_null() );
  CGAL_triangulation_precondition( !is_infinite(v));
    
  if  (number_of_vertices() == 1)     remove_first(v);
  else if (number_of_vertices() == 2) remove_second(v);
  else   if ( dimension() == 1) remove_1D(v);
  else  remove_2D(v);
  return;
}

template <class Gt, class Tds >
void inline
Triangulation_2<Gt, Tds>::
remove_1D(Vertex_handle v)
{
  _tds.remove_1D(&(*v));
}


template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  //test the dimensionality of the resulting triangulation
  //it goes down to 1 iff
  // 1) any finite face is incident to v
  // 2) all vertices are colinear
  bool  dim1 = true; 
  Face_iterator fit = finite_faces_begin();
  while (dim1==true && fit != finite_faces_end()) {
    dim1 = dim1 && fit->has_vertex(v);
    fit++;
  }
  Face_circulator fic = v->incident_faces();
  while (is_infinite(fic)) {++fic;}
  Face_circulator done(fic);
  Face_handle start(fic); int iv = start->index(v);
  Point p = start->vertex(cw(iv))->point(); 
  Point q = start->vertex(ccw(iv))->point();
  while ( dim1 && ++fic != done) {
    iv = fic->index(v);
    if (fic->vertex(ccw(iv)) != infinite_vertex()) {
      dim1 = dim1 &&
	geom_traits().orientation(p, q, fic->vertex(ccw(iv))->point()) 
	== CGAL_COLLINEAR; 
    }
  }
	       
  if (dim1) { 
    _tds.remove_dim_down(&(*v));
  }
  else {
    list<Edge> hole;
    make_hole(v, hole);
    fill_hole(v, hole);
    delete &(*v);
    set_number_of_vertices(number_of_vertices()-1);
  }
  return;       
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
make_hole ( Vertex_handle v, list<Edge> & hole)
{
   
  list<Edge>::iterator hit;
  list<Face_handle> to_delete;

  Face_handle  f, ff, fn;
  int i =0,ii =0, in =0;
  Vertex_handle  vv;
      
  Face_circulator fc = v->incident_faces();
  Face_circulator done(fc);
  do {
    f = (*fc).handle(); fc++;
    i = f->index(v);
    fn = f->neighbor(i);
    in = fn->index(f);
    vv = f->vertex(cw(i));
    if( vv->face()==  f) vv->set_face(fn);
    vv = f->vertex(ccw(i));
    if( vv->face()== f) vv->set_face(fn);
    fn->set_neighbor(in, NULL);
    hole.push_back(Edge(fn,in));
    to_delete.push_back(f);
  }
  while(fc != done);

  while (! to_delete.empty()){
    (to_delete.front()).Delete();
    to_delete.pop_front();
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::						     
fill_hole ( Vertex_handle v, list< Edge > & hole )
  {
    typedef list< Edge > Hole;
	
   Point p = v->point();
   Face_handle  f, ff, fn;
   int i =0,ii =0, in =0; 
   Vertex_handle v0, v1, v2;
   Point   p0, p1, p2;
   Triangle t;
   CGAL_Bounded_side side;
   CGAL_Orientation or;
   int nhole;

   if( hole.size() != 3) {
     //  find the first edge  v0v1 on the hole boundary
     //  such that
     //  v0, v1 and the next vertex v2 are all finite
     //   v0v1v2 is a left turn and
     //  triangle v0v1v2 does not contain the removed point

     // if found create face v0v1v2
     // stop when no more faces can be created that way
                  
     nhole= hole.size(); 
     // nhole decount the number of hole edges passed
     // from the last created edges
     while (hole.size()>3 && nhole>0) {           
       //ff  = (Face *) ( (hole.front()).first);
       ff =  hole.front().first;
       ii  = (hole.front()).second;
       hole.pop_front();

       v0 = ff->vertex(cw(ii));
       v1 = ff->vertex(ccw(ii));
       if( !is_infinite(v0) && !is_infinite(v1)) {
	 //fn  = (Face *) ( (hole.front()).first);
	  fn =  hole.front().first;
	  in  = (hole.front()).second;
            
	  v2 = fn->vertex(ccw(in));
	  if( !is_infinite(v2)) {
	    p0 = v0->point();
	    p1 = v1->point();
	    p2 = v2->point();
	    or = geom_traits().orientation(p0,p1,p2);
	    if ( or  == CGAL_LEFTTURN) {
		  side = bounded_side(p0,p1,p2, p);
		  if( side == CGAL_ON_UNBOUNDED_SIDE || 
		      (side == CGAL_ON_BOUNDARY && collinear_between(p0, p, p2)) ) {
		    //create face
		    Face_handle  newf = new Face(v0,v1,v2);
		    newf->set_neighbor(2,ff);
		    newf->set_neighbor(0,fn);
		    ff->set_neighbor(ii, newf);
		    fn->set_neighbor(in,newf);
		    hole.pop_front();
		    //hole.push_front(Hole_neighbor(&(*newf),1));
		    hole.push_front(Edge(newf,1));
		    nhole = hole.size();
		    continue;
		  }
		}
	      }
	    }

	    // not possible to create face v0,v1,v2;
            //hole.push_back(Hole_neighbor(&(*ff),ii));
            hole.push_back( Edge (ff,ii));
	    nhole--;
	      // either the hole has only three edges
      // or all its finite vertices are reflex or flat
      // except may be one vertex whose corresponding ear 
      // includes the vertex being removed

      // deal with the last leftturn if any
      if(hole.size() != 3) {
	nhole = hole.size();
	while ( nhole>0) {
	  //ff = (Face *) ((hole.front()).first);
	  ff = ((hole.front()).first);
	  ii = (hole.front()).second;
	  hole.push_back(hole.front());
	  hole.pop_front();
	  nhole--;
	  
	  v0 = ff->vertex(cw(ii));
	  v1 = ff->vertex(ccw(ii));
	  if(is_infinite(v0) || is_infinite(v1))  continue;

	  //fn = (Face *) ((hole.front()).first);
	  fn = ((hole.front()).first);
	  in = (hole.front()).second;
	  v2 = fn->vertex(ccw(in));
          if( is_infinite(v2) ) continue;
	  p0 = v0->point();
	  p1 = v1->point();
	  p2 = v2->point();
	  CGAL_Orientation or = geom_traits().orientation(p0,p1,p2);
          if ( or  == CGAL_LEFTTURN) {
	    Face_handle  newf = new Face(v0,v1,v2);
	    newf->set_neighbor(2,ff);
	    newf->set_neighbor(0,fn);
	    ff->set_neighbor(ii, newf);
	    fn->set_neighbor(in,newf);
	    hole.pop_back();
	    hole.pop_front();
	    hole.push_front(Edge(newf,1));
	    break;
	  }
	}
      }


      if(hole.size() != 3) {
	// look for infinite vertex
	ff = (hole.front()).first;
	ii = (hole.front()).second;
	while ( ! is_infinite(ff->vertex(cw(ii)))){
	  hole.push_back(hole.front());
	  hole.pop_front();
	  //ff = (Face *)((hole.front()).first);
	  ff = (hole.front()).first;
	  ii = (hole.front()).second;
          }
	//create faces
          while(hole.size() != 3){
            //ff = (Face *)((hole.front()).first);
	    ff = (hole.front()).first;
	    //ff = ((hole.front()).first)->handle();
            ii = (hole.front()).second;
            hole.pop_front();
            //fn = (Face *)((hole.front()).first);
	    fn = (hole.front()).first;
	    //fn = ((hole.front()).first)->handle();
            in = (hole.front()).second;
            hole.pop_front();
            Face_handle  newf = new Face(infinite_vertex(),fn->vertex(cw(in)),
                                                     fn->vertex(ccw(in)));
            ff->set_neighbor(ii,newf);
            fn->set_neighbor(in,newf);
            newf->set_neighbor(0,fn);
            newf->set_neighbor(2,ff);
            //hole.push_front(Hole_neighbor(&(*newf),1));
	    hole.push_front(Edge(newf,1));
            }
      }
    
 // now hole has three edges
      list<Edge>::iterator hit;
      Face_handle  newf = new Face();
      hit = hole.begin();
      for(int j = 0;j<3;j++) {
	//ff = (Face *)((*hit).first);
	ff = (*hit).first;
	//ff = ((*hit).first)->handle();
	ii = (*hit).second;
	hit++;
	ff->set_neighbor(ii,newf);
	newf->set_neighbor(j,ff);
	newf->set_vertex(newf->ccw(j),ff->vertex(ff->cw(ii)));
      }
  }


  // POINT LOCATION

template <class Gt, class Tds >    
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_1D(const Point& t,
                    Locate_type& lt,
                    int& li) const
    
    {
        Face_handle ff = infinite_face();
	int iv = ff->index(infinite_vertex());
	Face_handle f = ff->neighbor(iv);
	CGAL_Orientation pqt = geom_traits().orientation(f->vertex(0)->point(), 
							 f->vertex(1)->point(),
							  t);
        if(pqt == CGAL_RIGHTTURN || pqt == CGAL_LEFTTURN) {
            lt = OUTSIDE_AFFINE_HULL;
	    return f;
        }
	
	Vertex_handle  u = f->vertex(0);
	Vertex_handle  v = f->vertex(1);
	if (collinear_between(t,u->point(), v->point())) {
	  lt = OUTSIDE_CONVEX_HULL;
	  li = 1;
	  return ff;
	}
	if(geom_traits().compare(t,u->point())){
	  lt = VERTEX;
	  li=0;
	  return f;
	}
	ff = ff->neighbor(1-iv); //the other infinite face
	f = ff->neighbor(ff->index(infinite_vertex()));
	u = f->vertex(0);
	v = f->vertex(1);
	if (collinear_between(u->point(), v->point(), t)) {
	  lt = OUTSIDE_CONVEX_HULL;
	  li = 0;
	  return ff;
	}
	
	Edge_iterator eit= edges_begin(); 
	for( ; eit != edges_end() ; eit++) {
	  u = (*eit).first->vertex(0);
	  v = (*eit).first->vertex(1);
	  if(geom_traits().compare(t,v->point())){
	    lt = VERTEX;
	    li = 1;
	    return (*eit).first;
	  }
	  if(collinear_between(u->point(), t, v->point())){
	    lt = EDGE;
	    li =  2;
	    return (*eit).first;
	  }
	}
	CGAL_triangulation_assertion(false);
    }

template <class Gt, class Tds >    
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_2D(const Face_handle& start,
		const Point& t,
		Locate_type& lt,
		int& li) const
{
     //    CGAL_triangulation_precondition( ! is_infinite(start) );
        Triangulation_2 *ncthis = (Triangulation_2 *)this;
    
        Point p(start->vertex(0)->point());
        if(geom_traits().compare_x(t,p) == EQUAL &&  
	   geom_traits().compare_y(t,p) == EQUAL){
            lt = VERTEX;
            li = 0;
            return start;
        }

        Line_face_circulator lfc(start->vertex(0),
                                 ncthis,
                                 t);
	
        if(lfc.collinear_outside()) {
            // point t lies outside or on the convex hull
            // we walk clockwise on the hull to decide
            int i = lfc->index(infinite_vertex());
            p = lfc->vertex(ccw(i))->point();
            if(geom_traits().compare_x(t,p) == EQUAL &&  
	       geom_traits().compare_y(t,p) == EQUAL){
                lt = VERTEX;
                li = ccw(i);
                return lfc;
            }
         Point q(lfc->vertex(cw(i))->point());
	 Orientation pqt;
         Face_handle f(lfc);
         while(1){
           if(geom_traits().compare_x(t,q) == EQUAL &&  
	      geom_traits().compare_y(t,q) == EQUAL){
             lt = VERTEX;
             li = cw(i);
             return f;
           }
	   pqt = geom_traits().orientation(p,q,t);
	   if (pqt == COLLINEAR && collinear_between(p, t, q)){
	     lt = EDGE;
	     li = i;
	     return f;
	   }
	   if (pqt == CGAL_LEFTTURN){
	     lt = OUTSIDE_CONVEX_HULL;
	     return f ;
	   }
	   	       
           // go to the next face
           f = f->neighbor(ccw(i));
           i = f->index(infinite_vertex());
           p = q;
           q = f->vertex(cw(i))->point();
     	 }
        }

        while(! lfc.locate(t, lt, li) ){
	  lfc.increment();
        }
        return lfc;
    }

    
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt,Tds>::
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle start = Face_handle()) const
{
  if( dimension() == 0) {
    if(number_of_vertices() == 0) {
      lt = OUTSIDE_AFFINE_HULL;
    } else { // number_of_vertices() == 1
      lt = geom_traits().compare(p,finite_vertex()->point()) ? 
	VERTEX : OUTSIDE_AFFINE_HULL;
    }
    return NULL;
  }
  if(dimension() == 1){
    return march_locate_1D(p, lt, li);
  }
    
  if(start.is_null()){
    start = infinite_face()->
      neighbor(infinite_face()->index(infinite_vertex()));
  }else if(is_infinite(start)){
    start = start->neighbor(start->index(infinite_vertex()));
  }
  return march_locate_2D(start, p, lt, li);
}


template <class Gt, class Tds >
Triangulation_2<Gt, Tds>:: Face_handle
Triangulation_2<Gt, Tds>::
locate(const Point &p,
       Face_handle start = Face_handle()) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, start);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_faces_iterator
Triangulation_2<Gt, Tds>::
finite_faces_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds> *)this;
  return Finite_faces_iterator(ncthis); 
} 

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_faces_iterator
Triangulation_2<Gt, Tds>::
finite_faces_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return Finite_faces_iterator(ncthis,1);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_vertices_iterator
Triangulation_2<Gt, Tds>::
finite_vertices_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return Finite_vertices_iterator(ncthis);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_vertices_iterator
Triangulation_2<Gt, Tds>::
finite_vertices_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return Finite_vertices_iterator(ncthis,1);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_edges_iterator
Triangulation_2<Gt, Tds>::
finite_edges_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return Finite_edges_iterator(ncthis);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Finite_edges_iterator
Triangulation_2<Gt, Tds>::
finite_edges_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return Finite_edges_iterator(ncthis,1);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_faces_iterator
Triangulation_2<Gt, Tds>::
all_faces_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds> *)this;
  return All_faces_iterator(ncthis); 
} 

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_faces_iterator
Triangulation_2<Gt, Tds>::
all_faces_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return All_faces_iterator(ncthis,1);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_vertices_iterator
Triangulation_2<Gt, Tds>::
all_vertices_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return All_vertices_iterator(ncthis);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_vertices_iterator
Triangulation_2<Gt, Tds>::
all_vertices_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return All_vertices_iterator(ncthis,1);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_edges_iterator
Triangulation_2<Gt, Tds>::
all_edges_begin() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return All_edges_iterator(ncthis);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::All_edges_iterator
Triangulation_2<Gt, Tds>::
all_edges_end() const
{
  Triangulation_2<Gt, Tds>* 
    ncthis = (Triangulation_2<Gt, Tds>*)this;
  return All_edges_iterator(ncthis,1);
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Face_circulator
Triangulation_2<Gt, Tds>::
incident_faces( const Vertex_handle& v) const
{
  return v->incident_faces();
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Face_circulator
Triangulation_2<Gt, Tds>::
incident_faces( const Vertex_handle& v, const Face_handle& f) const
{
  return v->incident_faces(f);
}  

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Vertex_circulator
Triangulation_2<Gt, Tds>::
incident_vertices(const Vertex_handle& v) const
{
  return v->incident_vertices();
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Vertex_circulator
Triangulation_2<Gt, Tds>::  
incident_vertices(const Vertex_handle& v,
		  const Face_handle& f) const
{
  return v->incident_vertices(f);
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Edge_circulator
Triangulation_2<Gt, Tds>::  
incident_edges(const Vertex_handle& v) const
{
  return v->incident_edges();
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Edge_circulator
Triangulation_2<Gt, Tds>::    
incident_edges(const Vertex_handle& v,
	       const Face_handle& f) const
{
  return v->incident_edges(f);
}





template <class Gt, class Tds >
Triangulation_2<Gt, Tds>:: Line_face_circulator  
Triangulation_2<Gt, Tds>::    
line_walk(const Point& p,
              const Point& q,
              Face_handle f) 
{
        CGAL_triangulation_precondition( (dimension() == 2) && 
					 (! geom_traits().compare(p,q)) );
    
        Line_face_circulator lfc = (f.is_null())
                                   ? Line_face_circulator(p, q, this)
                                   : Line_face_circulator(p, q, f, this);
    
        if( (!lfc.is_empty()) && is_infinite( lfc )){
            return Line_face_circulator();
        }
        return lfc;
}
   
template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
oriented_side(const Point &p0, const Point &p1,
                  const Point &p2, const Point &p) const
{
        // depends on the orientation of the vertices
        Orientation o1 = geom_traits().orientation(p0, p1, p),
                         o2 = geom_traits().orientation(p1, p2, p),
                          o3 = geom_traits().orientation(p2, p0, p),
                         ot = geom_traits().orientation(p0, p1, p2);

        if (o1 == COLLINEAR ||
            o2 == COLLINEAR ||
            o3 == COLLINEAR)  {
            if ((o1 == COLLINEAR &&
                 collinear_between(p0, p, p1)) ||
                (o2 == COLLINEAR &&
                 collinear_between(p1, p, p2)) ||
                (o3 == COLLINEAR &&
                 collinear_between(p2, p, p0)))
                {
                return  ON_ORIENTED_BOUNDARY;
            }
                // for ot == ON_ORIENTED_BOUNDARY the following also
                // does the right thing:
                return (ot == LEFTTURN) ? ON_POSITIVE_SIDE
                                             : ON_NEGATIVE_SIDE;
            }
        if (ot == RIGHTTURN)
            {
                if (o1 == RIGHTTURN &&
                    o2 == RIGHTTURN &&
                    o3 == RIGHTTURN)
                    {
                        return ON_NEGATIVE_SIDE;
                    }
                return ON_POSITIVE_SIDE;
            }
        if (o1 == LEFTTURN &&
            o2 == LEFTTURN &&
            o3 == LEFTTURN)
            {
                return ON_POSITIVE_SIDE;
            }
        return ON_NEGATIVE_SIDE;
    }

template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{
        Comparison_result c_pr = geom_traits().compare_x(p, r);
        Comparison_result c_pq;
        Comparison_result c_qr;
        if(c_pr == EQUAL) {
            c_pr = geom_traits().compare_y(p, r);
            c_pq = geom_traits().compare_y(p, q);
            c_qr = geom_traits().compare_y(q, r);
        } else {
            c_pq = geom_traits().compare_x(p, q);
            c_qr = geom_traits().compare_x(q, r);
        }
        return ( (c_pq == SMALLER) && (c_qr == SMALLER) ) ||
            ( (c_qr == LARGER) && (c_pq == LARGER) );
    
}

template <class Gt, class Tds >
Bounded_side
Triangulation_2<Gt, Tds>::
bounded_side(const Point &p0, const Point &p1,
                 const Point &p2, const Point &p) const
{
      Orientation o1 = geom_traits().orientation(p0, p1, p),
        o2 = geom_traits().orientation(p1, p2, p),
        o3 = geom_traits().orientation(p2, p0, p),
        ot = geom_traits().orientation(p0, p1, p2);

      if(o1 == COLLINEAR ||
         o2 == COLLINEAR ||
         o3 == COLLINEAR)
        {
          if ((o1 == COLLINEAR &&
               collinear_between(p0, p, p1)) ||
              (o2 == COLLINEAR &&
               collinear_between(p1, p, p2)) ||
              (o3 == COLLINEAR &&
               collinear_between(p2, p, p0)))
            {
              return  ON_BOUNDARY;
            }
          return ON_UNBOUNDED_SIDE;
        }
      if (ot == RIGHTTURN)
        {
          if(o1==RIGHTTURN &&
             o2==RIGHTTURN &&
             o3==RIGHTTURN)
            {
              return ON_BOUNDED_SIDE;
            }
          return ON_UNBOUNDED_SIDE;
        }
      if (o1 == LEFTTURN &&
          o2 == LEFTTURN &&
          o3 == LEFTTURN)
        {
          return ON_BOUNDED_SIDE;
        }
      return ON_UNBOUNDED_SIDE;
    }

template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
    {
        Comparison_result c_pr = geom_traits().compare_x(p, r);
        Comparison_result c_pq;
        Comparison_result c_qr;
        if(c_pr == EQUAL) {
            c_pr = geom_traits().compare_y(p, r);
            c_pq = geom_traits().compare_y(p, q);
            c_qr = geom_traits().compare_y(q, r);
        } else {
            c_pq = geom_traits().compare_x(p, q);
            c_qr = geom_traits().compare_x(q, r);
        }
        return ( (c_pq == SMALLER) && (c_qr == SMALLER) ) ||
            ( (c_qr == LARGER) && (c_pq == LARGER) );
    
    }


template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
oriented_side(const Face_handle& f, const Point &p) const
    {
      CGAL_triangulation_precondition ( dimension()==2); 
      return oriented_side(f->vertex(0)->point(),
                             f->vertex(1)->point(),
                             f->vertex(2)->point(),
                             p);
    }
 
template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_all()
{
  cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<endl;
  typename Tds::Face_iterator fi = _tds.faces_begin();
  cerr<<"***"<<endl;
  while(fi != _tds.faces_end()) {
    show_face(fi);
    ++fi;
  }
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_face( typename Tds::Face_iterator fi)
{
  cerr << "face : "<<(void*)&(*fi)<<" => "<<endl;
  int i = fi->dimension(); 
  switch(i){
  case 0:
    cerr <<"point :"<<(fi->vertex(0)->point())<<" / voisin "<<&(*(fi->neighbor(0)))
	 <<"["<<(fi->neighbor(0))->vertex(0)->point()<<"]"
      	<<endl;
    break;
  case 1:
     cerr <<"point :"<<(fi->vertex(0)->point())<<" / voisin "<<&(*(fi->neighbor(0)))
				<<"["<<(fi->neighbor(0))->vertex(0)->point()
				<<"/"<<(fi->neighbor(0))->vertex(1)->point()<<"]"
  	<<endl;
     cerr <<"point :"<<(fi->vertex(1)->point())<<" / voisin "<<&(*(fi->neighbor(1)))
				<<"["<<(fi->neighbor(1))->vertex(0)->point()
				<<"/"<<(fi->neighbor(1))->vertex(1)->point()<<"]"
				<<endl;
     break;
  case 2:
  cerr <<"point :"<<(fi->vertex(0)->point())<<" / voisin "<<&(*(fi->neighbor(0)))
				<<"["<<(fi->neighbor(0))->vertex(0)->point()
				<<"/"<<(fi->neighbor(0))->vertex(1)->point()
				<<"/"<<(fi->neighbor(0))->vertex(2)->point()<<"]"
					<<endl;
  cerr <<"point :"<<(fi->vertex(1)->point())<<" / voisin "<<&(*(fi->neighbor(1)))
				<<"["<<(fi->neighbor(1))->vertex(0)->point()
				<<"/"<<(fi->neighbor(1))->vertex(1)->point()
				<<"/"<<(fi->neighbor(1))->vertex(2)->point()<<"]"
				<<endl;
  cerr <<"point :"<<(fi->vertex(2)->point())<<" / voisin "<<&(*(fi->neighbor(2)))
				<<"["<<(fi->neighbor(2))->vertex(0)->point()
				<<"/"<<(fi->neighbor(2))->vertex(1)->point()
				<<"/"<<(fi->neighbor(2))->vertex(2)->point()<<"]"
				<<endl;
  }
  return;
}


template <class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Triangulation_2<Gt, Tds> &tr)
{
  // to debug
  //operator<<(os, tr._tds);
  //os << tr._tds;

  std::map< void*, int, std::less<void*> > V;
  std::map< void*, int, std::less<void*> > F;
  typename Triangulation_2<Gt, Tds>::Vertex_handle  v;

    int n = tr.number_of_vertices() + 1;
    int m = tr.number_of_faces();
    if(is_ascii(os)){
        os << n << ' ' << m << ' ' << tr.dimension() << endl;
    } else {
        os << n << m << tr.dimension();
    }

    // write the vertex at infinity
    int i = 0;
    v = tr.infinite_vertex();
    V[&(*v)] = i;
    os << v->point();
    if(is_ascii(os)){
        os << ' ';
    }
    if(n == 1){
        return os;
    }

    // write the finite vertices
    {
        typename Triangulation_2<Gt, Tds>::Vertex_iterator
          it = tr.vertices_begin();

        while(it != tr.vertices_end()){
            V[&(*it)] = ++i;
            os << it->point();
            if(is_ascii(os)){
                os << ' ';
            }
            ++it;
        }
    }
    CGAL_triangulation_assertion( (i+1) == n );
    if(is_ascii(os)){ os << "\n";}

    if(n == 2){
        return os;
    }

    i = 0;
    // vertices of the finite faces
    {
        typename Triangulation_2<Gt, Tds>::Face_iterator
          it = tr.faces_begin();

        while(it != tr.faces_end()){
	  F[&(*it)] = i++;
	  for(int j = 0; j < 3; j++){
	    os << V[&(*it->vertex(j))];
	    if(is_ascii(os)){
	      if(j==2) {
		os << "\n";
	      } else {
		os <<  ' ';
	      }
	    }
	  }
	  ++it;
        }
    }

    // vertices of the infinite faces
    {
        typename Triangulation_2<Gt, Tds>::Face_circulator
            fc = tr.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            F[&(*fc)] = i++;
            for(int j = 0; j < 3; j++){
                os << V[&(*fc->vertex(j))];
                if(is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
        }while(++fc != done);
    }
    CGAL_triangulation_assertion( i == m );

    // neighbor pointers of the finite faces
    {
        typename Triangulation_2<Gt, Tds>::Face_iterator
            it = tr.faces_begin();
        while(it != tr.faces_end()){
            for(int j = 0; j < 3; j++){
                os << F[&(*it->neighbor(j))];
                if(is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
            ++it;
        }
    }

    // neighbor pointers of the infinite faces
    {
        typename Triangulation_2<Gt, Tds>::Face_circulator
            fc = tr.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            for(int j = 0; j < 3; j++){
                os << F[&(*fc->neighbor(j))];
                if(is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
        }while(++fc != done);
    }

    
    return os ;
}



template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_2<Gt, Tds> &tr)
{
  
  return operator>>(is, tr._tds);
}
 
CGAL_END_NAMESPACE
    

#endif //CGAL_TRIANGULATION_2_H

