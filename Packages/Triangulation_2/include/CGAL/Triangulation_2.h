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
// release_date  : 
//
// file          : include/CGAL/Triangulation_2.h
// package       : Triangulation 
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
#include <CGAL/circulator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Triangulation_face_2.h>
#include <CGAL/Triangulation_vertex_2.h>
#include <CGAL/Triangulation_handles_2.h>
#include <CGAL/Triangulation_iterators_2.h>
#include <CGAL/Triangulation_circulators_2.h>
#include <CGAL/Triangulation_line_face_circulator_2.h>



CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds > class Triangulation_2;
template < class Gt, class Tds > std::istream& operator>>
    (std::istream& is, Triangulation_2<Gt,Tds> &tr);
template < class Gt, class Tds >  std::ostream& operator<<
  (std::ostream& os, const Triangulation_2<Gt,Tds> &tr);
  

template < class Gt, 
           class Tds = Triangulation_data_structure_using_list_2 <
                       Triangulation_vertex_base_2<Gt>,
		       Triangulation_face_base_2<Gt> > >
class Triangulation_2
  : public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
                (std::istream& is, Triangulation_2<Gt,Tds> &tr);
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
                (std::ostream& os, const Triangulation_2<Gt,Tds> &tr);
  friend class Triangulation_all_faces_iterator_2<Gt,Tds>;
  friend class Triangulation_all_edges_iterator_2<Gt,Tds>;
  friend class Triangulation_all_vertices_iterator_2<Gt,Tds>;
  friend class Triangulation_finite_faces_iterator_2<Gt,Tds>;
  friend class Triangulation_finite_edges_iterator_2<Gt,Tds>;
  friend class Triangulation_finite_vertices_iterator_2<Gt,Tds>;

public:
  typedef Tds Triangulation_data_structure;
  typedef Triangulation_2<Gt,Tds> Triangulation;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point_2       Point;
  typedef typename Geom_traits::Segment_2     Segment;
  typedef typename Geom_traits::Triangle_2    Triangle;
  typedef typename Geom_traits::Orientation_2 Orientation_2;
  typedef typename Geom_traits::Compare_x_2   Compare_x;
  typedef typename Geom_traits::Compare_y_2   Compare_y;
 
  typedef Triangulation_face_2<Gt,Tds> Face;
  typedef Triangulation_vertex_2<Gt,Tds> Vertex;

  typedef Triangulation_face_handle_2<Gt,Tds> Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef std::pair<Face_handle, int>                Edge;

  typedef Triangulation_face_circulator_2<Gt,Tds>      Face_circulator;
  typedef Triangulation_edge_circulator_2<Gt,Tds>      Edge_circulator;
  typedef Triangulation_vertex_circulator_2<Gt,Tds>    Vertex_circulator;

  typedef Triangulation_all_faces_iterator_2<Gt,Tds>    All_faces_iterator;
  typedef Triangulation_all_edges_iterator_2<Gt,Tds>    All_edges_iterator;
  typedef Triangulation_all_vertices_iterator_2<Gt,Tds> All_vertices_iterator;
 
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds>    
                                                    Finite_faces_iterator;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds>    
                                                    Finite_edges_iterator;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> 
                                                    Finite_vertices_iterator;

  //for compatibility with previous version
  typedef Triangulation_finite_faces_iterator_2<Gt,Tds>    Face_iterator;
  typedef Triangulation_finite_edges_iterator_2<Gt,Tds>    Edge_iterator;
  typedef Triangulation_finite_vertices_iterator_2<Gt,Tds> Vertex_iterator;

  typedef Triangulation_line_face_circulator_2<Gt,Tds>  Line_face_circulator;

  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  typedef Project_point<Vertex>                           Proj_point;
  typedef Iterator_project<Vertex_iterator, 
                           Proj_point,
	                   const Point&, 
                           const Point*,
                           std::ptrdiff_t,
                           std::bidirectional_iterator_tag>  Point_iterator;

  typedef Point value_type; // to have a back_inserter
  typedef const value_type&    const_reference; 
  

  enum Locate_type {VERTEX=0, EDGE, FACE, 
		    OUTSIDE_CONVEX_HULL,
		    OUTSIDE_AFFINE_HULL};

protected:
//  //Helping classes

// The following nested template classes do not compile on bcc
//  //  to be used as adaptators from iterators with Edge value_type
//  //   to an iterator with Tds::Edge as value type
//   template<class It>
//   class To_tds_edge_iterator : public It {
//   public:
//     typedef typename Triangulation_data_structure::Edge  Tds_Edge;
//     To_tds_edge_iterator() {}
//     To_tds_edge_iterator(It i) : It(i) {} 
//     Tds_Edge  operator*() {
//       Edge e = It::operator*();
//       return Tds_Edge( &*(e.first), e.second);
//     }
//   };

//  //  to be used as adaptators from iterators with Face_jhandle value_type
//  //   to an iterator with Tds::Face* as value type
//   template<class It>
//   class To_tds_face_iterator : public It {
//   public:
//     typedef typename Triangulation_data_structure::Face  Tds_Face;
//     To_tds_face_iterator() {}
//     To_tds_face_iterator(It i) : It(i) {} 
//     Tds_Face* operator*() { return  &*(It::operator*() ); }
//   };

protected:
  Gt  _gt;
  Tds _tds;
  Vertex_handle _infinite_vertex;

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
  void set_number_of_vertices(int n) {_tds.set_number_of_vertices(n+1);}
  void set_infinite_vertex(const Vertex_handle& v) {_infinite_vertex=v;}
  void set_dimension(int n) {_tds.set_dimension(n);}

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
 

  // TEST INFINITE FEATURES AND OTHER FEATURES
  bool is_infinite(const Face_handle& f) const;
  bool is_infinite(const Vertex_handle& v) const; 
  bool is_infinite(const Face_handle& f, int i) const;
  bool is_infinite(const Edge& e) const;
  bool is_infinite(const Edge_circulator& ec) const;
  bool is_infinite(const All_edges_iterator& ei) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr,
	       int & i) const;
  bool includes_edge(Vertex_handle va, Vertex_handle vb,
		     Vertex_handle& vbb,
		     Face_handle& fr, int & i) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
       Face_handle &fr) const;

 // GEOMETRIC FEATURES AND CONSTRUCTION
  Triangle triangle(const Face_handle& f) const;
  Segment segment(const Face_handle& f, int i) const;
  Segment segment(const Edge& e) const;
  Segment segment(const Edge_circulator& ec) const;
  Segment segment(const All_edges_iterator& ei) const;
  Segment segment(const Finite_edges_iterator& ei) const;
  Point circumcenter(Face_handle  f) const; 
  Point circumcenter(const Point& p0, 
		     const Point& p1, 
		     const Point& p2) const;
  

  //INSERTION - DELETION - Flip
public:
  void   flip(Face_handle f, int i);
  
  Vertex_handle insert_first(const Point& p);
  Vertex_handle insert_second(const Point& p);
  Vertex_handle insert_in_edge(const Point& p, Face_handle f,int i);
  Vertex_handle insert_in_face(const Point& p, Face_handle f);
  Vertex_handle insert_outside_convex_hull(const Point& p, Face_handle f);
  Vertex_handle insert_outside_affine_hull(const Point& p);
  Vertex_handle insert(const Point &p, Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
//   template < class InputIterator >
//   int insert(InputIterator first, InputIterator last);
  Vertex_handle push_back(const Point& a);
 
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
  Point_iterator points_begin() const;
  Point_iterator points_end() const;

  All_faces_iterator all_faces_begin() const;
  All_faces_iterator all_faces_end() const;
  All_vertices_iterator all_vertices_begin() const;
  All_vertices_iterator all_vertices_end() const;
  All_edges_iterator all_edges_begin() const;
  All_edges_iterator all_edges_end() const; 

  //for compatibility with previous versions
  Face_iterator faces_begin() const {return finite_faces_begin();}
  Face_iterator faces_end() const {return finite_faces_end();}
  Edge_iterator edges_begin() const {return finite_edges_begin();}
  Edge_iterator edges_end() const {return finite_edges_end();}
  Vertex_iterator vertices_begin() const {return finite_vertices_begin();}
  Vertex_iterator vertices_end() const {return finite_vertices_end();}

  Face_circulator incident_faces( Vertex_handle v, 
				  Face_handle f = Face_handle()) const;
  Vertex_circulator incident_vertices(Vertex_handle v,
				      Face_handle f = Face_handle()) const;
  Edge_circulator incident_edges(Vertex_handle v,
				 Face_handle f = Face_handle()) const;
 
  Line_face_circulator    line_walk(const Point& p,
				    const Point& q,
				    Face_handle f = Face_handle()) const;

 // TO DEBUG
 void show_all();
 void show_face( Face_handle fh);

  // IO
// template < class Stream >
// Stream&  draw_triangulation(Stream& os) const;

 //PREDICATES
 Oriented_side
 oriented_side(const Point &p0, const Point &p1,
	       const Point &p2, const Point &p) const;
    
 Bounded_side
 bounded_side(const Point &p0, const Point &p1,
	      const Point &p2, const Point &p) const;
    
 Oriented_side
 oriented_side(const Face_handle& f, const Point &p) const;

 Oriented_side
 side_of_oriented_circle(Face_handle f, const Point & p) const; 

 bool 
 collinear_between(const Point& p, const Point& q, const Point& r)
   const;

  Comparison_result compare_x(const Point& p, const Point& q) const;
  Comparison_result compare_y(const Point& p, const Point& q) const;
  bool               xy_equal(const Point& p, const Point& q) const;
  Orientation orientation(const Point& p, 
			  const Point& q, 
			  const Point& r) const;


protected:
  void remove_1D(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  bool test_dim_down(Vertex_handle v);
  void fill_hole(Vertex_handle v, std::list<Edge> & hole);
  void fill_hole_delaunay(std::list<Edge> & hole);

public:
  void make_hole(Vertex_handle v, std::list<Edge> & hole);
//   template<class EdgeIt>
//   Vertex_handle star_hole( Point p, 
// 			      EdgeIt edge_begin,
// 			      EdgeIt edge_end);

//   template<class EdgeIt, class FaceIt>
//   Vertex_handle star_hole( Point p, 
// 			      EdgeIt edge_begin,
// 			      EdgeIt edge_end,
// 			      FaceIt face_begin,
// 			      FaceIt face_end);

  Face_handle create_face(Face_handle f1, int i1,
			  Face_handle f2, int i2,
			  Face_handle f3, int i3);
  Face_handle create_face(Face_handle f1, int i1,
			  Face_handle f2, int i2);
  Face_handle create_face(Face_handle f, int i, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3,
			  Face_handle f1, Face_handle f2, Face_handle f3);
  Face_handle create_face();
  Face_handle create_face(Face_handle); //calls copy constructor of Face
  void delete_face(Face_handle f);

  Vertex_handle file_input(std::istream& is);
  void file_output(std::ostream& os) const;

private:
  Vertex_handle insert_outside_convex_hull_1(const Point& p, Face_handle f);
  Vertex_handle insert_outside_convex_hull_2(const Point& p, Face_handle f);
  
  // template members
public:
template < class Stream >
Stream&  draw_triangulation(Stream& os) const 
{
  Finite_edges_iterator it = finite_edges_begin();
  for( ;it != finite_edges_end() ; ++it) {
    os << segment(it);
  }
  return os;
}

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

public:
  template<class EdgeIt>
  Vertex_handle star_hole( Point p, 
			   EdgeIt edge_begin,
			   EdgeIt edge_end) {
    std::list<Face_handle> empty_list;
    return star_hole(p, 
		     edge_begin, 
		     edge_end, 
		     empty_list.begin(),
		     empty_list.end());
  }

  template<class EdgeIt, class FaceIt>
  Vertex_handle star_hole( Point p, 
			   EdgeIt edge_begin,
			   EdgeIt edge_end,
			   FaceIt face_begin,
			   FaceIt face_end) {
    typedef typename Triangulation_data_structure::Edge  Tds_Edge;
    typedef typename Triangulation_data_structure::Face  Tds_Face;
    typedef To_tds_edge_iterator<EdgeIt, Tds_Edge> Tds_ei;
    typedef To_tds_face_iterator<FaceIt, Tds_Face> Tds_fi;
    Vertex_handle v = static_cast<Vertex*> 
                       (_tds.star_hole( Tds_ei(edge_begin), 
					Tds_ei(edge_end),
					Tds_fi(face_begin),
					Tds_fi(face_end)) );
    v->set_point(p);
    return v;
  }

};

// CONSTRUCTORS
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::
Triangulation_2(const Geom_traits& geom_traits) 
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
  _tds.clear();
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
Triangulation_2<Gt,Tds>::
infinite_face() const
{
  return infinite_vertex()->face();
}


template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  if (dimension() <= 0 ||
      (dimension()==1 && number_of_vertices() == 2 ) ) return result;

    if (dimension() == 1) {
    Finite_vertices_iterator it1 = finite_vertices_begin(),
                             it2(it1), it3(it1);
    ++it2;
    ++it3; ++it3;
    while( it3 != finite_vertices_end()) {
     Orientation s = orientation(it1->point(),
				 it2->point(),
				 it3->point()); 
     result = result && s == COLLINEAR ;
     CGAL_triangulation_assertion(result);
     ++it1 ; ++it2; ++it3;
    }
  }    

  else { //dimension() == 2
    for(Finite_faces_iterator it=finite_faces_begin(); 
	it!=finite_faces_end(); it++) {
      CGAL_triangulation_assertion( ! is_infinite(it));
      Orientation s = orientation(it->vertex(0)->point(),
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
      Orientation s = orientation(pc->point(),
				  qc->point(),
				  rc->point());
      CGAL_triangulation_assertion( s != LEFTTURN );
      result = result && ( s != LEFTTURN );
      ++pc ; ++qc ; ++rc;
    }while(pc != start);
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
is_infinite(const All_edges_iterator& ei) const 
{
  return is_infinite(*ei);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb) const
{
  return _tds.is_edge( &(*(va)), &(*(vb)));
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i) const
{
  typename Tds::Face* f ;
  bool b = _tds.is_edge( &(*(va)), &(*(vb)), f, i);
  fr = Face_handle(static_cast<Face*>(f));
  return b;
}

template <class Gt, class Tds >
bool 
Triangulation_2<Gt, Tds>::
includes_edge(Vertex_handle va, Vertex_handle vb,
	      Vertex_handle & vbb,
	      Face_handle& fr, int & i) const
  // returns true if the line segment ab contains an edge e of t 
  // incident to a, false otherwise
  // if true, vbb becomes the vertex of e distinct from a
  // fr is the face incident to e and e=(fr,i)
  // fr is on the right side of a->b
{
  Vertex_handle v;
  Orientation orient;
  int indv;
  Edge_circulator ec = va->incident_edges(), done(ec);
  if (ec != 0) {
    do { 
      //find the index of the other vertex of *ec
      indv = 3 - ((*ec).first)->index(va) - (*ec).second ; 
      v = ((*ec).first)->vertex(indv);
      if (!is_infinite(v)) {
	if (v==vb) {
	  vbb = vb;
	  fr=(*ec).first;
	  i= (*ec).second;
	  return true;
	}
	else {
	  orient = orientation(va->point(),
			   vb->point(),
			   v->point()); 
	  if((orient==COLLINEAR) && 
	     (collinear_between (va->point(),
				 v->point(),
				 vb->point()))) {
	    vbb = v;
	    fr=(*ec).first;
	    i= (*ec).second;
	    return true;
	  }
	}
      }
    } while (++ec != done);
  }
  return false;
}

template <class Gt, class Tds >
inline bool 
Triangulation_2<Gt, Tds>::
is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
{
  return _tds.is_face( &(*(v1)), &(*(v2)), &(*(v3)) );
}

template <class Gt, class Tds >
inline bool 
Triangulation_2<Gt, Tds>::
is_face(Vertex_handle v1, 
	Vertex_handle v2, 
	Vertex_handle v3,
	Face_handle &fr) const
{
  typename Tds::Face* f ;
  bool b = _tds.is_face( &(*(v1)), &(*(v2)), &(*(v3)), f);
  fr = Face_handle(static_cast<Face*>(f));
  return b;
}


template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Triangle
Triangulation_2<Gt, Tds>::
triangle(const Face_handle& f) const
{
  CGAL_triangulation_precondition( ! is_infinite(f) );
  typename Gt::Construct_triangle_2 
     construct_triangle = geom_traits().construct_triangle_2_object();
  return construct_triangle(f->vertex(0)->point(),
			    f->vertex(1)->point(),
			    f->vertex(2)->point());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Face_handle& f, int i) const
{
  CGAL_triangulation_precondition( ! is_infinite(f,i));
  typename Gt::Construct_segment_2 
     construct_segment = geom_traits().construct_segment_2_object();
  return construct_segment(f->vertex(ccw(i))->point(),
			   f->vertex(cw(i))->point());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge& e) const
{
  CGAL_triangulation_precondition(! is_infinite(e)); 
  typename Gt::Construct_segment_2 
     construct_segment = geom_traits().construct_segment_2_object();
  return construct_segment(e.first->vertex(ccw(e.second))->point(),
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
segment(const Finite_edges_iterator& ei) const
{
  return segment(*ei);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const All_edges_iterator& ei) const
{
  return segment(*ei);
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
flip(Face_handle f, int i)
{
  CGAL_triangulation_precondition ( ! f.is_null() );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2); 
    
  CGAL_triangulation_precondition( !is_infinite(f) && 
				   !is_infinite(f->neighbor(i)) );
  CGAL_triangulation_precondition( 
                  orientation(f->vertex(i)->point(),
			      f->vertex(cw(i))->point(),
			      f->mirror_vertex(i)->point()) == RIGHTTURN &&
                  orientation(f->vertex(i)->point(),
			      f->vertex(ccw(i))->point(),
			      f->mirror_vertex(i)->point()) == LEFTTURN); 
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
   return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_in_edge(const Point& p, Face_handle f,int i)
{
  CGAL_triangulation_precondition( orientation(f->vertex(cw(i))->point(), 
					    p,
					    f->vertex(ccw(i))->point()) 
                                   == COLLINEAR &&
      collinear_between(f->vertex(cw(i))->point(), p,
			f->vertex(ccw(i))->point()) );
  Vertex_handle v = static_cast<Vertex*>
                     (_tds.insert_in_edge(&(*f), i));
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_in_face(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition(oriented_side(f,p) ==   ON_POSITIVE_SIDE);
  Vertex_handle v= static_cast<Vertex*>
                     (_tds.insert_in_face( &(*f)));
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
  else   v=insert_outside_convex_hull_2(p, f);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull_1(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition( is_infinite(f) && dimension()==1);
  CGAL_triangulation_precondition(  
    orientation(
	     f->mirror_vertex(f->index(infinite_vertex()))->point(),
	     f->vertex(1- f->index(infinite_vertex()))->point(),
	     p) == COLLINEAR &&
    collinear_between( 
	     f->mirror_vertex(f->index(infinite_vertex()))->point(),
	     f->vertex(1- f->index(infinite_vertex()))->point(),
	     p) );
   Vertex_handle v=static_cast<Vertex*>(_tds.insert_in_edge( &(*f), 2));
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
  CGAL_triangulation_precondition( orientation(p,q,r) == LEFTTURN);

  std::list<Face_handle> ccwlist;
  std::list<Face_handle> cwlist;
    
  Face_circulator fc = infinite_vertex()->incident_faces(f);
  bool done = false;
  while(! done) {
    fc--;
    li = fc->index(infinite_vertex());
    q = fc->vertex(ccw(li))->point();
    r = fc->vertex(cw(li))->point();
    if(orientation(p,q,r) == LEFTTURN ) { ccwlist.push_back(&(*fc)); }
    else {done=true;}
  }

  fc= infinite_vertex()->incident_faces(f);
  done = false;
  while(! done){
    fc++;
    li = fc->index(infinite_vertex());
    q = fc->vertex(ccw(li))->point();
    r = fc->vertex(cw(li))->point();
    if(orientation(p,q,r) == LEFTTURN ) { cwlist.push_back(&(*fc));}
    else {done=true;}
  }

  Vertex_handle v =static_cast<Vertex*>(_tds.insert_in_face( &(*f)));
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
  Face_handle f = (*finite_edges_begin()).first;
  Orientation orient = orientation( f->vertex(0)->point(),
				    f->vertex(1)->point(),
				    p);
  CGAL_triangulation_precondition(orient != COLLINEAR);
  bool conform = ( orient == COUNTERCLOCKWISE);

  Vertex_handle v = static_cast<Vertex*>
    (_tds.insert_dim_up( &(*infinite_vertex()), conform));
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate (p, lt, li, start);
  return insert(p, lt, loc, li);
}

 
template <class Gt, class Tds >
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point& p, Locate_type lt, Face_handle loc, int li)
  // insert a point p, whose localisation is known (lt, f, i)
{
  if(number_of_vertices() == 0) {
    return(insert_first(p));
  }

  if(number_of_vertices() == 1) {
    if (lt == VERTEX) return finite_vertex();
    else return(insert_second(p));
  }

  switch(lt) {
  case FACE:
    return insert_in_face(p,loc);
  case EDGE:
    return insert_in_edge(p,loc,li);
  case OUTSIDE_CONVEX_HULL:
    return  insert_outside_convex_hull(p,loc);
  case OUTSIDE_AFFINE_HULL:
   return  insert_outside_affine_hull(p);
  case VERTEX:
    return loc->vertex(li);
  }
  CGAL_triangulation_assertion(false);  // locate step failed
  return Vertex_handle();
}


template <class Gt, class Tds >
inline
Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Gt, class Tds >
inline void 
Triangulation_2<Gt,Tds>::
remove_degree_3(Vertex_handle  v, Face_handle f)
{
  if (f == Face_handle()) f=v->face();
  _tds.remove_degree_3(&(*v), &(*f));
  return;
}

template <class Gt, class Tds >
inline void
Triangulation_2<Gt,Tds>::
remove_first(Vertex_handle  v)
{
  _tds.remove_second(&(*v));
  return;
}

template <class Gt, class Tds >
inline void 
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
inline void
Triangulation_2<Gt, Tds>::
remove_1D(Vertex_handle v)
{
  _tds.remove_1D(&(*v));
}

template <class Gt, class Tds >
bool
Triangulation_2<Gt,Tds>::
test_dim_down(Vertex_handle v)
{
  //test the dimensionality of the resulting triangulation
  //upon removing of vertex v
  //it goes down to 1 iff
  // 1) any finite face is incident to v
  // 2) all vertices are colinear
  CGAL_triangulation_precondition(dimension() == 2);
  bool  dim1 = true; 
  Finite_faces_iterator fit = finite_faces_begin();
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
	orientation(p, q, fic->vertex(ccw(iv))->point()) == COLLINEAR; 
    }
  }
  return dim1;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::
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

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
make_hole ( Vertex_handle v, std::list<Edge> & hole)
{
  std::list<Face_handle> to_delete;

  Face_handle  f, fn;
  int i, in ;
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
    delete_face(to_delete.front());
    to_delete.pop_front();
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::						     
fill_hole ( Vertex_handle v, std::list< Edge > & hole )
{
  typedef std::list<Edge> Hole;

  Face_handle  ff, fn;
  int ii , in; 
  Vertex_handle v0, v1, v2;
  Bounded_side side;
  

  //stack algorithm to create faces
  // create face v0,v1,v2
  //if v0,v1,v2 are finite vertices
  // and form a leftturn
  // and triangle v0v1v2 does not contain v->point()
  if( hole.size() != 3) {
    typename Hole::iterator hit = hole.begin();
    typename Hole::iterator next= hit; 
    while( hit != hole.end() && hole.size() != 3) {
      ff = (*hit).first;  
      ii = (*hit).second;
      v0 = ff->vertex(cw(ii));
      v1 = ff->vertex(ccw(ii));
      if( !is_infinite(v0) && !is_infinite(v1)) {
	next=hit; next++;
	if(next == hole.end()) next=hole.begin();
	fn = (*next).first; 
	in = (*next).second;
	v2 = fn->vertex(ccw(in));	
	if ( !is_infinite(v2) &&
	     orientation(v0->point(), v1->point(), v2->point()) == LEFTTURN ) {
	  side =  bounded_side(v0->point(), 
			       v1->point(), 
			       v2->point(),
			       v->point());

	  if( side == ON_UNBOUNDED_SIDE || 
	      (side == ON_BOUNDARY && orientation(v0->point(),
					     v->point(),
					     v2->point()) == COLLINEAR &&
	       collinear_between(v0->point(),v->point(),v2->point()) )) 
	    {
	      //create face
	      Face_handle  newf = create_face(ff,ii,fn,in); 
	      typename Hole::iterator tempo=hit;
	      hit = hole.insert(hit,Edge(newf,1)); //push newf
	      hole.erase(tempo); //erase ff
	      hole.erase(next); //erase fn
	      if (hit != hole.begin() ) --hit;
	      continue;
	    }
	}
      }
      ++hit; 
    } 
  }

  // either the hole has only three edges
  // or all its finite vertices are reflex or flat
  // except may be one vertex whose corresponding ear 
  // includes the vertex being removed

  // deal with the last leftturn if any
  if(hole.size() != 3) {
    typename Hole::iterator hit=hole.begin();
    while(hit != hole.end()) {
      ff = (*hit).first;  ii = (*hit).second;
      hit++;
      if(hit != hole.end()) { fn = (*hit).first; in = (*hit).second;}
      else { fn = ((hole.front()).first); in = (hole.front()).second;}
      if ( !is_infinite(ff->vertex(cw(ii))) &&
	   !is_infinite(fn->vertex(cw(in))) &&
	   !is_infinite(fn->vertex(ccw(in))) &&
	   orientation(ff->vertex(cw(ii))->point(),
		       fn->vertex(cw(in))->point(),
		       fn->vertex(ccw(in))->point()) == LEFTTURN) {
	  create_face(ff,ii,fn,in);
	  break;
	}
    }
  }

  // deal with a reflex chain of convex hull edges
  if(hole.size() != 3) {
    // look for infinite vertex
    ff = (hole.front()).first;
    ii = (hole.front()).second;
    while ( ! is_infinite(ff->vertex(cw(ii)))){
      hole.push_back(hole.front());
      hole.pop_front();
      ff = (hole.front()).first;
      ii = (hole.front()).second;
    }
    //create faces
    while(hole.size() != 3){
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      hole.pop_front();
      fn = (hole.front()).first;
      in = (hole.front()).second;
      hole.pop_front();
      Face_handle  newf = create_face(ff,ii,fn,in);
      hole.push_front(Edge(newf,1));
    }
  }
    
  // now hole has three edges
  typename Hole::iterator hit;
  hit = hole.begin();
  //  I don't know why the following yelds a segmentation fault
  //    create_face( (*hit).first, (*hit).second,
  // 	     (* ++hit).first, (*hit).second,
  // 	     (* ++hit).first, (*hit).second);
  ff = (*hit).first;      ii = (*hit).second;
  fn = (* ++hit).first;   in = (*hit).second;
  Face_handle f3 = (* ++hit).first;
  int i3 = (*hit).second;
  create_face(ff,ii,fn,in,f3,i3);
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
fill_hole_delaunay(std::list<Edge> & first_hole)
{
  typename Gt::Side_of_oriented_circle_2
    in_circle = geom_traits().side_of_oriented_circle_2_object();

  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;
  
  Face_handle  f, ff, fn;
  int i, ii, in;
  Hole_list hole_list;
  Hole hole;
        
  hole_list.push_front(first_hole);
  
  while( ! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      typename Hole::iterator hit = hole.begin();
  
      // if the hole has only three edges, create the triangle
      if (hole.size() == 3) {
	hit = hole.begin();
	f = (*hit).first;        i = (*hit).second;
	ff = (* ++hit).first;    ii = (*hit).second;
	fn = (* ++hit).first;    in = (*hit).second;
	create_face(f,i,ff,ii,fn,in);
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
  
      typename Hole::iterator hdone = hole.end();
      hit =  hole.begin();
      typename Hole::iterator cut_after(hit);
  
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
	  if (orientation(p0,p1,p) == COUNTERCLOCKWISE) {
	    if(is_infinite(v2)) { v2=vv; p2=p; cut_after=hit;}
	    else{
	      if( in_circle(p0,p1,p2,p) ==  ON_POSITIVE_SIDE){
		v2=vv; p2=p; cut_after=hit;}
	    }
	  }
	}
	++hit;
      }
 
      // create new triangle and update adjacency relations
       Face_handle newf;
    
      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes
  
      fn = (hole.front()).first;
      in = (hole.front()).second;
      if (fn->has_vertex(v2, i) && i == fn->ccw(in)) {
	newf = create_face(ff,ii,fn,in);
	hole.pop_front();
	hole.push_front(Edge( &(*newf),1));
	hole_list.push_front(hole);
      }
      else{
	fn = (hole.back()).first;
	in = (hole.back()).second;
	if (fn->has_vertex(v2, i) && i== fn->cw(in)) {
	  newf = create_face(fn,in,ff,ii);
	  hole.pop_back();
	  hole.push_back(Edge(&(*newf),1));
	  hole_list.push_front(hole);
	}
	else{
	  // split the hole in two holes
	  newf = create_face(ff,ii,v2);
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
  
template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2,
	 Face_handle f3, int i3)
{
  return static_cast<Face*>(_tds.create_face(&(*f1),i1, &(*f2),i2, &(*f3),i3));
}

template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2)
{
  return static_cast<Face*>(_tds.create_face(&(*f1),i1, &(*f2),i2));
}

template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f, int i, Vertex_handle v)
{
  return static_cast<Face*>(_tds.create_face(&(*f),i, &(*v)));
}

template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  return static_cast<Face*>(_tds.create_face(&(*v1), &(*v2), &(*v3)));
}

template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
	    Face_handle f1, Face_handle f2,  Face_handle f3)
{
  return static_cast<Face*>(_tds.create_face(&(*v1), &(*v2), &(*v3),
					     &(*f1), &(*f2), &(*f3)));
}

template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle fh)
{
  return static_cast<Face*>(_tds.create_face(&(*fh)));
}



template <class Gt, class Tds >    
inline
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face()
{
  return static_cast<Face*>(_tds.create_face());
}

template <class Gt, class Tds >    
inline
void
Triangulation_2<Gt, Tds>::
delete_face(Face_handle f)
{
  _tds.delete_face(&(*f));
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
  Orientation pqt = orientation(f->vertex(0)->point(), 
				f->vertex(1)->point(),
				t);
  if(pqt == RIGHTTURN || pqt == LEFTTURN) {
    lt = OUTSIDE_AFFINE_HULL;
    li = 4 ;// should not be used
    return Face_handle();
  }

  int i= f->index(ff);
  if (collinear_between(t,f->vertex(1-i)->point(),f->vertex(i)->point())) {
    lt = OUTSIDE_CONVEX_HULL;
    li = iv;
    return ff;
  }

  if( xy_equal(t,f->vertex(1-i)->point()) ){
    lt = VERTEX;
    li=1-i;
    return f;
  }

  ff = ff->neighbor(1-iv); //the other infinite face
  iv = ff->index(infinite_vertex());
  f = ff->neighbor(iv);
  i = f->index(ff);
  if (collinear_between(t,f->vertex(1-i)->point(),f->vertex(i)->point())) {
    lt = OUTSIDE_CONVEX_HULL;
    li = iv;
    return ff;
  }
  if( xy_equal(t,f->vertex(1-i)->point()) ){
      lt = VERTEX;
      li=1-i;
    return f;
  } 
	
  Finite_edges_iterator eit= finite_edges_begin();
  Vertex_handle u,v;
  for( ; eit != finite_edges_end() ; eit++) {
    u = (*eit).first->vertex(0);
    v = (*eit).first->vertex(1);
    if(xy_equal(t,v->point())){
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
  return Face_handle();
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
  Point p(start->vertex(0)->point());
  if(xy_equal(t,p)) {
    lt = VERTEX;
    li = 0;
    return start;
  }

  Line_face_circulator lfc(start->vertex(0), this, t);

  if(lfc==0 || lfc.collinear_outside()){
    // point t lies outside or on the convex hull
    // we walk on the convex hull to find it out
    Face_circulator fc = infinite_vertex()->incident_faces();
    Face_circulator done(fc);
    int ic = fc->index(infinite_vertex());
    if (xy_equal(t,fc->vertex(cw(ic))->point())){
      lt = VERTEX;
      li = cw(ic);
      return fc;
     }
    Orientation ori;
    do{ // walking ccw around convex hull
      ic = fc->index(infinite_vertex());
      if (xy_equal(t,fc->vertex(ccw(ic))->point())){
	lt = VERTEX;
	li = ccw(ic);
	return fc;
      }
      ori = orientation( fc->vertex(cw(ic))->point(),
			 fc->vertex(ccw(ic))->point(), t);
      if (ori == RIGHTTURN) {
	lt = OUTSIDE_CONVEX_HULL;
	li = ic;
	return fc;
      }
      if (ori == COLLINEAR &&
	  collinear_between(fc->vertex(cw(ic))->point(),
			    t, 
			    fc->vertex(ccw(ic))->point()) ) {
	lt = EDGE;
	li = ic;
	return fc;
      }
    } while (--fc != done);
    //should not arrive there;
    CGAL_triangulation_assertion(fc != done);
  }
	  
    while(! lfc.locate(t, lt, li) ){
    ++lfc;
  }
  return lfc;
}    
      
//   bool collinear_convex_hull_edge= false;
//   Face_handle collinear_fh;
//   int collinear_ih;
//   Line_face_circulator lfc(start->vertex(0), this, t);
//   if (lfc == 0) {
//     // point t lies outside or on the convex hull
//     Face_circulator fc = start->vertex(0).begin();
//     while(!is_infinite(fc)) fc++;
//     int i = fc->index(infinite_vertex());
//     Orientation ori = orientation(t,
// 				  fc->vertex(ccw(i))->point(),
// 				  fc->vertex(ccw(i))->point());
//     if (ori= RIGHTTURN) {
//       fc++;
//       i = fc->index(infinite_vertex());
//       ori = orientation(t,
// 			fc->vertex(ccw(i))->point(),
// 			fc->vertex(ccw(i))->point());
//     }
//     //ori is now LEFTTURN or COLLINEAR
//     CGAL_Triangulation_assertion( ori == LEFTTURN || ori == COLLINEAR );
//     if (ori == LEFTTURN) {
//       lt = OUTSIDE_CONVEX_HULL;
//       li = i;
//       return fc;
//     }
//     else { 
//       // segment [start->vertex(0)t] is collinear to a convex hull
//       //edge
//       collinear_convex_hull_edge = true;
//       collinear_fh = fc;
//       collinear_ih = i;
//     }
//   }

//   if (lfc.collinear_outside() {
//     collinear_convex_hull_edge = true;
//     collinear_fh = lfc->neighbor(ccw(i));
//     collinear_ih = collinear_fh->index(infinite_vertex());
//   }
	
//   if(collinear_convex_hull_edge) {
//     // point t lies outside or on the convex hull
//     // we walk on the hull to decide
//     Face_circulator fc =
//       incident_vertex()->incident_faces(collinear_fh);
//     Point p,q;
//     int ic;
//     do {
//       ic = fc->index(infinite_vertex());
//       r = fc->vertex(ccw(ic))->point();
//       q = fc->vertex(cw(ic))->point();
//       fc = fc++;
//       ic = fc->index(infinite_vertex());
//       p = fc->vertex(cw(ic))->point();
//     } while(orientation(p, q, r) == COLLINEAR);
//     if (orientation(p, q, t) == LEFTTURN){
//       	lt = OUTSIDE_CONVEX_HULL;
// 	li = ic;
// 	return fc ;
//     }
    
//     fc--;
//     ic = fc->index(infinite_vertex());
//     p = fc->vertex(cw(ic))->point();
//     q = fc->vertex(ccw(ic))->point();
//     if(xy_equal(t,p)){
//       lt = VERTEX;
//       li = cw(ic);
//       return lfc;
//     }
//     Orientation pqt;
//     while(1) {
//       if(xy_equal(t,q)){
// 	lt = VERTEX;
// 	li = ccw(i);
// 	return f;
//       }
//       pqt = orientation(p,q,t);
//       if (pqt == COLLINEAR && collinear_between(p, t, q)){
// 	lt = EDGE;
// 	li = ic;
// 	return fc;
//       }
//       if (pqt == RIGHTTURN){
// 	lt = OUTSIDE_CONVEX_HULL;
// 	li = ic;
// 	return fc ;
//       }
	   	       
//       // go to the next face
//       fc--;
//       ic = fc->index(infinite_vertex());
//       p = q;
//       q = fc->vertex(ccw(i))->point();
//     }
//   }

//   while(! lfc.locate(t, lt, li) ){
//     ++lfc;
//   }
//   return lfc;
//}

    
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt,Tds>::
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle start) const
{
  if( dimension() <= 0) {
    if(number_of_vertices() == 0) {
      lt = OUTSIDE_AFFINE_HULL;
      li = 4; // li should not be used in this case
    } else { // number_of_vertices() == 1
      if (xy_equal(p,finite_vertex()->point())){
	lt = VERTEX ;
      }
      else{
	lt = OUTSIDE_AFFINE_HULL;
      }
      li = 4; // li should not be used in this case
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
       Face_handle start) const
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
Triangulation_2<Gt, Tds>::Point_iterator
Triangulation_2<Gt, Tds>::
points_begin() const
{
  return Point_iterator(finite_vertices_begin());
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Point_iterator
Triangulation_2<Gt, Tds>::
points_end() const
{
  return Point_iterator(finite_vertices_end());
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
incident_faces(Vertex_handle v, Face_handle f) const
{
  return Face_circulator(v,f);
}  

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Vertex_circulator
Triangulation_2<Gt, Tds>::  
incident_vertices(Vertex_handle v,
		  Face_handle f) const
{
  return Vertex_circulator(v,f);
}

template <class Gt, class Tds >
inline
Triangulation_2<Gt, Tds>::Edge_circulator
Triangulation_2<Gt, Tds>::    
incident_edges(Vertex_handle v,
	       Face_handle f) const
{
  return Edge_circulator(v,f);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Line_face_circulator  
Triangulation_2<Gt, Tds>::    
line_walk(const Point& p, const Point& q,  Face_handle f) const
{
  CGAL_triangulation_precondition( (dimension() == 2) && 
				! xy_equal(p,q));
  Line_face_circulator lfc = (f.is_null())
    ? Line_face_circulator(p, q, this)
    : Line_face_circulator(p, q, f, this);
    
  // the following lines may be useless :
  //  Line_face_circulator(p,q...) returns either a null circulator 
  //  or a pointer to a finite face (to be checked)
  if( (!lfc.is_empty()) && is_infinite( lfc )){
    do {      ++lfc ;} 
    while (is_infinite(lfc));
  }
  return lfc;
}
   
template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
oriented_side(const Point &p0, const Point &p1,
	      const Point &p2, const Point &p) const
{
  // return position of point p with respect to the oriented triangle p0p1p2
  // depends on the orientation of the vertices
  Bounded_side bs=bounded_side(p0,p1,p2,p);
  if (bs == ON_BOUNDARY) return ON_ORIENTED_BOUNDARY;
  Orientation      ot = orientation(p0, p1, p2);
  if (bs == ON_BOUNDED_SIDE)
    return (ot == LEFTTURN) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  // bs == ON_UNBOUNDED_SIDE
  return (ot == LEFTTURN) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
}



template <class Gt, class Tds >
Bounded_side
Triangulation_2<Gt, Tds>::
bounded_side(const Point &p0, const Point &p1,
	     const Point &p2, const Point &p) const
{
  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition( orientation(p0, p1, p2) != COLLINEAR);
  Orientation o1 = orientation(p0, p1, p),
              o2 = orientation(p1, p2, p),
              o3 = orientation(p2, p0, p);
    
  if (o1 == COLLINEAR){
    if (o2 == COLLINEAR ||  o3 == COLLINEAR) return ON_BOUNDARY;
    if (collinear_between(p0, p, p1))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o2 == COLLINEAR){
    if (o3 == COLLINEAR)                     return ON_BOUNDARY;
    if (collinear_between(p1, p, p2))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o3 == COLLINEAR){
    if (collinear_between(p2, p, p0))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  // from here none ot, o1, o2 and o3 are known to be non null
    if (o1 == o2 && o2 == o3)  return ON_BOUNDED_SIDE;
    return ON_UNBOUNDED_SIDE;
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

template < class Gt, class Tds >
Oriented_side
Triangulation_2<Gt,Tds>::
side_of_oriented_circle(Face_handle f, const Point & p) const
{
  if ( ! is_infinite(f) ) {
    typename Gt::Side_of_oriented_circle_2 
      in_circle = geom_traits().side_of_oriented_circle_2_object();
    return in_circle(f->vertex(0)->point(),
		     f->vertex(1)->point(),
		     f->vertex(2)->point(),p);
  }

  int i = f->index(infinite_vertex());
  Orientation o = orientation(f->vertex(ccw(i))->point(),
			      f->vertex(cw(i))->point(),
			      p);
  return (o == NEGATIVE) ? ON_NEGATIVE_SIDE :
                (o == POSITIVE) ? ON_POSITIVE_SIDE :
                      ON_ORIENTED_BOUNDARY;
}



template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{
  // return true if point q is strictly between p and r
  // p,q and r are supposed to be collinear points
  Comparison_result c_pr = compare_x(p, r);
  Comparison_result c_pq;
  Comparison_result c_qr;
  if(c_pr == EQUAL) {
    //c_pr = compare_y(p, r);
    c_pq = compare_y(p, q);
    c_qr = compare_y(q, r);
  } else {
    c_pq = compare_x(p, q);
    c_qr = compare_x(q, r);
  }
  return ( (c_pq == SMALLER) && (c_qr == SMALLER) ) ||
         ( (c_pq == LARGER)  && (c_qr == LARGER) );
    
}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_2<Gt, Tds>::
compare_x(const Point& p, const Point& q) const
{
  return geom_traits().compare_x_2_object()(p,q);
}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_2<Gt, Tds>::
compare_y(const Point& p, const Point& q) const
{
  return geom_traits().compare_y_2_object()(p,q);
}

template <class Gt, class Tds >
inline
bool
Triangulation_2<Gt, Tds>::
xy_equal(const Point& p, const Point& q) const
{
  return compare_x(p,q)== EQUAL && compare_y(p,q)== EQUAL ;
}

template <class Gt, class Tds >
inline
Orientation
Triangulation_2<Gt, Tds>::
orientation(const Point& p, const Point& q,const Point& r ) const
{
  return geom_traits().orientation_2_object()(p,q,r);
}

template<class Gt, class Tds>
inline
Triangulation_2<Gt,Tds>::Point
Triangulation_2<Gt,Tds>::
circumcenter (const Point& p0, const Point& p1, const Point&  p2) const
{
  return 
    geom_traits().construct_circumcenter_2_object()(p0,p1,p2);
}


template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Point
Triangulation_2<Gt, Tds>::
circumcenter(Face_handle  f) const
{
  CGAL_triangulation_precondition (dimension()==2);
  // typename Gt::Construct_circumcenter_2
//     circumcenter = geom_traits().construct_circumcenter_2_object();
  return circumcenter((f->vertex(0))->point(), 
		      (f->vertex(1))->point(), 
		      (f->vertex(2))->point());
}

 
template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_all()
{
  std::cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<std::endl;
  All_faces_iterator fi = all_faces_begin();
  std::cerr<<"***"<<std::endl;
  while(fi != all_faces_end()) {
    show_face(fi);
    ++fi;
  }
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_face(Face_handle fh)
{
  std::cerr << "face : "<<(void*)&(*fh)<<" => "<<std::endl;
  int i = fh->dimension(); 
  switch(i){
  case 0:
    std::cerr <<"point :"<<(fh->vertex(0)->point())
	      <<" / voisin "<<&(*(fh->neighbor(0)))
	      <<"["<<(fh->neighbor(0))->vertex(0)->point()<<"]"
	      <<std::endl;
    break;
  case 1:
     std::cerr <<"point :" <<(fh->vertex(0)->point())
	       <<" / voisin "<<&(*(fh->neighbor(0)))
	       <<"["<<(fh->neighbor(0))->vertex(0)->point()
	       <<"/"<<(fh->neighbor(0))->vertex(1)->point()<<"]"
	       <<std::endl;
     std::cerr <<"point :"<<(fh->vertex(1)->point())
	       <<" / voisin "<<&(*(fh->neighbor(1)))
	       <<"["<<(fh->neighbor(1))->vertex(0)->point()
	       <<"/"<<(fh->neighbor(1))->vertex(1)->point()<<"]"
	       <<std::endl;
     break;
  case 2:
  std::cerr <<"point :"<<(fh->vertex(0)->point())
	    <<" / voisin "<<&(*(fh->neighbor(0)))
	    <<"["<<(fh->neighbor(0))->vertex(0)->point()
	    <<"/"<<(fh->neighbor(0))->vertex(1)->point()
	    <<"/"<<(fh->neighbor(0))->vertex(2)->point()<<"]"
	    <<std::endl;
  std::cerr <<"point :"<<(fh->vertex(1)->point())
	    <<" / voisin "<<&(*(fh->neighbor(1)))
	    <<"["<<(fh->neighbor(1))->vertex(0)->point()
	    <<"/"<<(fh->neighbor(1))->vertex(1)->point()
	    <<"/"<<(fh->neighbor(1))->vertex(2)->point()<<"]"
	    <<std::endl;
  std::cerr <<"point :"<<(fh->vertex(2)->point())
	    <<" / voisin "<<&(*(fh->neighbor(2)))
	    <<"["<<(fh->neighbor(2))->vertex(0)->point()
	    <<"/"<<(fh->neighbor(2))->vertex(1)->point()
	    <<"/"<<(fh->neighbor(2))->vertex(2)->point()<<"]"
	    <<std::endl;
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
file_output(std::ostream& os) const
{
  _tds.file_output(os, &(*infinite_vertex()), true);
}

template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::Vertex_handle
Triangulation_2<Gt, Tds>::
file_input(std::istream& is)
{
  clear();
  Vertex_handle v= static_cast<Vertex*>(_tds.file_input(is, true));
  set_infinite_vertex(v);
  return v;
}

template <class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Triangulation_2<Gt, Tds> &tr)
{
  tr.file_output(os);
  return os ;
}



template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_2<Gt, Tds> &tr)
{
  tr.file_input(is);
  CGAL_triangulation_assertion(tr.is_valid());
  return is;
}
 
CGAL_END_NAMESPACE
    

#endif //CGAL_TRIANGULATION_2_H

