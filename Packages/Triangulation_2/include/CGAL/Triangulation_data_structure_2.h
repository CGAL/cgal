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
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Triangulation_data_structure_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_2_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_2_H

#include <CGAL/basic.h>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <algorithm>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Trivial_iterator.h>
#include <CGAL/DS_Container.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_ds_face_2.h>
#include <CGAL/Triangulation_ds_vertex_2.h>
#include <CGAL/Triangulation_ds_handles_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>

CGAL_BEGIN_NAMESPACE 
template < class Vb, class Fb>
class Triangulation_data_structure_2;

template < class Vb, class Fb>
std::istream& operator>>
(std::istream& is, Triangulation_data_structure_2<Vb,Fb>&
  tds);

template < class Vb, class Fb>
std::ostream& operator<< 
( std::ostream& os, 
  const Triangulation_data_structure_2<Vb,Fb>& tds);

//for backward compatibility
template < class Gt , class Vb, class Fb>
class Triangulation_default_data_structure_2 
  : public Triangulation_data_structure_2<Vb,Fb>
{
public:
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;
  typedef Triangulation_default_data_structure_2<Gt,Vb,Fb> Tdds;
  typedef Gt                           Geom_traits; 

  Triangulation_default_data_structure_2(const Geom_traits&
					 gt=Geom_traits())
    : Tds() {}
 
  Triangulation_default_data_structure_2(const Tdds &tdds)
    : Tds(tdds) {}
};

//for backward compatibility
template <class Vb, class Fb>
class Triangulation_data_structure_using_list_2
  :public Triangulation_data_structure_2<Vb, Fb>
{
public:
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;
  typedef Triangulation_data_structure_using_list_2<Vb,Fb>  Tdsul;

  Triangulation_data_structure_using_list_2(): Tds() {} 
  Triangulation_data_structure_using_list_2(const Tdsul &tdsul)
    : Tds(tdsul) {}
};

 

template <class Vb, class Fb>
class Triangulation_data_structure_2 
  :public Triangulation_cw_ccw_2
{
public:
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;

  friend class Triangulation_ds_face_iterator_2<Tds>;
  friend class Triangulation_ds_edge_iterator_2<Tds>;
  friend class Triangulation_ds_vertex_iterator_2<Tds>;
                   
  typedef Vb                                         Vertex_base;
  typedef Fb                                         Face_base;
  typedef Triangulation_ds_vertex_2<Tds>             Vertex;
  typedef Triangulation_ds_face_2<Tds>               Face;
  typedef CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER_POINTER(Face)   Face_handle;
  typedef CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER_POINTER(Vertex)
                                                     Vertex_handle;
  typedef std::pair<Face_handle, int>                Edge;
  
  typedef DS_Container<Face>                         Face_container;
  typedef DS_Container<Vertex>                       Vertex_container;
  typedef typename Face_container::iterator          Iterator_base;
  typedef Triangulation_ds_face_iterator_2<Tds>      Face_iterator;
  typedef Triangulation_ds_edge_iterator_2<Tds>      Edge_iterator;
  typedef Triangulation_ds_vertex_iterator_2<Tds>    Vertex_iterator;

  typedef Triangulation_ds_face_circulator_2<Tds>    Face_circulator;
  typedef Triangulation_ds_vertex_circulator_2<Tds>  Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_2<Tds>    Edge_circulator;
  typedef std::list<Edge> List_edges;

protected:
  int _number_of_vertices; 
  int _dimension;
  Face_container   _face_container;
  Vertex_container _vertex_container;

  //CREATORS - DESTRUCTORS
public:
  Triangulation_data_structure_2(); 
  Triangulation_data_structure_2(const Tds &tds);
  ~Triangulation_data_structure_2();
  Tds& operator= (const Tds &tds);
  void swap(Tds &tds);

  //ACCESS FUNCTIONS
private:
  Face_container& face_container()             { return _face_container;}
  const Face_container& face_container() const { return  _face_container;}
  Vertex_container& vertex_container()         {return _vertex_container;}
  const Vertex_container& vertex_container() const
                                               {return  _vertex_container;}

public:
  int  dimension() const { return _dimension;  }
  int number_of_vertices() const {return _number_of_vertices;}
  int number_of_faces() const ;
  int number_of_edges() const;
  int number_of_full_dim_faces() const; //number of faces stored by tds
  
  // TEST FEATURES
  bool is_vertex(const Vertex_handle v) const;
  bool is_edge(const Vertex_handle va, const Vertex_handle vb) const;
  bool is_edge(const Vertex_handle va, const Vertex_handle vb, 
	       Face_handle& fr,  int& i) const;
  bool is_face(const Vertex_handle v1, 
	       const Vertex_handle v2, 
	       const Vertex_handle v3) const;
  bool is_face(const Vertex_handle v1, 
	       const Vertex_handle v2, 
	       const Vertex_handle v3,
	       Face_handle& fr) const;

  // ITERATORS AND CIRCULATORS
public:
// The Iterator_base  gives the possibility to iterate over all
// in the container  independently of the dimension.
  // public for the need of file_ouput() of Constrained triangulation
  // should be made public later

  Iterator_base iterator_base_begin() const    {
    return face_container().begin();
  }
  Iterator_base iterator_base_end() const    {
    return face_container().end();
  }

public:
  Face_iterator faces_begin() const {
    if (dimension() < 2) return faces_end();
    return Face_iterator(this);
  }
    
  Face_iterator faces_end() const {
    return Face_iterator(this, 1);
  }

  Vertex_iterator vertices_begin() const  {
    return Vertex_iterator(this);
  }

  Vertex_iterator vertices_end() const {
    return Vertex_iterator(this,1);
  }
  
  Edge_iterator edges_begin() const {
    return Edge_iterator(this);
  }

  Edge_iterator edges_end() const {
    return Edge_iterator(this,1);
  }
  
  Face_circulator incident_faces(Vertex_handle v, 
				 Face_handle f =  Face_handle(NULL)) const{
    return v->incident_faces(f);
  }
  Vertex_circulator incident_vertices(Vertex_handle v, 
				      Face_handle f = Face_handle(NULL)) const
  {    
    return v->incident_vertices(f);  
  }

  Edge_circulator incident_edges(Vertex_handle v, 
				 Face_handle f = Face_handle(NULL)) const{
    return v->incident_edges(f);
  }

  // MODIFY
  void flip(Face_handle f, int i);
 
  Vertex_handle insert_first();
  Vertex_handle insert_second();
  Vertex_handle insert_in_face(Face_handle f);
  Vertex_handle insert_in_edge(Face_handle f, int i);
  Vertex_handle insert_dim_up(Vertex_handle w = Vertex_handle(NULL), 
			      bool orient=true);

  void remove_degree_3(Vertex_handle v, Face_handle f = Face_handle(NULL));
  void remove_1D(Vertex_handle v); 
   
  void remove_second(Vertex_handle v);
  void remove_first(Vertex_handle v);
  void remove_dim_down(Vertex_handle v);

  Vertex_handle star_hole(List_edges& hole);
  void    star_hole(Vertex_handle v, List_edges& hole);
  void    make_hole(Vertex_handle v, List_edges& hole);

//   template< class EdgeIt>
//   Vertex_handle star_hole(EdgeIt edge_begin,EdgeIt edge_end);
 
//   template< class EdgeIt>
//   void  star_hole(Vertex_handle v, EdgeIt edge_begin,  EdgeIt edge_end);

//   template< class EdgeIt, class FaceIt>
//   Vertex_handle star_hole(EdgeIt edge_begin, 
// 		    EdgeIt edge_end,
// 		    FaceIt face_begin,
// 		    FaceIt face_end);
 
//   template< class EdgeIt, class FaceIt>
//   void  star_hole(Vertex_handle v,
// 		  EdgeIt edge_begin, 
// 		  EdgeIt edge_end,
// 		  FaceIt face_begin,
// 		  FaceIt face_end);
  
  Vertex_handle create_vertex();
  Face_handle create_face(Face_handle f1, int i1, 
			  Face_handle f2, int i2, 
			  Face_handle f3, int i3);
  Face_handle create_face(Face_handle f1, int i1, 
			  Face_handle f2, int i2);
  Face_handle create_face(Face_handle f1, int i1, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, 
			  Vertex_handle v2, 
			  Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, 
			  Vertex_handle v2, 
			  Vertex_handle v3,
			  Face_handle f1, 
			  Face_handle f2, 
			  Face_handle f3);
  Face_handle create_face(Face_handle f); //calls copy constructor of Face
  Face_handle create_face();
  void set_adjacency(Face_handle f0, int i0, Face_handle f1, int i1) const;
  void delete_face(Face_handle);
  void delete_vertex(Vertex_handle);

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  
  // HELPING
private:
  typedef std::pair<Vertex_handle,Vertex_handle> Vh_pair;
  void  set_adjacency(Face_handle fh, 
		      int ih, 
		      std::map< Vh_pair, Edge>& edge_map);
  void reorient_faces();

public:
  void clear();
  Vertex_handle copy_tds(const Tds &tds, Vertex_handle = Vertex_handle(NULL));
  
  // I/O
  Vertex_handle file_input(std::istream& is, bool skip_first=false);
  void file_output(std::ostream& os,
		   Vertex_handle v = Vertex_handle(NULL),
		   bool skip_first=false) const;
  Vertex_handle off_file_input(std::istream& is, bool verbose=false);
  void  vrml_output(std::ostream& os,
		    Vertex_handle v = Vertex_handle(NULL),
		    bool skip_first=false) const;


  // SETTING (had to make them public for use in remove from Triangulations)
  void set_number_of_vertices(int n) {_number_of_vertices = n;}
  void set_dimension (int n) {_dimension = n ;}

  // template members definition
public:
  template< class EdgeIt>
  Vertex_handle star_hole(EdgeIt edge_begin, EdgeIt edge_end)
  // creates a new vertex 
  // and stars from it
  // the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
  {
     Vertex_handle newv = create_vertex();
     star_hole(newv, edge_begin, edge_end);
     return newv;
  }
 
  template< class EdgeIt>
  void  star_hole(Vertex_handle v, EdgeIt edge_begin,  EdgeIt edge_end)
  // uses vertex v
  // to star the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // the hole is supposed to be ccw oriented
  { 
    std::list<Face_handle> empty_list;
    star_hole(v, 
	      edge_begin, 
	      edge_end, 
	      empty_list.begin(),
	      empty_list.end());
    return;    
  }


  template< class EdgeIt, class FaceIt>
  Vertex_handle star_hole(EdgeIt edge_begin, 
		    EdgeIt edge_end,
		    FaceIt face_begin,
		    FaceIt face_end)
  // creates a new vertex 
  // and stars from it
  // the hole described by the range [edge_begin,edge_end[
    // reusing the faces in the range [face_begin,face_end[
  // the triangulation is assumed to have dim=2
  // the hole is supposed to be ccw oriented
  {
    Vertex_handle newv = create_vertex();
    star_hole(newv, edge_begin, edge_end, face_begin, face_end);
    return newv;
  }
 
  template< class EdgeIt, class FaceIt>
  void  star_hole(Vertex_handle newv,
		  EdgeIt edge_begin, 
		  EdgeIt edge_end,
		  FaceIt face_begin,
		  FaceIt face_end)
    // uses vertex v
    // to star the hole described by the range [edge_begin,edge_end[
    // reusing the faces in the range [face_begin,face_end[
    // the triangulation is assumed to have dim=2
    // hole is supposed to be ccw oriented
  {
    CGAL_triangulation_precondition(dimension() == 2);
    EdgeIt eit = edge_begin;
    FaceIt fit = face_begin;

    Face_handle fn = (*eit).first;
    int in = (*eit).second;
    fn->vertex(cw(in))->set_face(fn);
    Face_handle first_f =  reset_or_create_face(fn, in , newv, fit, face_end);
    Face_handle previous_f=first_f, next_f;
    ++eit; 

    for( ; eit != edge_end ; eit++) {
      fn = (*eit).first;
      in = (*eit).second;
      fn->vertex(cw(in))->set_face(fn);
      next_f = reset_or_create_face(fn, in , newv, fit, face_end);
      next_f->set_neighbor(1, previous_f);
      previous_f->set_neighbor(0, next_f);
      previous_f=next_f;
    }
    
    next_f->set_neighbor(0, first_f);
    first_f->set_neighbor(1, next_f);
    newv->set_face(first_f);
    return;    
  }

private:
  template< class FaceIt>
  Face_handle  reset_or_create_face(Face_handle fn, 
			      int in, 
			      Vertex_handle v,
			      FaceIt& fit,
			      const FaceIt& face_end)
  {
    if (fit == face_end) return create_face(fn, in, v);
    (*fit)->set_vertices(fn->vertex(cw(in)), fn->vertex(ccw(in)), v);
    (*fit)->set_neighbors(0,0,fn);
    fn->set_neighbor(in, *fit);
    return &(**fit++);    
  }

};

template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
Triangulation_data_structure_2() 
  :  _number_of_vertices(0),_dimension(-2)
{ }

template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
Triangulation_data_structure_2(const Tds &tds)
{
  copy_tds(tds);
}

template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
~Triangulation_data_structure_2()
{
  clear();
}

//assignement  
template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb>&
Triangulation_data_structure_2<Vb,Fb> ::
operator= (const Tds &tds)
{
  copy_tds(tds);
  return *this;
}  

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
clear()
{
  face_container().clear();
  vertex_container().clear();

  set_number_of_vertices(0);
  set_dimension(-2);

  return;
}

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
swap(Tds &tds)
{
  CGAL_triangulation_expensive_precondition(tds.is_valid() && is_valid());
  std::swap(_dimension, tds._dimension);
  std::swap(_number_of_vertices, tds._number_of_vertices);
  face_container().swap(tds.face_container());
  vertex_container().swap(tds.vertex_container());
  return;
}

//ACCESS FUNCTIONS
template <class Vb, class Fb>
inline int 
Triangulation_data_structure_2<Vb,Fb> ::
number_of_faces() const 
{
  if (dimension() < 2) return 0;
  return face_container().size();
}

template <class Vb, class Fb>
int
Triangulation_data_structure_2<Vb,Fb>::
number_of_edges() const
{
  switch (dimension()) {
  case 1:  return number_of_vertices();
  case 2:  return 3*number_of_faces()/2;
  default: return 0;
  }
}
      
template <class Vb, class Fb>
int
Triangulation_data_structure_2<Vb,Fb>::
number_of_full_dim_faces() const
{
//   switch (dimension()) {
//   case -2: return 0;
//   case -1: return 1;
//   case 0:  return 2;
//   case 1:  return number_of_edges();
//   case 2:  return number_of_faces();
//   default: return 0;
//   }
  return face_container().size();
}

template <class Vb, class Fb>
inline bool
Triangulation_data_structure_2<Vb,Fb>::
is_vertex(const Vertex_handle v) const
{
  return vertex_container().is_element(&*v);
}

template <class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_edge(const Vertex_handle va, const Vertex_handle vb) const
// returns true (false) if the line segment ab is (is not) an edge of t
//It is assumed that va is a vertex of t
{
  Vertex_circulator vc= va->incident_vertices(), done(vc);
  if ( vc == 0) return false;
  do {
    if( vb == vc->handle() ) {return true;} 
  } while (++vc != done);
  return false;
}
 

template <class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_edge(const Vertex_handle va, const Vertex_handle vb, 
	Face_handle &fr,  int & i) const
// assume va is a vertex of t
// returns true (false) if the line segment ab is (is not) an edge of t
// if true is returned (fr,i) is the edge ab
// with face fr on the right of a->b
{
  Face_handle fc = va->face(); 
  Face_handle start = fc;
  if (fc == 0) return false;
  int inda, indb;
  do {
    inda=fc->index(va);
    indb = (dimension() == 2 ? cw(inda) : 1-inda);
    if(fc->vertex(indb) == vb) {
      fr=fc;
      i = 3 - inda - indb; //works in dim 1 or 2
      return true;
    }
    fc=fc->neighbor(indb); //turns ccw around va
  } while (fc != start);
  return false;
}


template <class Vb, class Fb>
inline bool 
Triangulation_data_structure_2<Vb,Fb>::
is_face(const Vertex_handle v1, 
	const Vertex_handle v2, 
	const Vertex_handle v3) const
{
  Face_handle f;
  return is_face(v1,v2,v3,f);
}

template <class Vb, class Fb>
bool 
Triangulation_data_structure_2<Vb,Fb>::
is_face(const Vertex_handle v1, 
	const Vertex_handle v2, 
	const Vertex_handle v3,
	Face_handle &f) const
{
  if (dimension() != 2) return false;
  int i;
  bool b = is_edge(v1,v2,f,i);
  if (!b) return false;
  else if (v3== f->vertex(i)) return true;
  f = f-> neighbor(i);
  int ind1= f->index(v1);
  int ind2= f->index(v2);
  if (v3 == f->vertex(3-ind1-ind2)) { return true;}
  return false;  
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
flip(Face_handle f, int i)
{
  CGAL_triangulation_precondition( dimension()==2);
  Face_handle n  = f->neighbor(i);
  int ni = f->mirror_index(i); //ni = n->index(f);
    
  Vertex_handle  v_cw = f->vertex(cw(i));
  Vertex_handle  v_ccw = f->vertex(ccw(i));

  // bl == bottom left, tr == top right
  Face_handle tr = f->neighbor(ccw(i));
  int tri =  f->mirror_index(ccw(i));  
  Face_handle bl = n->neighbor(ccw(ni));
  int bli =  n->mirror_index(ccw(ni)); 
      
  f->set_vertex(cw(i), n->vertex(ni));
  n->set_vertex(cw(ni), f->vertex(i));
    
  // update the neighborhood relations
  f->set_neighbor(i, bl);
  bl->set_neighbor(bli, f);
    
  f->set_neighbor(ccw(i), n);
  n->set_neighbor(ccw(ni), f);
    
  n->set_neighbor(ni, tr);
  tr->set_neighbor(tri, n);
    
  if(v_cw->face() == f) {
    v_cw->set_face(n);
  }
    
  if(v_ccw->face() == n) {
    v_ccw->set_face(f);
  }
}
  
template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_first( )
{
  CGAL_triangulation_precondition( number_of_vertices() == 0 &&
				   dimension()==-2 );
  return insert_dim_up();
}

template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle 
Triangulation_data_structure_2<Vb,Fb>::
insert_second()
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 &&
				   dimension()==-1 );
  return insert_dim_up();

}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_in_face(Face_handle f)
  // New vertex will replace f->vertex(0) in face f
{
  CGAL_triangulation_precondition( f != NULL && dimension()== 2);
  Vertex_handle  v = create_vertex();

  Vertex_handle v0 = f->vertex(0);
  Vertex_handle v2 = f->vertex(2);
  Vertex_handle v1 = f->vertex(1);
    
  Face_handle n1 = f->neighbor(1);
  Face_handle n2 = f->neighbor(2);
    
  Face_handle f1 = create_face(v0, v, v2, f, n1, NULL);
  Face_handle f2 = create_face(v0, v1, v, f, NULL, n2);

  f1->set_neighbor(2, f2);
  f2->set_neighbor(1, f1);
  if (n1 != NULL) {
    int i1 = f->mirror_index(1); //int i1 = n1->index(f);
    n1->set_neighbor(i1,f1);
  }
  if (n2 != NULL) {
    int i2 = f->mirror_index(2);//int i2 = n2->index(f);
    n2->set_neighbor(i2,f2);}

  f->set_vertex(0, v);
  f->set_neighbor(1, f1);
  f->set_neighbor(2, f2);

  if( v0->face() == f  ) {  v0->set_face(f2); }
  v->set_face(f);

  return v;
}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_in_edge(Face_handle f, int i)
  //insert in the edge opposite to vertex i of face f
{
  CGAL_triangulation_precondition(f != NULL && dimension() >= 1); 
  if (dimension() == 1) {CGAL_triangulation_precondition(i == 2);}
  if (dimension() == 2) {CGAL_triangulation_precondition(i == 0 || 
							 i == 1 || 
							 i == 2);}
  Vertex_handle v;
  if (dimension() == 1) {
    v = create_vertex();
    Face_handle ff = f->neighbor(0);
    Vertex_handle vv = f->vertex(1);
    Face_handle g = create_face(v,vv,NULL,ff, f, NULL);
    f->set_vertex(1,v);f->set_neighbor(0,g);
    ff->set_neighbor(1,g);
    v->set_face(g);
    vv->set_face(ff);
  }

    else { //dimension() ==2
    Face_handle n = f->neighbor(i);
    int in = f->mirror_index(i); //n->index(f);
    v = insert_in_face(f);
    flip(n,in); 
    }

  return v;
}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_dim_up(Vertex_handle w,  bool orient)
{
  // the following function insert 
  // a vertex  v which is outside the affine  hull of Tds
  // The triangulation will be starred from  v and w 
  // ( geometrically w=  // the infinite vertex )
  // w=NULL for first and second insertions
  // orient governs the orientation of the resulting triangulation

  Vertex_handle v = create_vertex();
  set_dimension(dimension() + 1);
  Face_handle f1;
  Face_handle f2;
    
  switch (dimension()) { //it is the resulting dimension
  case -1:
    f1 = create_face(v,Vertex_handle(NULL),Vertex_handle(NULL));
    v->set_face(f1);
    break;
  case 0 :
    f1 = &(*iterator_base_begin());
    f2 = create_face(v,Vertex_handle(NULL),Vertex_handle(NULL));
    f1->set_neighbor(0,f2);
    f2->set_neighbor(0,f1);
    v->set_face(f2);
    break;
  case 1 :
  case 2 :
    {
      std::list<Face_handle> faces_list;
      Iterator_base ib= iterator_base_begin(); 
      Iterator_base ib_end = iterator_base_end();
      for (; ib != ib_end ; ++ib){
	faces_list.push_back( & (*ib));
      }
      
      std::list<Face_handle>  to_delete;
      typename std::list<Face_handle>::iterator lfit = faces_list.begin();
      int i = dimension(); // maximun non NULL index in faces 
      Face_handle f, g;

      for ( ; lfit != faces_list.end() ; ++lfit) {
	f = * lfit;
	g = create_face(f); //calls copy constructor of face
	f->set_vertex(i,v); f->set_neighbor(i,g);
	g->set_vertex(i,w); g->set_neighbor(i,f);
	if (f->has_vertex(w)) to_delete.push_back(g); // flat face to delete
      }

      lfit = faces_list.begin();
      for ( ; lfit != faces_list.end() ; ++lfit) {
	f = * lfit;
	g = f->neighbor(i);
	for(int j = 0; j < i ; ++j) {
	  g->set_neighbor(j, f->neighbor(j)->neighbor(i));
	}
      }

      // couldn't unify the code for reorientation mater
      lfit = faces_list.begin() ; 
      if (dimension() == 1){
	if (orient) {
	  (*lfit)->reorient(); ++lfit ;  (*lfit)->neighbor(1)->reorient();
	}
	else {
	  (*lfit)->neighbor(1)->reorient(); ++lfit ; (*lfit)->reorient(); 
	}
      }
      else { // dimension == 2
	for( ;lfit  != faces_list.end(); ++lfit ) {
	  if (orient) {(*lfit)->neighbor(2)->reorient();}
	  else { (*lfit)->reorient();}
	}
      }

      lfit = to_delete.begin();
      int i1, i2;
      for ( ;lfit  != to_delete.end(); ++lfit){
	f = *lfit ;
	int j ;
	if (f->vertex(0) == w) {j=0;}
	else {j=1;}
	f1= f->neighbor(i); i1= f->mirror_index(i); //f1->index(f);
	f2= f->neighbor(j); i2= f->mirror_index(j); //f2->index(f);
	f1->set_neighbor(i1,f2);
	f2->set_neighbor(i2,f1);
	delete_face(f);
      }
    
      v->set_face( *(faces_list.begin()));
    }
    break;
  default:
    CGAL_triangulation_assertion(false);
    break;
  }
  return v;
}


template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
remove_degree_3(Vertex_handle v, Face_handle f)
// remove a vertex of degree 3
{
  CGAL_triangulation_precondition(v != NULL);
  CGAL_triangulation_precondition(v->degree() == 3);

  if (f == NULL) {f= v->face();}
  else { CGAL_triangulation_assertion( f->has_vertex(v));}
      
  int i = f->index(v);
  Face_handle left = f->neighbor(cw(i));
  int li = f->mirror_index(cw(i)); 
  Face_handle right = f->neighbor(ccw(i));
  int ri = f->mirror_index(ccw(i)); 

  Face_handle ll, rr;
  Vertex_handle q = left->vertex(li);
  CGAL_triangulation_assertion( left->vertex(li) == right->vertex(ri));
    
  ll = left->neighbor(cw(li));
  if(ll != NULL) {
    int lli = left->mirror_index(cw(li)); 
    ll->set_neighbor(lli, f);
  } 
  f->set_neighbor(cw(i), ll);
  if (f->vertex(ccw(i))->face() == left) f->vertex(ccw(i))->set_face(f);    
        
  rr = right->neighbor(ccw(ri));
  if(rr != NULL) {
    int rri =  right->mirror_index(ccw(ri)); //rr->index(right);
    rr->set_neighbor(rri, f);
  } 
  f->set_neighbor(ccw(i), rr);
  if (f->vertex(cw(i))->face() == right) f->vertex(cw(i))->set_face(f);  
        
  f->set_vertex(i, q);
  if (q->face() == right || q->face() == left) {
    q->set_face(f);
  }
  delete_face(right);
  delete_face(left);
        
  delete_vertex(v);
} 


  
template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
remove_dim_down(Vertex_handle v)
{
  Face_handle f;
  switch( dimension()){
  case -1: 
    delete_face(v->face());
    break;
  case 0:
    f = v->face();
    f->neighbor(0)->set_neighbor(0,NULL);
    delete_face(v->face());
    break;
  case 1:
  case 2:
    CGAL_triangulation_precondition ( 
	     (dimension() == 1 &&  number_of_vertices() == 3) ||
	     (dimension() == 2 && number_of_vertices() > 3) );
    // the faces incident to v are down graded one dimension
    // the other faces are deleted
    std::list<Face_handle > to_delete;
    std::list<Face_handle > to_downgrade;
    Iterator_base ib = iterator_base_begin();
    for( ; ib != iterator_base_end(); ++ib ){
      if ( ! ib->has_vertex(v) ) { to_delete.push_back(&(*ib));}
      else { to_downgrade.push_back(&(*ib));}
    }

    typename std::list<Face_handle>::iterator lfit = to_downgrade.begin();
    int j;
    for( ; lfit !=  to_downgrade.end() ; ++lfit) {
      f = *lfit; j = f->index(v);
      if (dimension() == 1) {
	if (j == 0) 	f->reorient();
	f->set_vertex(1,NULL);
	f->set_neighbor(1,NULL);
      }
      else { //dimension() == 2
	if (j == 0) f->cw_permute();
	else if(j == 1) f->ccw_permute();
	f->set_vertex(2,NULL);
	f->set_neighbor(2,NULL);
      }
      f->vertex(0)->set_face(f);
    }

    lfit = to_delete.begin();
    for( ; lfit !=  to_delete.end() ; ++lfit) {
      delete_face(*lfit);
    }
  }  
  delete_vertex(v);
  //set_number_of_vertices(number_of_vertices() -1);
  set_dimension(dimension() -1);
  return;
}

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::  
remove_1D(Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 1 &&
				   number_of_vertices() > 3);
  Face_handle f = v->face();
  int i = f->index(v);
  if (i==0) {f = f->neighbor(1);}
  CGAL_triangulation_assertion( f->index(v) == 1);
  Face_handle g= f->neighbor(0);
  f->set_vertex(1, g->vertex(1));
  f->set_neighbor(0,g->neighbor(0));
  g->neighbor(0)->set_neighbor(1,f);
  g->vertex(1)->set_face(f);
  delete_face(g);
  delete_vertex(v);
  return;
}



template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
remove_second(Vertex_handle v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 2 &&
 				  dimension() == 0);
  remove_dim_down(v);
  return;
}

    
template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
remove_first(Vertex_handle v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 1 && 
 				  dimension() == -1);
  remove_dim_down(v);
  return; 
}

template <class Vb, class Fb>
inline
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
star_hole(List_edges& hole)
{
  Vertex_handle newv = create_vertex();
  star_hole(newv, hole);
  return newv;
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
star_hole(Vertex_handle newv, List_edges& hole)
  // star the hole represented by hole around newv
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
{
   
  star_hole(newv, hole.begin(), hole.end());
  return;	    
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
make_hole(Vertex_handle v, List_edges& hole)
  // delete the faces incident to v and v
  // and return the dscription of the hole in hole
{
 CGAL_triangulation_precondition(dimension() == 2);
 std::list<Face_handle> to_delete;  

 Face_handle  f, fn;
 int i =0, in =0;
 Vertex_handle  vv;

 Face_circulator fc = v->incident_faces();
 Face_circulator done(fc);
 do {
   f = &(*fc);
   i = f->index(v);
   fn = f->neighbor(i);
   in = f->mirror_index(i); //fn->index(f);
   vv = f->vertex(cw(i));
   if( vv->face()==  f) vv->set_face(fn);
   vv = fc->vertex(ccw(i));
   if( vv->face()== f) vv->set_face(fn);
   fn->set_neighbor(in, NULL);
   hole.push_back(Edge(fn,in));
   to_delete.push_back(f);
 }
  while(++fc != done);

  while (! to_delete.empty()){
    delete_face(to_delete.front());
    to_delete.pop_front();
  }
  delete_vertex(v);
  return;
}


template <class Vb, class Fb>
inline
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
create_vertex()
{
  ++_number_of_vertices;
  return vertex_container().get_new_element();
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, 
	    Face_handle f2, int i2, 
	    Face_handle f3, int i3)
{
  Face_handle newf = create_face();
  newf->set_vertices(f1->vertex(cw(i1)),
		     f2->vertex(cw(i2)),
		     f3->vertex(cw(i3)));
  newf->set_neighbors(f2, f3, f1);
  // new Face(f1->vertex(cw(i1)),
// 			f2->vertex(cw(i2)),
// 			f3->vertex(cw(i3)),
// 			f2, f3, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  f3->set_neighbor(i3,newf);
  //  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, Face_handle f2, int i2)
{
  Face_handle newf = create_face();
  newf->set_vertices(f1->vertex(cw(i1)),
		     f2->vertex(cw(i2)),
		     f2->vertex(ccw(i2)));
  newf->set_neighbors(f2, NULL, f1);
//   new Face(f1->vertex(cw(i1)),
// 			f2->vertex(cw(i2)),
// 			f2->vertex(ccw(i2)),
// 			f2, NULL, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  //  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, Vertex_handle v)
{
  Face_handle newf = create_face();
  newf->set_vertices(f1->vertex(cw(i1)), f1->vertex(ccw(i1)), v);
  newf->set_neighbors(NULL, NULL, f1);
//     new Face(f1->vertex(cw(i1)), f1->vertex(ccw(i1)),v,
// 			NULL, NULL, f1);
  f1->set_neighbor(i1,newf);
  //  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face()
{
  Face_handle newf = face_container().get_new_element();
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face( Face_handle f)
{
  Face_handle newf = create_face();
  *newf = *f;
  //new Face(*f); //calls copy constructor
  //add_face(newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  Face_handle newf = create_face();
  newf->set_vertices(v1, v2, v3);
    //new Face(v1, v2, v3);
    //add_face(newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
	    Face_handle f1, Face_handle f2, Face_handle f3)
{
 Face_handle newf = create_face();
 newf->set_vertices(v1, v2, v3);
 newf->set_neighbors(f1, f2, f3);
 //new Face( v1, v2, v3, f1, f2, f3);
 //add_face(newf);
 return(newf);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
set_adjacency(Face_handle f0, int i0, Face_handle f1, int i1) const
{
  CGAL_triangulation_assertion(i0 >= 0 && i0 <= dimension());
  CGAL_triangulation_assertion(i1 >= 0 && i1 <= dimension());
  CGAL_triangulation_assertion(f0 != f1);
  f0->set_neighbor(i0,f1);
  f1->set_neighbor(i1,f0);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
delete_face(Face_handle f)
{
  CGAL_triangulation_expensive_precondition( face_container().is_element(&*f));
  face_container().release_element(&*f); 
//   f->previous()->set_next(f->next());
//   f->next()->set_previous(f->previous());
//   delete &*f;
}



template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
delete_vertex(Vertex_handle v)
{
  CGAL_triangulation_expensive_precondition( is_vertex(v) );
  --_number_of_vertices;
  vertex_container().release_element(&*v);
  // delete &*v;
//   set_number_of_vertices(number_of_vertices() -1);
//   return;
}

// CHECKING
template <  class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_valid(bool verbose, int level) const
{
  if(number_of_vertices() == 0){ 
    return (dimension() == -2);
  }

      
  bool result = (dimension()>= -1);
  CGAL_triangulation_assertion(result);

  //count and test the validity of the faces (for positive dimensions)
  Iterator_base ib = iterator_base_begin(); 
  Iterator_base ib_end = iterator_base_end();
  int count_stored_faces =0;
  for ( ; ib != ib_end ; ++ib){
    count_stored_faces += 1;
    if (dimension()>= 0) {
      result = result && ib->is_valid(verbose,level);
      CGAL_triangulation_assertion(result);
    }
  }
  
  result = result && (count_stored_faces == number_of_full_dim_faces());
  CGAL_triangulation_assertion(
		 count_stored_faces == number_of_full_dim_faces());
 
  // vertex count
  int vertex_count = 0;
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end();
      ++vit) {
    CGAL_triangulation_assertion( vit->face() != NULL);
    result = result && vit->is_valid(verbose,level);
    CGAL_triangulation_assertion( result );
    ++vertex_count;
  }
  result = result && (number_of_vertices() == vertex_count);
  CGAL_triangulation_assertion( number_of_vertices() == vertex_count );
    
  //edge count
  int edge_count = 0;
  for(Edge_iterator eit = edges_begin(); eit != edges_end(); ++eit) { 
    ++edge_count;
  }

  // face count
  int face_count = 0;
  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit) {
    ++face_count;
  }
        
  switch(dimension()) {
  case -1: 
    result = result && vertex_count == 1 && face_count == 0
      && edge_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 0:
    result = result && vertex_count == 2 && face_count == 0
      && edge_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 1:
    result = result &&  edge_count == vertex_count;
    CGAL_triangulation_assertion(result);
    result = result &&  face_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 2:
    result = result &&  edge_count == 3*face_count/2 ;
    CGAL_triangulation_assertion(edge_count == 3*face_count/2);
    break;
  default:
    result = false;
    CGAL_triangulation_assertion(result);
  }
  return result;
}



template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
copy_tds(const Tds &tds, Vertex_handle vh)
{
  if (this == &tds) return Vertex_handle(NULL);
  if (vh != NULL) 
    CGAL_triangulation_precondition( tds.is_vertex(vh));
  _number_of_vertices = tds.number_of_vertices();
  _dimension = tds.dimension();
  _face_container = tds.face_container();
  _vertex_container = tds.vertex_container();

  if(tds.number_of_vertices() == 0){return Vertex_handle(NULL);}
  //initializes maps
  std::map<Vertex_handle,Vertex_handle> vmap;
  std::map<Face_handle,Face_handle> fmap;
  Iterator_base it1 = tds.iterator_base_begin();
  Iterator_base it2 = iterator_base_begin();
  for( ; it1 != tds.iterator_base_end(); ++it1,++it2) {
    fmap[it1->handle()] = it2->handle();
  }
  Vertex_iterator vit1 = tds.vertices_begin();
  Vertex_iterator vit2 = vertices_begin();
  for( ; vit1 != tds.vertices_end(); ++vit1,++vit2) {
    vmap[vit1->handle()] = vit2->handle();
  }

  //update pointers
  it2 = iterator_base_begin();
  int j;
  if (dimension() == -1)  it2->set_vertex(0, vmap[it2->vertex(0)]);
  else {  // dimension() >= 0
    for ( ; it2 != iterator_base_end(); ++it2) {
      for (j = 0; j < dimension()+ 1; ++j) {
	it2->set_vertex(j, vmap[it2->vertex(j)]);
	it2->set_neighbor(j, fmap[it2->neighbor(j)]);
      }
    }
  }
  vit2 = vertices_begin();
  for ( ; vit2 != vertices_end(); vit2++) {
    vit2->set_face(fmap[vit2->face()]);
  }

  if (vh == Vertex_handle(NULL)) return Vertex_handle(NULL);
  return vmap[vh];
}
 

template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
file_output( std::ostream& os, Vertex_handle v, bool skip_first) const
{
  // ouput to a file
  // if non NULL, v is the vertex to be output first
  // if skip_first is true, the point in the first vertex is not output
  // (it may be for instance the infinite vertex of the triangulation)
  
  int n = number_of_vertices();
  int m = number_of_full_dim_faces();
  if(is_ascii(os))  os << n << ' ' << m << ' ' << dimension() << std::endl;
  else     os << n << m << dimension();
  if (n==0) return;

  std::map<Vertex_handle,int> V;
  std::map<Face_handle,int> F;

  // first vertex 
  int inum = 0;
  if ( v != NULL) {
    V[v] = inum++;
    if( ! skip_first){
    os << v->point();
    if(is_ascii(os))  os << std::endl;
    }
  }
  
  // other vertices
  for( Vertex_iterator vit= vertices_begin(); vit != vertices_end() ; ++vit) {
    if ( vit->handle() != v) {
	V[vit->handle()] = inum++;
	os << vit->point();
	if(is_ascii(os)) os << std::endl;
    }
  }
  //if(is_ascii(os)) os << "\n";

  // vertices of the faces
  inum = 0;
  int dim = (dimension() == -1 ? 1 :  dimension() + 1);
  for( Iterator_base ib = iterator_base_begin();
       ib != iterator_base_end(); ++ib) {
    F[&(*ib)] = inum++;
    for(int j = 0; j < dim ; ++j){
      os << V[ib->vertex(j)];
      if(is_ascii(os)){
	if(j== dimension())   os << "\n";
	else     os <<  ' ';
      }
    }
  }
    
  // neighbor pointers of the  faces
  for( Iterator_base it = iterator_base_begin();
       it != iterator_base_end(); ++it) {
    for(int j = 0; j < dimension()+1; ++j){
      os << F[&(*(it->neighbor(j)))];
      if(is_ascii(os)){
	if(j== dimension()) os << "\n";
	else    os <<  ' ';
      }
    }
  }

  return ;
}


template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
file_input( std::istream& is, bool skip_first)
{
  //input from file
  //return a pointer to the first input vertex
  // if skip_first is true, a first vertex is added (infinite_vertex)
  //set this  first vertex as infinite_Vertex
  if(number_of_vertices() != 0)    clear();
  
  int n, m, d;
  is >> n >> m >> d;

  if (n==0){ return NULL;}

  //set_number_of_vertices(n);
  set_dimension(d);

  std::vector<Vertex_handle > V(n);
  std::vector<Face_handle> F(m);

  // read vertices
  int i = 0;
  if(skip_first){
    V[0] = create_vertex();
    ++i;
  }
  for( ; i < n; ++i) {
    typename Vertex_base::Point p;
    is >> p;
    V[i] = create_vertex();
    V[i]->set_point(p);
  }
  
  // Creation of the faces
  int index;
  int dim = (dimension() == -1 ? 1 :  dimension() + 1);
  {
    for(i = 0; i < m; ++i) {
      F[i] = create_face() ;
      for(int j = 0; j < dim ; ++j){
	is >> index;
	F[i]->set_vertex(j, V[index]);
	// The face pointer of vertices is set too often,
	// but otherwise we had to use a further map
	V[index]->set_face(F[i]);
      }
    }
  }

  // Setting the neighbor pointers 
  {
    for(i = 0; i < m; ++i) {
      for(int j = 0; j < dimension()+1; ++j){
	is >> index;
	F[i]->set_neighbor(j, F[index]);
      }
    }
  }
  
  return V[0];
}


template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
vrml_output( std::ostream& os, Vertex_handle v, bool skip_infinite) const
{
  // ouput to a vrml file style
  // Point are assumed to be 3d points with a stream oprator <<
  // if non NULL, v is the vertex to be output first
  // if skip_inf is true, the point in the first vertex is not output
  // and the faces incident to v are not output
  // (it may be for instance the infinite vertex of the terrain)
  os << "#VRML V2.0 utf8" << std::endl;
  os << "Shape {" << std::endl;
  os << "\tgeometry IndexedFaceSet {" << std::endl;
  os << "\t\tcoord Coordinate {" << std::endl;
  os << "\t\t\tpoint [" << std::endl;

  std::map<Vertex_handle,int> vmap;
  Vertex_iterator vit;
  Face_iterator fit;

  //first vertex
  int inum = 0;
  if ( v != NULL) {
    vmap[v] = inum++;
    if( ! skip_infinite)  os << "\t\t\t\t" << vit->point() << std::endl;
  }

  //other vertices
  for( vit= vertices_begin(); vit != vertices_end() ; ++vit) {
    if ( vit->handle() != v) {
      vmap[vit->handle()] = inum++;
      os << "\t\t\t\t" << vit->point() << std::endl;
    }
  }

   os << "\t\t\t]" << std::endl;
   os << "\t\t}" << std::endl;
   os << "\t\tcoordIndex [" << std::endl;

   // faces
   for(fit= faces_begin(); fit != faces_end(); ++fit) {
     if (!skip_infinite || !fit->has_vertex(v)) {
   	os << "\t\t\t";
	os << vmap[(*fit).vertex(0)] << ", ";
	os << vmap[(*fit).vertex(1)] << ", ";
	os << vmap[(*fit).vertex(2)] << ", ";
	os << "-1, " << std::endl;  
     }
   }
   os << "\t\t]" << std::endl;
   os << "\t}" << std::endl;
   os << "}" << std::endl;
   return;
}

template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
off_file_input( std::istream& is, bool verbose)
{
  // input from an OFF file
  // assume a dimension 2 triangulation
  // create an infinite-vertex and  infinite faces with the
  // boundary edges if any.
  // return the infinite vertex if created
  Vertex_handle vinf(NULL);
  File_scanner_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr << " " << std::endl;
	 std::cerr << "TDS::off_file_input" << std::endl;
	 std::cerr << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }

  if(number_of_vertices() != 0)    clear();
  int dim = 2;
  set_dimension(dim);

  std::vector<Vertex_handle > vvh(scanner.size_of_vertices());
  std::map<Vh_pair, Edge> edge_map;
  typedef typename Vb::Point   Point;

  // read vertices
  int i;
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Point p;
    file_scan_vertex( scanner, p);
    vvh[i] = create_vertex();
    vvh[i]->set_point(p);
    scanner.skip_to_next_vertex( i);
  }
  if ( ! is ) {
    is.clear( std::ios::badbit);
    return vinf;
  }
  //vinf = vvh[0];

  // create the facets
  for ( i = 0; i < scanner.size_of_facets(); i++) {
    Face_handle fh = create_face();
    Integer32 no;
    scanner.scan_facet( no, i);
    if( ! is || no != 3) {
      if ( scanner.verbose()) {
	std::cerr << " " << std::endl;
	std::cerr << "TDS::off_file_input" << std::endl;
	std::cerr << "facet " << i << "does not have  3 vertices." 
		  << std::endl;
      }
      is.clear( std::ios::badbit);
      return vinf;
    }

    for ( int j = 0; j < no; ++j) {
      Integer32 index;
      scanner.scan_facet_vertex_index( index, i);
      fh->set_vertex(j, vvh[index]);
      vvh[index]->set_face(fh);
    }

    for (int ih  = 0; ih < no; ++ih) {
	set_adjacency(fh, ih, edge_map);
    }
  }

  // deal with  boundaries
  if ( !edge_map.empty()) {
    vinf = create_vertex();
    std::map<Vh_pair, Edge> inf_edge_map;
   while (!edge_map.empty()) {
     Face_handle fh = edge_map.begin()->second.first;
     int ih = edge_map.begin()->second.second;
     Face_handle fn = create_face( vinf, 
				   fh->vertex(cw(ih)), 
				   fh->vertex(ccw(ih)));
     vinf->set_face(fn);
     set_adjacency(fn, 0, fh, ih);
     set_adjacency(fn, 1, inf_edge_map);
     set_adjacency(fn, 2, inf_edge_map);
     edge_map.erase(edge_map.begin());
   }
   CGAL_triangulation_assertion(inf_edge_map.empty());
  }
  
  
  // coherent orientation
  reorient_faces();
  return vinf;
}


template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
set_adjacency(Face_handle fh, 
	      int ih, 
	      std::map< Vh_pair, Edge>& edge_map)
{
  // set adjacency to (fh,ih) using the the map edge_map
  // or insert (fh,ih) in edge map
  Vertex_handle vhcw  =  fh->vertex(cw(ih));
  Vertex_handle vhccw =  fh->vertex(ccw(ih)); 
  Vh_pair  vhp =  vhcw < vhccw ?  
                  std::make_pair(vhcw, vhccw) 
                : std::make_pair(vhccw, vhcw) ;
  typename std::map<Vh_pair, Edge>::iterator emapit = edge_map.find(vhp);
  if (emapit == edge_map.end()) {// not found, insert edge
    edge_map.insert(std::make_pair(vhp, Edge(fh,ih)));
  }
  else { //found set adjacency and erase
    Edge e = emapit->second;
    set_adjacency( fh,ih, e.first, e.second);
    edge_map.erase(emapit);
  } 
}



template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
reorient_faces()
{
  // reorient the faces of a triangulation 
  // needed for example in off_file_input
  // because the genus is not known, the number of faces 
  std::set<Face_handle> oriented_set;
  std::stack<Face_handle>  st;
  Face_iterator fit = faces_begin();
  int nf  = std::distance(faces_begin(),faces_end());

  while (static_cast<int>(oriented_set.size()) != nf) {
    while ( oriented_set.find(fit->handle()) != oriented_set.end()){
      ++fit; // find a germ for  non oriented components 
    }
    // orient component
    oriented_set.insert(fit->handle());
    st.push(fit->handle());
    while ( ! st.empty()) {
      Face_handle fh = st.top();
      st.pop();
      for(int ih = 0 ; ih < 3 ; ++ih){
	Face_handle fn = fh->neighbor(ih);
	if (oriented_set.find(fn) == oriented_set.end()){
	  int in = fn->index(fh);
	  if (fn->vertex(cw(in)) != fh->vertex(ccw(ih))) fn->reorient();
	  oriented_set.insert(fn);
	  st.push(fn);
	}
      }
    }

  }
  return;
}
	  

template <  class Vb, class Fb>
std::istream&
operator>>(std::istream& is,  
	   Triangulation_data_structure_2<Vb,Fb>& tds) 
{
  tds.file_input(is);
  return is;
}


template <  class Vb, class Fb>
std::ostream&
operator<<(std::ostream& os, 
	   const Triangulation_data_structure_2<Vb,Fb>  &tds) 
{
   tds.file_output(os);
   return os;
}

CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_DATA_STRUCTURE_2_H
