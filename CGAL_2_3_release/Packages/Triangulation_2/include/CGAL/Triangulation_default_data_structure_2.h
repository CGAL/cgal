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
// file          : include/CGAL/Triangulation_default_data_structure_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_DEFAULT_DATA_STRUCTURE_2_H
#define CGAL_TRIANGULATION_DEFAULT_DATA_STRUCTURE_2_H

#include <utility>
#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_ds_face_2.h>
#include <CGAL/Triangulation_ds_vertex_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>


CGAL_BEGIN_NAMESPACE 
template < class Gt , class Vb, class Fb>
class Triangulation_default_data_structure_2;

template < class Gt , class Vb, class Fb>
std::istream& operator>>
(std::istream& is, Triangulation_default_data_structure_2<Gt,Vb,Fb>&
  tds);

template < class Gt , class Vb, class Fb>
std::ostream& operator<< 
( std::ostream& os, 
  const Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds);
 
  


template < class Gt , class Vb, class Fb>
class Triangulation_default_data_structure_2 
  :public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
     ( std::istream& is, 
       Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds);
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
     ( std::ostream& os, 
      const Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds);
  friend class Triangulation_ds_iterator_base_2<
                 Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
  friend class Triangulation_ds_face_iterator_2<
                   Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
  friend class Triangulation_ds_edge_iterator_2<
                   Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
  friend class Triangulation_ds_vertex_iterator_2<
                   Triangulation_default_data_structure_2<Gt,Vb,Fb> >;

public:
  typedef Gt Geom_traits;
  typedef Vb  Vertex_base;
  typedef Fb  Face_base;

  typedef Triangulation_ds_vertex_2<Vertex_base,Face_base> Vertex;
  typedef Triangulation_ds_face_2<Vertex_base,Face_base> Face;
  typedef std::pair<Face*, int>  Edge;

  typedef Triangulation_default_data_structure_2<Gt,Vertex_base,Face_base> Tds;
  typedef Triangulation_ds_iterator_base_2<Tds> Iterator_base;
  typedef Triangulation_ds_face_iterator_2<Tds> Face_iterator;
  typedef Triangulation_ds_vertex_iterator_2<Tds> Vertex_iterator;
  typedef Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;

  typedef Triangulation_ds_face_circulator_2<Vertex,Face> 
							Face_circulator;
  typedef Triangulation_ds_vertex_circulator_2<Vertex,Face> 
							Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_2<Vertex,Face> 
							Edge_circulator;
  typedef std::list<Edge> List_edges;

  
private:
  Geom_traits _geom_traits;
  Vertex* _infinite_vertex; // this is for the iterator only
  int _number_of_vertices; 
  int _dimension;

public:
  //CREATORS - DESTRUCTORS
  Triangulation_default_data_structure_2(const Geom_traits& gt=Geom_traits()); 
  Triangulation_default_data_structure_2(const Tds &tds);
  ~Triangulation_default_data_structure_2();
  Tds& operator= (const Tds &tds);
  void swap(Tds &tds);
  void clear();

  //ACCESS FUNCTIONS
  int  dimension() const { return _dimension;  }
  int number_of_vertices() const {return _number_of_vertices;}
  int number_of_faces() const ;
  int number_of_edges() const;
  int number_of_full_dim_faces() const; //number of faces stored by tds
  const Geom_traits& geom_traits() const {return _geom_traits;}

  // TEST FEATURES
  bool is_vertex(const Vertex* v) const;
  bool is_edge(const Vertex* va, const Vertex* vb) const;
  bool is_edge(const Vertex* va, const Vertex* vb, Face* &fr,  int &i) const;
  bool is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3) const;
  bool is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3,
      Face* &fr) const;

  // ITERATORS AND CIRCULATORS
  inline Iterator_base iterator_base_begin() const    {
    return Iterator_base(this);
  }
  inline Iterator_base iterator_base_end() const    {
    return Iterator_base(this,1);
  }
  inline Face_iterator faces_begin() const {
    Tds* ncthis = (Tds *)this;
    return Face_iterator(ncthis);
  }
  inline Face_iterator faces_end() const {
    Tds* ncthis = (Tds *)this;
    return Face_iterator(ncthis, 1);
  }
  inline Vertex_iterator vertices_begin() const  {
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis);
  }
  inline Vertex_iterator vertices_end() const {
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis,1);
  }
  inline Edge_iterator edges_begin() const {
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis);
  }
  inline Edge_iterator edges_end() const {
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis,1);
  }
  
  //Face_circulator incident_faces(Vertex* v, Face* f = NULL) const;
  //Vertex_circulator incident_vertices(Vertex* v, Face* f = NULL) const;
  //Edge_circulator incident_edges(Vertex* v, Face* f = NULL) const;

  // MODIFY
  void flip(Face* f, int i);
 
  Vertex* insert_first();
  Vertex* insert_second();
  Vertex* insert_in_face(Face* f);
  Vertex* insert_in_edge(Face* f, int i);
  Vertex* insert_dim_up(Vertex *w, bool orient=true);
  
  void remove_degree_3(Vertex* v, Face* f = NULL);
  void remove_1D(Vertex* v); 
   
  void remove_second(Vertex* v);
  void remove_first(Vertex* v);
  void remove_dim_down(Vertex* v);

  Vertex* star_hole(List_edges& hole);
  void    star_hole(Vertex* v, List_edges& hole);
  void    make_hole(Vertex* v, List_edges& hole);

//   template< class EdgeIt>
//   Vertex* star_hole(EdgeIt edge_begin,EdgeIt edge_end);
 
//   template< class EdgeIt>
//   void  star_hole(Vertex* v, EdgeIt edge_begin,  EdgeIt edge_end);

//   template< class EdgeIt, class FaceIt>
//   Vertex* star_hole(EdgeIt edge_begin, 
// 		    EdgeIt edge_end,
// 		    FaceIt face_begin,
// 		    FaceIt face_end);
 
//   template< class EdgeIt, class FaceIt>
//   void  star_hole(Vertex* v,
// 		  EdgeIt edge_begin, 
// 		  EdgeIt edge_end,
// 		  FaceIt face_begin,
// 		  FaceIt face_end);


  Vertex* create_vertex();
  Face* create_face(Face* f1, int i1, Face* f2, int i2, Face* f3, int i3);
  Face* create_face(Face* f1, int i1, Face* f2, int i2);
  Face* create_face(Face* f1, int i1, Vertex* v);
  Face* create_face(Vertex* v1, Vertex* v2, Vertex* v3);
  Face* create_face(Vertex* v1, Vertex* v2, Vertex* v3,
		    Face* f1, Face* f2, Face* f3);
  Face* create_face(Face* f); // calls copy constructor of Face
  Face* create_face();
  void  delete_face(Face*);
  void  delete_vertex(Vertex*);

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  
  // HELPING
  void copy_tds(const Tds &tds);
  Vertex* copy_tds(const Tds &tds, const Vertex*);

  Vertex* file_input(std::istream& is, bool skip_first=false);
  void file_output(std::ostream& os,
		   Vertex* v = NULL,
		   bool skip_first=false) const;


  // SETTING (had to make them public for use in remove from Triangulations)
  void set_number_of_vertices(int n) {_number_of_vertices = n;}
  void set_dimension (int n) {_dimension = n ;}
  void set_infinite_vertex(Vertex*  v) { _infinite_vertex = v;}

private:
   // INFINITE FEATURES
  Vertex* infinite_vertex() const  {return _infinite_vertex;  }
  Face* infinite_face() const { return _infinite_vertex->face();}
  bool is_infinite(const Face* f) const;
  bool is_infinite(const Vertex* v) const;
 
  // template members definition
public:
  template< class EdgeIt>
  Vertex* star_hole(EdgeIt edge_begin, EdgeIt edge_end)
  // creates a new vertex 
  // and stars from it
  // the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
  {
     Vertex* newv = create_vertex();
     star_hole(newv, edge_begin, edge_end);
     return newv;
  }
 
  template< class EdgeIt>
  void  star_hole(Vertex* v, EdgeIt edge_begin,  EdgeIt edge_end)
  // uses vertex v
  // to star the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // the hole is supposed to be ccw oriented
  { 
    std::list<Face*> empty_list;
    return star_hole(v, 
		     edge_begin, 
		     edge_end, 
		     empty_list.begin(),
		     empty_list.end());
    
  }


  template< class EdgeIt, class FaceIt>
  Vertex* star_hole(EdgeIt edge_begin, 
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
    Vertex* newv = create_vertex();
    star_hole(newv, edge_begin, edge_end, face_begin, face_end);
    return newv;
  }
 
  template< class EdgeIt, class FaceIt>
  void  star_hole(Vertex* newv,
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

    Face* first_f =  reset_or_create_face((*eit).first, (*eit).second, newv,
					  fit, face_end);
    ++eit;
    Face* previous_f=first_f, *next_f;
    for( ; eit != edge_end ; eit++) {
    next_f = reset_or_create_face((*eit).first, (*eit).second, newv,
				  fit, face_end);
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
  Face*  reset_or_create_face(Face* fn, 
			      int in, 
			      Vertex* v,
			      FaceIt fit,
			      FaceIt face_end) {
    if (fit == face_end) return create_face(fn, in, v);
    (*fit)->set_vertices(fn->vertex(cw(in)), fn->vertex(ccw(in)), v);
    (*fit)->set_neighbors(0,0,fn);
    return &(**fit++);    
  }

};

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb> ::
Triangulation_default_data_structure_2(const Geom_traits& gt) 
  :  _geom_traits(gt),_infinite_vertex(NULL),
     _number_of_vertices(0),_dimension(-1)
{ }

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb> ::
Triangulation_default_data_structure_2(const Tds &tds)
{
  copy_tds(tds);
}

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb> ::
~Triangulation_default_data_structure_2()
{
  clear();
}

//assignement  
template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>&
Triangulation_default_data_structure_2<Gt,Vb,Fb> ::
operator= (const Tds &tds)
{
  copy_tds(tds);
  return *this;
}  

///ACCESS FUNCTIONS
template < class Gt , class Vb, class Fb>
inline int 
Triangulation_default_data_structure_2<Gt,Vb,Fb> ::
number_of_faces() const 
{
  return ( dimension() < 2) ? 0 : (2*number_of_vertices()- 4);
}

template < class Gt , class Vb, class Fb>
int
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
number_of_edges() const
{
  switch (dimension()) {
  case 1:  return number_of_vertices();
  case 2:  return 3*number_of_vertices()-6;
  default: return 0;
  }
}
      
template < class Gt , class Vb, class Fb>
int
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
number_of_full_dim_faces() const
{
  switch (dimension()) {
  case -1: return 0;
  case 0:  return number_of_vertices();
  case 1:  return number_of_edges();
  case 2:  return number_of_faces();
  default: return 0;
  }
}

template < class Gt , class Vb, class Fb>
bool
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_vertex(const Vertex* v) const
{
  if (v == infinite_vertex()) return true; //short cut for frequent  case
  if (number_of_vertices() == 0) return false;
  for (Vertex_iterator vit = vertices_begin();
                       vit != vertices_end(); ++vit) {
   if ( v == &(*vit)) return true;
  }
  return false;
}

template <class Gt , class Vb, class Fb>
bool
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_edge(const Vertex* va, const Vertex* vb) const
// returns true (false) if the line segment ab is (is not) an edge of t
{
  Vertex_circulator vc= va->incident_vertices(), done(vc);
  if ( vc == 0) return false;
  do {
    if( vb == &(*vc) ) {return true;} 
  } while (++vc != done);
  return false;
}
 

template <class Gt , class Vb, class Fb>
bool
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_edge(const Vertex* va, const Vertex* vb, Face* &fr,  int & i) const
// returns true (false) if the line segment ab is (is not) an edge of t
// (fr,i) is the edge ab
  // with face fr on the right of a->b
{
  Face* fc=va->face();
  Face* start = fc;
  if (fc == NULL) return false;
  int inda, indb;
  do {
    inda = fc->index(va);
    indb = (dimension() == 2 ? cw(inda) : 1-inda);
    if(fc->vertex(indb) == vb) {
      fr=fc;
      i = 3 - inda - indb ; //works in dim 1 or 2
      return true;
    }
    fc=fc->neighbor(indb); //turns ccw around va
  } while (fc != start);
  return false;
}

template <class Gt , class Vb, class Fb>
inline bool 
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3) const
{
  Face* f;
  return is_face(v1,v2,v3,f);
}

template <class Gt , class Vb, class Fb>
bool
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3,
      Face* &f) const
{
  if (dimension() != 2) return false;
  int i;
  bool b = is_edge(v1,v2,f,i);
  if (!b) return false;
  else if (v3== f->vertex(i)) return true;
  f = f-> neighbor(i);
  if (v3 == f->vertex(3 - f->index(v1)- f->index(v2))) {return true;}
  return false;  
}



template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
flip(Face* f, int i)
{
  CGAL_triangulation_precondition( dimension()==2);
  Face* n  = f->neighbor(i);
  int ni = n->index(f);
    
  Vertex*  v_cw = f->vertex(cw(i));
  Vertex*  v_ccw = f->vertex(ccw(i));

  // bl == bottom left, tr == top right
  Face* tr = f->neighbor(ccw(i));
  Face* bl = n->neighbor(ccw(ni));
  int bli, tri;
  bli = bl->index(n);
  tri = tr->index(f);
    
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
  
template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
insert_first( )
{
  CGAL_triangulation_precondition( number_of_vertices() == 0 &&
				   dimension()==-1 );
  Vertex* v1 = create_vertex();
  set_infinite_vertex(v1);
  return v1;
}

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex* 
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
insert_second()
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 &&
				   dimension()==-1 );
  Vertex* v2= create_vertex();
  Face* f1 = create_face( _infinite_vertex, NULL, NULL);
  Face* f2 = create_face( v2, NULL, NULL);
  f1->set_neighbor(0,f2);
  f2->set_neighbor(0,f1);
  _infinite_vertex->set_face(f1);
  v2->set_face(f2);
  set_dimension(0);
  return v2;
}


template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
insert_in_face(Face* f)
  // New vertex will replace f->vertex(0) in face f
{
  CGAL_triangulation_precondition( f != NULL && dimension()== 2);
  Vertex*  v = create_vertex();
  Vertex* v0 = f->vertex(0);
  Vertex* v2 = f->vertex(2);
  Vertex* v1 = f->vertex(1);
    
  Face* n1 = f->neighbor(1);
  Face* n2 = f->neighbor(2);
    
  Face* f1 = create_face(v0, v, v2, f, n1, NULL);
  Face* f2 = create_face(v0, v1, v, f, NULL, n2);

  f1->set_neighbor(2, f2);
  f2->set_neighbor(1, f1);
  if (n1 != NULL) {
    int i1 = n1->index(f);
    n1->set_neighbor(i1,f1);
  }
  if (n2 != NULL) {
    int i2 = n2->index(f);
    n2->set_neighbor(i2,f2);}

  f->set_vertex(0, v);
  f->set_neighbor(1, f1);
  f->set_neighbor(2, f2);

  if( v0->face() == f  ) {  v0->set_face(f2); }
  v->set_face(f);

  return v;
}


template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
insert_in_edge(Face* f, int i)
  //insert in the edge opposite to vertex i of face f
{
  CGAL_triangulation_precondition(f != NULL && dimension() >= 1); 
  if (dimension() == 1) {CGAL_triangulation_precondition(i == 2);}
  if (dimension() == 2) {CGAL_triangulation_precondition(i == 0 || 
							 i == 1 || 
							 i == 2);}
  Vertex * v;
  if (dimension() == 1) {
    v = create_vertex();
    Face * ff = f->neighbor(0);
    Vertex * vv = f->vertex(1);
    Face * g = create_face(v,vv,NULL,ff, f, NULL);
    f->set_vertex(1,v);f->set_neighbor(0,g);
    ff->set_neighbor(1,g);
    v->set_face(g);
    vv->set_face(ff);
  }

    else { //dimension() ==2
    Face* n = f->neighbor(i);
    int in = n->index(f);
    v = insert_in_face(f);
    flip(n,in); 
    }

  return v;
}


template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
insert_dim_up(Vertex *w, bool orient)
{
  // the following function insert 
  // a vertex  v which is outside the affine  hull of a 1 dim or 0 dim Tds
  // w is the infinite vertex of the triangulation
  // orient governs the orientation of the resulting triangulation
  CGAL_triangulation_precondition(dimension()==0 || dimension()== 1);

  //Vertex * v = new Vertex;
  Vertex* v=create_vertex();

  std::list<Face *> faces_list;
  Iterator_base ib= Iterator_base(this); 
  Iterator_base ib_end = Iterator_base(this,1);
  for ( ; ib != ib_end ; ++ib){
    faces_list.push_back( & (*ib));
  }

  std::list<Face *>  to_delete;
  typename std::list<Face *>::iterator lfit = faces_list.begin();
  int i = dimension()+1; // maximun non NULL index in faces after the insertion
  Face *f, *g;

  for ( ; lfit != faces_list.end() ; ++lfit) {
    f = * lfit;
    g = create_face(f);
    f->set_vertex(i,v); f->set_neighbor(i,g);
    g->set_vertex(i,w); g->set_neighbor(i,f);
    if (f->has_vertex(w)) to_delete.push_back(g); // flat face to be deleted 
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
  if (dimension() == 0){
    if (orient) {
      (*lfit)->reorient(); ++lfit ;  (*lfit)->neighbor(1)->reorient();
    }
    else {
      (*lfit)->neighbor(1)->reorient(); ++lfit ; (*lfit)->reorient(); 
    }
  }
  else { // dimension == 1
    for( ;lfit  != faces_list.end(); ++lfit ) {
     if (orient) {(*lfit)->neighbor(2)->reorient();}
      else { (*lfit)->reorient();}
     }
  }

  lfit = to_delete.begin();
  Face *f1, *f2;
  int i1, i2;
  for ( ;lfit  != to_delete.end(); ++lfit){
    f = *lfit ;
    int j ;
    if (f->vertex(0) == w) {j=0;}
    else {j=1;}
    f1= f->neighbor(i); i1= f1->index(f);
    f2= f->neighbor(j); i2 = f2->index(f);
    f1->set_neighbor(i1,f2);
    f2->set_neighbor(i2,f1);
    delete f;
  }
    
  v->set_face( *(faces_list.begin()));
  set_dimension(dimension() + 1);
  return v;
}


template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
remove_degree_3(Vertex* v, Face* f)
   // remove a vertex of degree 3
  {
    CGAL_triangulation_precondition(v != NULL);
    CGAL_triangulation_precondition(v != _infinite_vertex);
    CGAL_triangulation_precondition(v->degree() == 3);

    if (f == NULL) {f= v->face();}
    else { CGAL_triangulation_assertion( f->has_vertex(v));}
      
    int i = f->index(v);
    Face* left = f->neighbor(cw(i));
    Face* right = f->neighbor(ccw(i));
    Face *ll, *rr;
        
    int li = left->index(f);
    int ri = right->index(f);
    Vertex* q = left->vertex(li);
    CGAL_triangulation_assertion( left->vertex(li) == right->vertex(ri));
    
    ll = left->neighbor(cw(li));
    if(ll != NULL) {
      int lli = ll->index(left);
      ll->set_neighbor(lli, f);
    } 
    f->set_neighbor(cw(i), ll);
    if (f->vertex(ccw(i))->face() == left) f->vertex(ccw(i))->set_face(f);    
        
    rr = right->neighbor(ccw(ri));
    if(rr != NULL) {
      int rri = rr->index(right);
      rr->set_neighbor(rri, f);
    } 
    f->set_neighbor(ccw(i), rr);
    if (f->vertex(cw(i))->face() == right) f->vertex(cw(i))->set_face(f);  
        
    f->set_vertex(i, q);
    if (q->face() == right || q->face() == left) {
	   q->set_face(f);
    }
    delete right;
    delete left;
        
    delete v;
    set_number_of_vertices( number_of_vertices() -1);
  } 


  
template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
remove_dim_down(Vertex* v)
{

  CGAL_triangulation_precondition ( 
			    (dimension() == 1 &&  number_of_vertices() == 3) ||
			    (dimension() == 2 && number_of_vertices() > 3) );
  // the faces incident to v are down graded one dimension
  // the other faces are deleted
  std::list<Face* > to_delete;
  std::list<Face* > to_downgrade;
  Iterator_base ib = iterator_base_begin();
  for( ; ib != iterator_base_end(); ++ib ){
    if ( ! ib->has_vertex(v) ) { to_delete.push_back(&(*ib));}
    else { to_downgrade.push_back(&(*ib));}
  }

  typename std::list<Face*>::iterator lfit = to_downgrade.begin();
  int j;
  Face * f;
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
    delete *lfit;
  }

    
  delete v;
  set_number_of_vertices(number_of_vertices() -1);
  set_dimension(dimension() -1);
  return;
}

template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::  
remove_1D(Vertex* v)
{
  CGAL_triangulation_precondition( dimension() == 1 &&
				   number_of_vertices() > 3);
  Face* f = v->face();
  int i = f->index(v);
  if (i==0) {f = f->neighbor(1);}
  CGAL_triangulation_assertion( f->index(v) == 1);
  Face* g= f->neighbor(0);
  f->set_vertex(1, g->vertex(1));
  f->set_neighbor(0,g->neighbor(0));
  g->neighbor(0)->set_neighbor(1,f);
  g->vertex(1)->set_face(f);
  delete g;
  delete v;
  set_number_of_vertices(number_of_vertices() -1);
  return;
}



template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
remove_second(Vertex* v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 2 &&
				  dimension() == 0);
  delete v->face()->neighbor(0);
  delete v->face();
  delete v;
  set_number_of_vertices(1);
  set_dimension(-1);
  infinite_vertex()->set_face(NULL);
  return;
}

    
template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
remove_first(Vertex* v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 1 && 
				  dimension() == -1 &&
				  v== infinite_vertex());
  delete v;
  set_number_of_vertices(0);
  set_infinite_vertex(NULL);
  return;
}


template <class Gt ,class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
star_hole(List_edges& hole)
{
  Vertex* newv = create_vertex();
  star_hole(newv, hole);
  return newv;
}

template <class Gt ,class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
star_hole(Vertex* newv, List_edges& hole)
  // star the hole represented by hole around newv
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
{
  CGAL_triangulation_precondition(dimension() == 2);
  typedef typename  List_edges::const_iterator Hole_it;
  Hole_it hit = hole.begin();

  Face* first_f = create_face(hit->first, hit->second, newv);
  ++hit;
  Face* previous_f=first_f, *next_f;
  for( ; hit != hole.end(); hit++) {
    next_f = create_face(hit->first, hit->second, newv);
    next_f->set_neighbor(1, previous_f);
    previous_f->set_neighbor(0, next_f);
    previous_f=next_f;
  }
  next_f->set_neighbor(0, first_f);
  first_f->set_neighbor(1, next_f);
  newv->set_face(first_f);
  return;    
}

template <class Gt ,class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
make_hole(Vertex* v, List_edges& hole)
  // delete the faces incident to v and v
  // and return the dscription of the hole in hole
{
 CGAL_triangulation_precondition(dimension() == 2);
 std::list<Face*> to_delete;  

 Face*  f, *fn;
 int i =0, in =0;
 Vertex*  vv;

 Face_circulator fc = v->incident_faces();
 Face_circulator done(fc);
 do {
   f = &(*fc);
   i = f->index(v);
   fn = f->neighbor(i);
   in = fn->index(f);
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

template <class Gt ,class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_vertex()
{
  Vertex* newv= new Vertex();
  set_number_of_vertices(number_of_vertices()+ 1);
  return newv;
}


template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Face* f1, int i1, Face* f2, int i2, Face* f3, int i3)
{
  Face* newf = new Face(f1->vertex(cw(i1)),
			f2->vertex(cw(i2)),
			f3->vertex(cw(i3)),
			f2, f3, f1);

  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  f3->set_neighbor(i3,newf);
  return newf;
}

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Face* f1, int i1, Face* f2, int i2)
{
  Face* newf = new Face(f1->vertex(cw(i1)),
			f2->vertex(cw(i2)),
			f2->vertex(ccw(i2)),
			f2, NULL, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  return newf;
}

template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Face* f1, int i1, Vertex* v)
{
  Face* newf = new Face(f1->vertex(cw(i1)), f1->vertex(ccw(i1)),v,
			NULL, NULL, f1);
  f1->set_neighbor(i1,newf);
  return newf;
}


template <class Gt, class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Face* f)
{
 Face* newf = new Face(*f);
 return newf;
}

template <class Gt, class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face()
{
 Face* newf= new Face();
 return newf;
}

template <class Gt, class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Vertex* v1, Vertex* v2, Vertex* v3)
{
 Face* newf= new Face(v1, v2, v3);
 return newf;
}

template <class Gt, class Vb, class Fb>
inline
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Face*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
create_face(Vertex* v1, Vertex* v2, Vertex* v3,
	    Face* f1, Face* f2, Face* f3)
{
  Face* newf= new Face( v1, v2, v3, f1, f2, f3);
  return(newf);
}

template <class Gt, class Vb, class Fb>
inline void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
delete_face(Face* f)
{
  delete f;
}

template <class Gt, class Vb, class Fb>
inline void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
delete_vertex(Vertex* v)
{
  delete v;
  set_number_of_vertices(number_of_vertices() -1);
  return;
}

// CHECKING
template < class Gt , class Vb, class Fb>
bool
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_valid(bool verbose, int level) const
{
  if(number_of_vertices() == 0){ 
    return (dimension() == -1  && infinite_vertex() == NULL);
  }

  if(number_of_vertices() == 1) {
    return (dimension() == -1 && infinite_vertex() != NULL );
  }
    
  bool result = (dimension()>= 0);
  CGAL_triangulation_assertion(result);

  //test the validity of the faces whatever dimension() is
  Iterator_base ib = Iterator_base(this); 
  Iterator_base ib_end = Iterator_base(this,1);
  int count_stored_faces =0;
  for ( ; ib != ib_end ; ++ib){
    count_stored_faces += 1;
    result = result && ib->is_valid(verbose,level);
    CGAL_triangulation_assertion(result);
  }
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
    result = result &&  edge_count == 3*vertex_count - 6 ;
    CGAL_triangulation_assertion(result);
    result = result && ( face_count == 2*vertex_count - 4 );
    CGAL_triangulation_assertion( face_count == 2*vertex_count - 4 );
    break;
  }
  return result;
}


template < class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
copy_tds(const Tds &tds, const Vertex* v)
{
  CGAL_triangulation_precondition( tds.is_vertex(v));
  _number_of_vertices = tds.number_of_vertices();
  _geom_traits = tds.geom_traits();
  //the class Geom_traits is required to have a pertinent operator=
  _dimension = tds.dimension();

  //if(tds.number_of_vertices() == 0){return static_cast<Vertex*>(NULL);}

  std::map< const void*, void*, std::less<const void*> > V;
  std::map< const void*, void*, std::less<const void*> > F;

  Vertex*  v2,* v_inf;
  Face* f2;

  // create the vertices
  for( Vertex_iterator it=tds.vertices_begin();
       it != tds.vertices_end(); ++it) {
    V[&(*it)] = new Vertex( *it );
  }
  V[0] = 0 ; //to cope with lower dimensional cases  when creating faces

  //set infinite_vertex
  v_inf = static_cast<Vertex*>(V[v]);
  set_infinite_vertex(v_inf);
  
  // create the faces
  for(Iterator_base ib = tds.iterator_base_begin();
      ib != tds.iterator_base_end(); ++ib) {
    f2 = create_face(&(*ib));
    F[&(*ib)]=  f2;
    f2->set_vertices((Vertex*) V[ib->vertex(0)],
		     (Vertex*) V[ib->vertex(1)],
		     (Vertex*) V[ib->vertex(2)] );
  }

  // link each vertex to a face
  for( Vertex_iterator vit = tds.vertices_begin();
       vit != tds.vertices_end() ; ++vit) {
    v2 = (Vertex*) V[&(*vit)];
    v2->set_face( (Face*) F[vit->face()] );
  }

  // hook neighbor of the  faces
 for( Iterator_base ibb = tds.iterator_base_begin();
      ibb != tds.iterator_base_end(); ++ibb){
   for(int j = 0; j <= tds.dimension(); ++j){
     f2 = (Face*) F[&(*ibb)];
     f2->set_neighbor(j, (Face*) F[ibb->neighbor(j)] );
   }
 }

  CGAL_triangulation_postcondition( is_valid() );
  return v_inf;
}

template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
copy_tds(const Tds &tds)
{
  _number_of_vertices = tds.number_of_vertices();
  _geom_traits = tds.geom_traits();
  _dimension = tds.dimension();

  if(tds.number_of_vertices() == 0){return;}
  copy_tds(tds, tds.infinite_vertex());
}
 
template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
swap(Tds &tds)
{
  Geom_traits   t  = geom_traits();
  Vertex*  iv = infinite_vertex();
  int     nv = number_of_vertices();
  int    dim = dimension();

  _geom_traits = tds.geom_traits(); 
  //the class Geom_traits is required to have a pertinent operator=
  set_infinite_vertex(tds.infinite_vertex());
  set_number_of_vertices(tds.number_of_vertices());
  set_dimension(tds.dimension());

  tds._geom_traits = t;
  tds.set_infinite_vertex(iv);
  tds.set_number_of_vertices(nv);
  tds.set_dimension(dim);
}

template < class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
clear()
{
  if(number_of_vertices() == 0) return;
  if(number_of_vertices() ==1) delete infinite_vertex();
  else{
    std::list<Face*> Faces;
    std::list<Vertex*> Vertices;

    for(Vertex_iterator it = vertices_begin();
	it != vertices_end(); ++it) {
      Vertices.push_front(&(*it));
    }
    for( Iterator_base ib = iterator_base_begin();
	 ib != iterator_base_end(); ++ib) {
      Faces.push_front(&(*ib));
    }
  
    for(typename std::list<Face*>::iterator lfit=Faces.begin();
	lfit != Faces.end(); ++lfit) {
      delete *lfit;
    }
    for( typename std::list<Vertex*>::iterator lvit=Vertices.begin();
	 lvit != Vertices.end(); ++lvit) {
      delete *lvit;
    }
  }  
  set_infinite_vertex(NULL);
  set_number_of_vertices(0);
  set_dimension(-1);
  return;
}

template <class Gt , class Vb, class Fb>
void
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
file_output( std::ostream& os, Vertex* v, bool skip_first) const
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

  std::map< void*, int, std::less<void*> > V;
  std::map< void*, int, std::less<void*> > F;

  // first vertex 
  int inum = 0;
  if ( v != NULL) {
    V[v] = inum++;
    if( ! skip_first){
    os << v->point();
    if(is_ascii(os))  os << ' ';
    }
  }
  
  // other vertices
  for( Vertex_iterator vit= vertices_begin(); vit != vertices_end() ; ++vit) {
    if ( &(*vit) != v) {
	V[&(*vit)] = inum++;
	os << vit->point();
	if(is_ascii(os)) os << ' ';
    }
  }
  if(is_ascii(os)) os << "\n";

  // vertices of the faces
  inum = 0;
  for( Iterator_base ib = iterator_base_begin();
       ib != iterator_base_end(); ++ib) {
    F[&(*ib)] = inum++;
    for(int j = 0; j < dimension()+1; ++j){
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


template <class Gt , class Vb, class Fb>
Triangulation_default_data_structure_2<Gt,Vb,Fb>::Vertex*
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
file_input( std::istream& is, bool skip_first)
{
  //input from file
  //return a pointer to the first input vertex
  // if no_first is true, a first vertex is added (infinite_vertex)
  //set this  first vertex as infinite_Vertex
  if(number_of_vertices() != 0)    clear();
  
  int n, m, d;
  is >> n >> m >> d;

  if (n==0){ return NULL;}

  set_number_of_vertices(n);
  set_dimension(d);

  std::vector<Vertex* > V(n);
  std::vector<Face*> F(m);

  // read vertices
  int i = 0;
  if(skip_first){
    V[0] = new Vertex();
    ++i;
  }
  for( ; i < n; ++i) {
    typename Vertex_base::Point p;
    is >> p;
    V[i] = new Vertex(p);
  }
  if(n>0)   set_infinite_vertex(V[0]);

  // Creation of the faces
  int index;
  {
    for(i = 0; i < m; ++i) {
      F[i] = new Face() ;
      for(int j = 0; j < dimension()+1; ++j){
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

template <class Gt , class Vb, class Fb>
inline bool 
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_infinite(const Face* f) const
{
  return f->has_vertex(infinite_vertex()); 
}

template <class Gt , class Vb, class Fb>
inline bool 
Triangulation_default_data_structure_2<Gt,Vb,Fb>::
is_infinite(const Vertex* v) const 
{
  return v == infinite_vertex();
}


template < class Gt , class Vb, class Fb>
std::istream&
operator>>(std::istream& is,  
	   Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds) 
{
  tds.file_input(is);
  return is;
}


template < class Gt, class Vb, class Fb>
std::ostream&
operator<<(std::ostream& os, 
	   const Triangulation_default_data_structure_2<Gt,Vb,Fb>  &tds) 
{
   tds.file_output(os, tds.infinite_vertex());
   return os;
}

CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_DEFAULT_DATA_STRUCTURE_2_H
