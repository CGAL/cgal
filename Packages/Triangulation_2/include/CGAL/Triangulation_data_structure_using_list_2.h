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
// file          : include/CGAL/Triangulation_data_structure_using_list_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_USING_LIST_2_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_USING_LIST_2_H

#include <CGAL/basic.h>
#include <iostream>
#include <list>
#include <map>
#include <vector>


#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_dsul_face_2.h>
#include <CGAL/Triangulation_dsul_vertex_2.h>
#include <CGAL/Triangulation_dsul_iterators_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>


CGAL_BEGIN_NAMESPACE 
template < class Vb, class Fb>
class Triangulation_data_structure_using_list_2;

template < class Vb, class Fb>
std::istream& operator>>
(std::istream& is, Triangulation_data_structure_using_list_2<Vb,Fb>&
  tds);

template < class Vb, class Fb>
std::ostream& operator<< 
( std::ostream& os, 
  const Triangulation_data_structure_using_list_2<Vb,Fb>& tds);
 
  


template <class Vb, class Fb>
class Triangulation_data_structure_using_list_2 
  :public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
     ( std::istream& is, 
       Triangulation_data_structure_using_list_2<Vb,Fb>& tds);
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
     ( std::ostream& os, 
      const Triangulation_data_structure_using_list_2<Vb,Fb>& tds);
  friend class Triangulation_dsul_iterator_base_2<
                 Triangulation_data_structure_using_list_2<Vb,Fb> >;
  friend class Triangulation_dsul_face_iterator_2<
                   Triangulation_data_structure_using_list_2<Vb,Fb> >;
  friend class Triangulation_dsul_edge_iterator_2<
                   Triangulation_data_structure_using_list_2<Vb,Fb> >;
  friend class Triangulation_dsul_vertex_iterator_2<
                   Triangulation_data_structure_using_list_2<Vb,Fb> >;

public:
  typedef Vb  Vertex_base;
  typedef Fb  Face_base;
  //typedef typename Vertex_base::Point                        Point;
  typedef Triangulation_dsul_vertex_2<Vertex_base,Face_base> Vertex;
  typedef Triangulation_dsul_face_2<Vertex_base,Face_base>   Face;
  typedef std::pair<Face*, int>                              Edge;

  typedef Triangulation_data_structure_using_list_2<Vertex_base,Face_base> Tds;
  typedef Triangulation_dsul_iterator_base_2<Tds>        Iterator_base;
  typedef Triangulation_dsul_face_iterator_2<Tds>        Face_iterator;
  typedef Triangulation_dsul_vertex_iterator_2<Tds>      Vertex_iterator;
  typedef Triangulation_dsul_edge_iterator_2<Tds>        Edge_iterator;

  typedef Triangulation_ds_face_circulator_2<Vertex,Face> 
							Face_circulator;
  typedef Triangulation_ds_vertex_circulator_2<Vertex,Face> 
							Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_2<Vertex,Face> 
							Edge_circulator;
  typedef std::list<Edge> List_edges;

protected:
  int _number_of_vertices; 
  int _dimension;
  Face  _dummy;

  // we maintain the list of cells to be able to traverse the triangulation
  // it starts with a "foo" element called _dummy
  // that will never be removed.
  // The list is circular, the foo element being used to recognize the end
  // of the list
 
public:
  //CREATORS - DESTRUCTORS
  Triangulation_data_structure_using_list_2(); 
  Triangulation_data_structure_using_list_2(const Tds &tds);
  ~Triangulation_data_structure_using_list_2();
  Tds& operator= (const Tds &tds);
  void swap(Tds &tds);
  void clear();

  //ACCESS FUNCTIONS
  int  dimension() const { return _dimension;  }
  int number_of_vertices() const {return _number_of_vertices;}
  int number_of_faces() const ;
  int number_of_edges() const;
  int number_of_full_dim_faces() const; //number of faces stored by tds
 
  
  // TEST FEATURES
  bool is_vertex(const Vertex* v) const;
  bool is_edge(const Vertex* va, const Vertex* vb) const;
  bool is_edge(const Vertex* va, const Vertex* vb, 
	       Face* &fr,  int &i) const;
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
  Vertex* insert_dim_up(Vertex *w = NULL, bool orient=true);

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
  Face* create_face(Face* f); //calls copy constructor of Face
  Face* create_face();

  void  delete_face(Face*);
  void  delete_vertex(Vertex*);

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  
  // HELPING
private:
  const Face& dummy() const {return _dummy;}
  Face& dummy() {return _dummy;}
  const Face* past_end() const {return &(_dummy);}
  Face* past_end() {return &(_dummy);}
  void  add_face(Face* newf);

public:
  void copy_tds(const Tds &tds);
  Vertex* copy_tds(const Tds &tds, const Vertex*);

  Vertex* file_input(std::istream& is, bool skip_first=false);
  void file_output(std::ostream& os,
		   Vertex* v = NULL,
		   bool skip_first=false) const;


  // SETTING (had to make them public for use in remove from Triangulations)
  void set_number_of_vertices(int n) {_number_of_vertices = n;}
  void set_dimension (int n) {_dimension = n ;}

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
    star_hole(v, 
	      edge_begin, 
	      edge_end, 
	      empty_list.begin(),
	      empty_list.end());
    return;    
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

template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb> ::
Triangulation_data_structure_using_list_2() 
  :  _number_of_vertices(0),_dimension(-2), _dummy()
{ }

template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb> ::
Triangulation_data_structure_using_list_2(const Tds &tds)
{
  copy_tds(tds);
}

template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb> ::
~Triangulation_data_structure_using_list_2()
{
  clear();
}

//assignement  
template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>&
Triangulation_data_structure_using_list_2<Vb,Fb> ::
operator= (const Tds &tds)
{
  copy_tds(tds);
  return *this;
}  

///ACCESS FUNCTIONS
template <class Vb, class Fb>
inline int 
Triangulation_data_structure_using_list_2<Vb,Fb> ::
number_of_faces() const 
{
  return ( dimension() < 2) ? 0 : (2*number_of_vertices()- 4);
}

template <class Vb, class Fb>
int
Triangulation_data_structure_using_list_2<Vb,Fb>::
number_of_edges() const
{
  switch (dimension()) {
  case 1:  return number_of_vertices();
  case 2:  return 3*number_of_vertices()-6;
  default: return 0;
  }
}
      
template <class Vb, class Fb>
int
Triangulation_data_structure_using_list_2<Vb,Fb>::
number_of_full_dim_faces() const
{
  switch (dimension()) {
  case -2: return 0;
  case -1: return 1;
  case 0:  return 2;
  case 1:  return number_of_edges();
  case 2:  return number_of_faces();
  default: return 0;
  }
}

template <class Vb, class Fb>
bool
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_vertex(const Vertex* v) const
{
  if (number_of_vertices() == 0) return false;
  for (Vertex_iterator vit = vertices_begin();
                       vit != vertices_end(); ++vit) {
   if ( v == &(*vit)) return true;
  }
  return false;
}

template <class Vb, class Fb>
bool
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_edge(const Vertex* va, const Vertex* vb) const
// returns true (false) if the line segment ab is (is not) an edge of t
//It is assumed that va is a vertex of t
{
  Vertex_circulator vc= va->incident_vertices(), done(vc);
  if ( vc == 0) return false;
  do {
    if( vb == &(*vc) ) {return true;} 
  } while (++vc != done);
  return false;
}
 

template <class Vb, class Fb>
bool
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_edge(const Vertex* va, const Vertex* vb, Face* &fr,  int & i) const
// assume va is a vertex of t
// returns true (false) if the line segment ab is (is not) an edge of t
// if true is returned (fr,i) is the edge ab
// with face fr on the right of a->b
{
  Face* fc = va->face(); 
  Face* start = fc;
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
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3) const
{
  Face* f;
  return is_face(v1,v2,v3,f);
}

template <class Vb, class Fb>
bool 
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_face(const Vertex* v1, const Vertex* v2, const Vertex* v3,
      Face* &f) const
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
Triangulation_data_structure_using_list_2<Vb,Fb>::
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
  
template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
insert_first( )
{
  CGAL_triangulation_precondition( number_of_vertices() == 0 &&
				   dimension()==-2 );
  return insert_dim_up();
}

template < class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex* 
Triangulation_data_structure_using_list_2<Vb,Fb>::
insert_second()
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 &&
				   dimension()==-1 );
  return insert_dim_up();

}


template <  class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
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


template <  class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
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


template <  class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
insert_dim_up(Vertex *w, bool orient)
{
  // the following function insert 
  // a vertex  v which is outside the affine  hull of Tds
  // The triangulation will be starred from  v and w 
  // ( geometrically w=  // the infinite vertex )
  // w=NULL for first and second insertions
  // orient governs the orientation of the resulting triangulation

  Vertex * v = create_vertex();
  set_dimension(dimension() + 1);
  //set_number_of_vertices(number_of_vertices() + 1);
  Face* f1;
  Face* f2;
    
  switch (dimension()) { //it is the resulting dimension
  case -1:
    f1 = create_face(v,NULL,NULL);
    v->set_face(f1);
    break;
  case 0 :
    f1 = &(*iterator_base_begin());
    f2 = create_face(v,NULL,NULL);
    f1->set_neighbor(0,f2);
    f2->set_neighbor(0,f1);
    v->set_face(f2);
    break;
  case 1 :
  case 2 :
    {
      std::list<Face *> faces_list;
      Iterator_base ib= Iterator_base(this); 
      Iterator_base ib_end = Iterator_base(this,1);
      for (; ib != ib_end ; ++ib){
	faces_list.push_back( & (*ib));
      }
      
      std::list<Face *>  to_delete;
      typename std::list<Face *>::iterator lfit = faces_list.begin();
      int i = dimension(); // maximun non NULL index in faces 
      Face *f, *g;

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
	f1= f->neighbor(i); i1= f1->index(f);
	f2= f->neighbor(j); i2 = f2->index(f);
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
Triangulation_data_structure_using_list_2<Vb,Fb>::
remove_degree_3(Vertex* v, Face* f)
   // remove a vertex of degree 3
  {
    CGAL_triangulation_precondition(v != NULL);
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
    delete_face(right);
    delete_face(left);
        
    delete_vertex(v);
    //set_number_of_vertices( number_of_vertices() -1);
  } 


  
template <class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
remove_dim_down(Vertex* v)
{
  Face* f;
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
    std::list<Face* > to_delete;
    std::list<Face* > to_downgrade;
    Iterator_base ib = iterator_base_begin();
    for( ; ib != iterator_base_end(); ++ib ){
      if ( ! ib->has_vertex(v) ) { to_delete.push_back(&(*ib));}
      else { to_downgrade.push_back(&(*ib));}
    }

    typename std::list<Face*>::iterator lfit = to_downgrade.begin();
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
Triangulation_data_structure_using_list_2<Vb,Fb>::  
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
  delete_face(g);
  delete_vertex(v);
  //set_number_of_vertices(number_of_vertices() -1);
  return;
}



template <class Vb, class Fb>
inline void
Triangulation_data_structure_using_list_2<Vb,Fb>::
remove_second(Vertex* v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 2 &&
 				  dimension() == 0);
  remove_dim_down(v);
  return;
}

    
template <class Vb, class Fb>
inline void
Triangulation_data_structure_using_list_2<Vb,Fb>::
remove_first(Vertex* v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 1 && 
 				  dimension() == -1);
  remove_dim_down(v);
  return; 
}

template <class Vb, class Fb>
inline
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
star_hole(List_edges& hole)
{
  Vertex* newv = create_vertex();
  star_hole(newv, hole);
  return newv;
}

template <class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
star_hole(Vertex* newv, List_edges& hole)
  // star the hole represented by hole around newv
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
{
   
  star_hole(newv, hole.begin(), hole.end());
  return;	    
}

template <class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
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


template <class Vb, class Fb>
inline
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_vertex()
{
  Vertex* newv= new Vertex();
  set_number_of_vertices(number_of_vertices()+ 1);
  return newv;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face(Face* f1, int i1, Face* f2, int i2, Face* f3, int i3)
{
  Face* newf = new Face(f1->vertex(cw(i1)),
			f2->vertex(cw(i2)),
			f3->vertex(cw(i3)),
			f2, f3, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  f3->set_neighbor(i3,newf);
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face(Face* f1, int i1, Face* f2, int i2)
{
  Face* newf = new Face(f1->vertex(cw(i1)),
			f2->vertex(cw(i2)),
			f2->vertex(ccw(i2)),
			f2, NULL, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face(Face* f1, int i1, Vertex* v)
{
  Face* newf = new Face(f1->vertex(cw(i1)), f1->vertex(ccw(i1)),v,
			NULL, NULL, f1);
  f1->set_neighbor(i1,newf);
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face()
{
  Face* newf= new Face();
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face( Face * f)
{
  Face* newf= new Face(*f); //calls copy constructor
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face(Vertex* v1, Vertex* v2, Vertex* v3)
{
  Face* newf= new Face(v1, v2, v3);
  add_face(newf);
  return newf;
}

template <class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Face*
Triangulation_data_structure_using_list_2<Vb,Fb>::
create_face(Vertex* v1, Vertex* v2, Vertex* v3,
	    Face* f1, Face* f2, Face* f3)
{
 Face* newf= new Face( v1, v2, v3, f1, f2, f3);
 add_face(newf);
 return(newf);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_using_list_2<Vb,Fb>::
delete_face(Face* f)
{
  f->previous()->set_next(f->next());
  f->next()->set_previous(f->previous());
  delete f;
}

template <class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
add_face(Face* newf)
{
  newf->set_next(dummy().next());
  newf->set_previous((dummy().next())->previous());
  dummy().set_next(newf);
  newf->next()->set_previous(newf);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_using_list_2<Vb,Fb>::
delete_vertex(Vertex* v)
{
  delete v;
  set_number_of_vertices(number_of_vertices() -1);
  return;
}

// CHECKING
template <  class Vb, class Fb>
bool
Triangulation_data_structure_using_list_2<Vb,Fb>::
is_valid(bool verbose, int level) const
{
  if(number_of_vertices() == 0){ 
    return (dimension() == -2);
  }

      
  bool result = (dimension()>= -1);
  CGAL_triangulation_assertion(result);

  //count and test the validity of the faces (for positive dimensions)
  Iterator_base ib = Iterator_base(this); 
  Iterator_base ib_end = Iterator_base(this,1);
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
    result = result &&  edge_count == 3*vertex_count - 6 ;
    CGAL_triangulation_assertion(result);
    result = result && ( face_count == 2*vertex_count - 4 );
    CGAL_triangulation_assertion( face_count == 2*vertex_count - 4 );
    break;
  default:
    result = false;
    CGAL_triangulation_assertion(result);
  }
  return result;
}


template <  class Vb, class Fb>
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
copy_tds(const Tds &tds, const Vertex* v) 
  // copy tds in this and return the vertex  corresponding to v
{
  CGAL_triangulation_precondition( tds.is_vertex(v));
  _number_of_vertices = tds.number_of_vertices();
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

  //vertex corresponding to v will be returned
  v_inf = static_cast<Vertex*>(V[v]);
    
  // create the faces
  //because faces are inserted at the head of the list
  //they have to be created in reverse order
  // to provide parallel faces and vertices iterators
//   for(Iterator_base ib = tds.iterator_base_begin();
//       ib != tds.iterator_base_end(); ++ib) {
  typedef  std::reverse_iterator<Iterator_base> RIB;
  for(RIB rib=RIB(tds.iterator_base_end());
      rib != RIB(tds.iterator_base_begin());
      ++rib) {
    f2 = create_face(&(*rib));
    F[&(*rib)]=  f2;
    f2->set_vertices((Vertex*) V[rib->vertex(0)],
		     (Vertex*) V[rib->vertex(1)],
		     (Vertex*) V[rib->vertex(2)] );
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

template <  class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
copy_tds(const Tds &tds)
{
  _number_of_vertices = tds.number_of_vertices();
  _dimension = tds.dimension();

  if(tds.number_of_vertices() == 0){return;}
  copy_tds(tds, &(*(tds.vertices_begin())));
}
 
template <  class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
swap(Tds &tds)
{
  int     nv = number_of_vertices();
  int    dim = dimension();
  Face*  tmp_next = dummy().next();
  Face*  tmp_prev = dummy().previous();

  //set dummy() next and previous 
  if (tds.number_of_full_dim_faces() == 0) {
    dummy().set_next(&dummy());
    dummy().set_previous(&dummy());
  }
  else {
    dummy().set_next(tds.dummy().next());
    dummy().set_previous(tds.dummy().previous());
    dummy().next()->set_previous(&dummy());
    dummy().previous()->set_next(&dummy());
  }

  //set tds.dummy() next and previous
  if(nv == 0 ){
    tds.dummy().set_next(& (tds.dummy()));
    tds.dummy().set_previous(& (tds.dummy()));
  }
  else {
    tds.dummy().set_next(tmp_next);
    tds.dummy().set_previous(tmp_prev);
    tds.dummy().next()->set_previous(& (tds.dummy()));
    tds.dummy().previous()->set_next(& (tds.dummy()));
  }

  set_number_of_vertices(tds.number_of_vertices());
  set_dimension(tds.dimension());

  tds.set_number_of_vertices(nv);
  tds.set_dimension(dim);
}

template <  class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
clear()
{
  if(number_of_vertices() == 0) return;
  
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
    delete_face(*lfit);
  }
  for( typename std::list<Vertex*>::iterator lvit=Vertices.begin();
       lvit != Vertices.end(); ++lvit) {
    delete_vertex(*lvit);
  }
  
  //set_number_of_vertices(0);
  set_dimension(-2);
  dummy().set_next(&dummy());
  dummy().set_previous(&dummy());
  return;
}

template < class Vb, class Fb>
void
Triangulation_data_structure_using_list_2<Vb,Fb>::
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
Triangulation_data_structure_using_list_2<Vb,Fb>::Vertex*
Triangulation_data_structure_using_list_2<Vb,Fb>::
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



template <  class Vb, class Fb>
std::istream&
operator>>(std::istream& is,  
	   Triangulation_data_structure_using_list_2<Vb,Fb>& tds) 
{
  tds.file_input(is);
  return is;
}


template <  class Vb, class Fb>
std::ostream&
operator<<(std::ostream& os, 
	   const Triangulation_data_structure_using_list_2<Vb,Fb>  &tds) 
{
   tds.file_output(os);
   return os;
}

CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_DATA_STRUCTURE_USING_LIST_2_H
