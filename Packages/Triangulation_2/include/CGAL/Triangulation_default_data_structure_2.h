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
// file     : Triangulatin/include/CGAL/Triangulation_default_data_structure_2.h
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

#include <pair.h>
#include <list.h>
#include <map.h>
#include <vector.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_ds_face_2.h>
#include <CGAL/Triangulation_ds_vertex_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

template <class Tds>
class CGAL_Triangulation_ds_iterator_base_2;

template <class Tds>
class CGAL_Triangulation_ds_face_iterator_2;

template <class Tds>
class CGAL_Triangulation_ds_vertex_iterator_2;

template <class Tds>
class CGAL_Triangulation_ds_edge_iterator_2;

template<class Vertex, class Face>
class CGAL_Triangulation_ds_face_circulator_2;

template<class Vertex, class Face>
class CGAL_Triangulation_ds_vertex_circulator_2;

template<class Vertex, class Face>
class CGAL_Triangulation_ds_edges_circulator_2;


template < class Gt , class Vb, class Fb>
class CGAL_Triangulation_default_data_structure_2
{
friend istream& operator>> CGAL_NULL_TMPL_ARGS
     (istream& is, CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds);
friend CGAL_Triangulation_ds_iterator_base_2<CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
friend CGAL_Triangulation_ds_face_iterator_2<CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
friend CGAL_Triangulation_ds_vertex_iterator_2<CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> >;
friend CGAL_Triangulation_ds_edge_iterator_2<CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> >;

public:
  typedef Gt Geom_traits;

  typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
  typedef pair<Face*, int>  Edge;

  typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
  typedef CGAL_Triangulation_ds_face_iterator_2<Tds> Face_iterator;
  typedef CGAL_Triangulation_ds_vertex_iterator_2<Tds> Vertex_iterator;
  typedef CGAL_Triangulation_ds_edge_iterator_2<Tds> Edge_iterator;


  typedef CGAL_Triangulation_ds_face_circulator_2<Vertex,Face> 
							Face_circulator;
  typedef CGAL_Triangulation_ds_vertex_circulator_2<Vertex,Face> 
							Vertex_circulator;
  typedef CGAL_Triangulation_ds_edge_circulator_2<Vertex,Face> 
							Edge_circulator;


  //creators
  CGAL_Triangulation_default_data_structure_2() 
    : _infinite_vertex(NULL),_number_of_vertices(0)
  { }

   CGAL_Triangulation_default_data_structure_2(const Geom_traits& gt) 
    : _geom_traits(gt), _infinite_vertex(NULL), 
      _number_of_vertices(0), _dimension(0)
      
  { }

  CGAL_Triangulation_default_data_structure_2(Vertex * v)
  {
    init(v);
    CGAL_triangulation_postcondition( is_valid() );
  }


  CGAL_Triangulation_default_data_structure_2(Vertex * v, 
					      const Geom_traits& gt)
    : _infinite_vertex(NULL), 
      _number_of_vertices(0), _geom_traits(gt)
  {
    init(v);
    CGAL_triangulation_postcondition( is_valid() );
  }

  CGAL_Triangulation_default_data_structure_2(const Tds &tds)
  {
    copy_tds(tds);
  }

   ~CGAL_Triangulation_default_data_structure_2()
  {
    clear();
  }
  
  //assignement
 Tds& operator= (const Tds &tds)
  {
    copy_tds(tds);
    return *this;
  }  
    
private:
  Geom_traits _geom_traits;
  Vertex* _infinite_vertex;
  int _number_of_vertices; 
  int _dimension;




public:
  // STATIC
  static int ccw(int i) {return (i+1) % 3;}
  static int cw(int i) {return (i+2) % 3;}
 
  //ACCESS FUNCTIONS
  int  dimension() const { return _dimension;  }
  int number_of_vertices() const {return _number_of_vertices;}
  int number_of_faces() const {
    return (number_of_vertices() <= 1) ? 0 : 2 * number_of_vertices() - 4;
  }
  const Geom_traits& geom_traits() const {return _geom_traits;}

public:
  //this should be private but
  // it is used by >> input operator
  // and apparently friend declaration does not work
  Vertex* infinite_vertex() const  {return _infinite_vertex;  }

private:
  Face* infinite_face() const
  {
    CGAL_triangulation_precondition( number_of_vertices() >= 2  &&
				     _infinite_vertex->face() != NULL );
    return _infinite_vertex->face();
  }

  // TEST IF INFINITE FEATURES
  bool is_infinite(const Face* f) const {
    return f->has_vertex(infinite_vertex());
  }

  bool is_infinite(const Vertex* v) const {
    return v == infinite_vertex();
  }

  bool is_infinite(const Face* f, int i) const {
    return ( is_infinite(f->vertex(ccw(i))) ||
	     is_infinite(f->vertex(cw(i))) );
  }

  bool is_infinite(Edge e) const {
      return is_infinite(e.first,e.second);
  }

  bool is_infinite(Edge_circulator& ec) const {
      return is_infinite(*ec);
  }

  bool is_infinite(Edge_iterator& ei) const {
    return is_infinite(*ei);
  }


  // SETTING
  // to be protected ?
public:
  void set_number_of_vertices(int n) {_number_of_vertices = n;}
  void set_dimension (int n) {_dimension = n ;}

private:
  void set_infinite_vertex(Vertex*  v) { _infinite_vertex = v;}

public:
  // MODIFY
   void flip(Face* f, int i)
    {
      CGAL_triangulation_precondition( dimension()==2);
      Face* n  = f->neighbor(i);
      int ni = n->index(f);
    
      Vertex*  v_cw = f->vertex(cw(i));
      Vertex*  v_ccw = f->vertex(ccw(i));

      CGAL_triangulation_assertion( f->vertex(cw(i)) == n->vertex(ccw(ni)));
      CGAL_triangulation_assertion( f->vertex(ccw(i)) == n->vertex(cw(ni)));	
      CGAL_triangulation_assertion( f->vertex(i) != n->vertex(ni));
      CGAL_triangulation_assertion( f == n->neighbor(ni) );
    
     
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
  


  void insert_second(Vertex* v)
  {
    CGAL_triangulation_precondition( number_of_vertices() == 1 &&
				     v != infinite_vertex());
    Face * f1 = new Face( _infinite_vertex, NULL, NULL);
    Face * f2 = new Face( v, NULL, NULL);
    f1->set_neighbor(0,f2);
    f2->set_neighbor(0,f1);
    _infinite_vertex->set_face(f1);
    v->set_face(f2);
    set_number_of_vertices(2);
    set_dimension(0);
    return;
  }



  void insert_in_face(Vertex* v, Face* f)
    //insert in face
    // vertex v will replace f->vertex(0)
  {
    CGAL_triangulation_precondition( v != NULL && f != NULL);
    
    Vertex* v0 = f->vertex(0);
    Vertex* v2 = f->vertex(2);
    Vertex* v1 = f->vertex(1);
    
    Face* n1 = f->neighbor(1);
    Face* n2 = f->neighbor(2);
    int i1,i2 ;
    if (n1 != NULL) {i1= cw(n1->index(f->vertex(cw(1))));}
    if (n2 != NULL) {i2 = cw(n2->index(f->vertex(cw(2))));}
    
    Face* f1 = new Face(v0, v, v2,
			f, n1, NULL);
    
    Face* f2 = new Face(v0, v1, v,
			f, NULL, n2);

    f1->set_neighbor(2, f2);
    f2->set_neighbor(1, f1);
    if (n1 != NULL) {n1->set_neighbor(i1,f1);}
    if (n2 != NULL) {n2->set_neighbor(i2,f2);}

    f->set_vertex(0, v);
    f->set_neighbor(1, f1);
    f->set_neighbor(2, f2);

    if( v0->face() == f  ) {  v0->set_face(f2); }
    v->set_face(f);

    set_number_of_vertices(number_of_vertices() +1);
  }

  void insert_in_edge(Vertex* v, Face* f, int i)
    //insert in the edge opposite to vertex i of face f
  {
    CGAL_triangulation_precondition(v != NULL && f != NULL); 
    CGAL_triangulation_precondition(dimension() >= 1);
    if (dimension() == 1) {CGAL_triangulation_precondition( i=3);}
    if (dimension() == 2) {CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);}
  
    if (dimension() == 1) {
      Face * g = new Face(v,f->vertex(1),NULL, NULL, NULL, NULL);
      g->set_neighbor(1,f); g->set_neighbor(0, f->neighbor(0));
      f->set_vertex(1,v); f->set_neighbor(0,g);
      set_number_of_vertices(number_of_vertices() +1);
    }

    else { //dimension() ==2
    Face* n = f->neighbor(i);
    int in = n->index(f);
    CGAL_triangulation_assertion( f->vertex(cw(i)) == n->vertex(ccw(in)) &&
				  f->vertex(ccw(i)) == n->vertex(cw(in)) );
    insert_in_face(v,f);
    flip(n,in); 
    }

    return;
  }




  void insert_outside_affine_hull(Vertex *v, Vertex *w, bool orient)
  {
  // the following function insert 
  // a vertex  v which is outside the affine  hull of a 1 dim or 0 dim triangulation
  // w is the infinite vertex of the triangulation
  // orient governs the orientation of the resulting triangulation

    list<Face *> faces_list;
    Face_iterator fit = faces_begin();
    for ( ; fit != faces_end() ; ++fit){
      faces_list.push_back( & (*fit));
    }

    list<Face *>  to_delete;
    list<Face *>::iterator lfit = faces_list.begin();
    int i = dimension()+1; // maximun non NULL index in faces after the insertion
    Face *f, *g;

    for ( ; lfit != faces_list.end() ; ++lfit) {
      f = * lfit;
      g = new Face( f);
      f->set_vertex(i,v); f->set_neighbor(i,g);
      g->set_vertex(i,w); g->set_neighbor(i,f);
      if (f->has_vertex(w)) to_delete.push_back(g); // flat face to be deleted later
    }

    lfit = faces_list.begin();
    for ( ; lfit != faces_list.end() ; ++lfit) {
      f = * lfit;
      g = f->neighbor(i);
      for(int j = 0; j < i ; j++) {
	g->set_neighbor(j, f->neighbor(j)->neighbor(i));
      }
    }

    // couldn't unify the code for reorientation mater
    if (dimension() == 0){
      lfit = faces_list.begin() ; ++lfit;
      f = *lfit;
      Vertex* vtemp = f->vertex(0); 
	f->set_vertex(0, f->vertex(1)) ; f->set_vertex(1,vtemp);
	Face* ftemp = f->neighbor(0); f->set_neighbor(0, f->neighbor(1)); 
	f->set_neighbor(1,ftemp);
    }
    else { // dimension == 1
      lfit = faces_list.begin();
      for( ;lfit  != faces_list.end(); ++lfit ){
	if (orient) {f = (*lfit)->neighbor(2);}
	else { f = *lfit;}
	Vertex* vtemp = f->vertex(0); 
	f->set_vertex(0, f->vertex(1)) ; f->set_vertex(1,vtemp);
	Face* ftemp = f->neighbor(0); f->set_neighbor(0, f->neighbor(1)); 
	f->set_neighbor(1,ftemp);
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
    set_number_of_vertices(number_of_vertices() + 1);
  }


  void remove_degree_3(Vertex* v, Face* f = NULL)
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


  void remove_down(Vertex* v)
  {

    CGAL_triangulation_precondition ( (dimension() == 1 && number_of_vertices() == 3) ||
				      (dimension() == 2 && number_of_vertices() > 3) );
    // the faces incident to v are down graded to one dimension
    // the other faces are deleted
    list<Face* > to_delete;
    list<Face* > to_downgrade;
    Face_iterator fit = faces_begin();
    if ( ! fit->has_vertex(v) ) { to_delete.push_back(&(*fit));}
    else { to_downgrade.push_back(&(*fit));}

    list<Face*>::iterator lfit = to_downgrade.begin();
    int j;
    Face * f;
    for( ; lfit !=  to_downgrade.end() ; lfit++) {
      f = *lfit; j = f->index(v);
      if (dimension() == 1){
	if (j == 0) {
	  f->set_vertex(0, f->vertex(1));
	  f->set_neighbor(0, f->neighbor(1));
	}
	f->set_vertex(1,NULL);
	f->set_neighbor(1,NULL);
      }
      else{ //dimension() == 2
	switch(j) {
	case 0 : 
	  f->set_vertex(0, f->vertex(1));
	  f->set_vertex(1, f->vertex(2));
	  f->set_neighbor(0, f->neighbor(1));
	  f->set_neighbor(2, f->neighbor(2));
	  break;
	case 1 :
	  f->set_vertex(1, f->vertex(0));
	  f->set_vertex(0, f->vertex(2));
	  f->set_neighbor(1, f->neighbor(0));
	  f->set_neighbor(0, f->neighbor(2));
	  break;
	case 2 :
	  break;
	}
	f->set_vertex(2,NULL);
	f->set_neighbor(2,NULL);
      }
      f->vertex(0)->set_face(f);
    }

    lfit = to_delete.begin();
    for( ; lfit !=  to_delete.end() ; lfit++) {
     delete *lfit;
    }

    
    delete v;
    set_number_of_vertices(number_of_vertices() -1);
    set_dimension(dimension() -1);
    return;
  }

  
  void remove_1D(Vertex* v)
  {
    CGAL_triangulation_precondition( dimension() == 1 &&
				     number_of_vertices() >= 3);
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



  void remove_second(Vertex* v)
  {
    CGAL_triangulation_precondition(number_of_vertices()== 2);
        
    delete v;
    set_number_of_vertices(1);
    infinite_vertex()->set_face(NULL);
    return;
  }
    

// ITERATOR METHODS
    Face_iterator faces_begin() const
    {
        Tds* ncthis = (Tds *)this;
        return Face_iterator(ncthis);
    }

    Face_iterator faces_end() const
    {
        Tds* ncthis = (Tds *)this;
        return Face_iterator(ncthis, 1);
    }

    Vertex_iterator vertices_begin() const
    {
        Tds* ncthis = (Tds*)this;
        return Vertex_iterator(ncthis);
    }

    Vertex_iterator vertices_end() const
    {
        Tds* ncthis = (Tds*)this;
        return Vertex_iterator(ncthis,1);
    }

    Edge_iterator edges_begin() const
    {
        Tds* ncthis = (Tds*)this;
        return Edge_iterator(ncthis);
    }

    Edge_iterator edges_end() const
    {
        Tds* ncthis = (Tds*)this;
        return Edge_iterator(ncthis,1);
    }





  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const
    {
        if (number_of_vertices() <= 3) return true;
    
        bool result = true;
	CGAL_triangulation_assertion( is_infinite(infinite_vertex()->face()));
	result = result && is_infinite(infinite_vertex()->face());			      
	// vertex count
        int vertex_count = 0;
        {
            Vertex_iterator it = vertices_begin();
    
            while(it != vertices_end()){
	        result = result && it->is_valid(verbose,level);
		CGAL_triangulation_assertion( it->is_valid(verbose, level) );
                ++vertex_count;
                ++it;
            }
        }
	result = result && (number_of_vertices() == vertex_count);
	CGAL_triangulation_assertion( number_of_vertices() == vertex_count );
    
        //edge count
	int edge_count = 0;
        {
            Edge_iterator it = edges_begin();
    
            while(it != edges_end()){
                ++edge_count;
                ++it;
            }
        }
    
        // face count
	int face_count = 0;
        {
            Face_iterator it = faces_begin();
    
            while(it != faces_end()) {
                ++face_count;
                result = result && it->is_valid(verbose, level);
                CGAL_triangulation_assertion( it->is_valid(verbose, level) );
		++it;
            }
        }
        
//         Face_circulator fc(infinite_vertex());
//         Face_circulator    fcdone(fc);
//         do {
// 	  result = result && fc->is_valid(verbose, level);
//           CGAL_triangulation_assertion( fc->is_valid(verbose, level) );
// 	  ++face_count;
// 	  ++edge_count;
//         }
//         while(++fc != fcdone);
    
	result = result && ( edge_count == 3*vertex_count - 6 );
        CGAL_triangulation_assertion( edge_count == 3*vertex_count -6 );
        result = result && ( face_count == 2*vertex_count - 4 );
        CGAL_triangulation_assertion( face_count == 2*vertex_count - 4 );
    
        return result;
    }

  //Helping functions
public:
   void init(Vertex*  v)
    {
        if( v == NULL ){
            return;
        } 
	
	set_infinite_vertex(v);

	if( v->face() == NULL ){ //only one vertex
	  set_number_of_vertices(1);
	  set_dimension(0);
	  return;
        } 
	
	if( v->face()->neighbor(1) == NULL ){ //two vertices
	  set_number_of_vertices(2);
	  set_dimension(0);
	  return;
        } 
	  
	if( v->face()->neighbor(2) == NULL ){  //One dimensional triangulation
	  set_dimension(1);
	}
	    
	else{
	  set_dimension(2);
	}
	  
	  //count number of vertices
	set_number_of_vertices(2);
	Vertex_iterator it = vertices_begin();
	while(it != vertices_end()){
	  ++_number_of_vertices;
	  ++it;
	}
	set_number_of_vertices(number_of_vertices()-2);
	return;
    }


public:
  void copy_tds(const Tds &tds)
  {
    map< void*, void*, less<void*> > V;
    map< void*, void*, less<void*> > F;
    Vertex*  v2;
    Face* f2;

    int n = tds.number_of_vertices();
    //non used?? int m = tds.number_of_faces();
    set_number_of_vertices(n);
    _geom_traits = tds.geom_traits();
    //the class Geom_traits is required to have a pertinent operator=

    // create the vertex at infinity
    v2 = new Vertex();
    V[tds.infinite_vertex()]=v2;
    set_infinite_vertex(v2);
            
    if(n == 0){  return ;    }

    if(n == 1) { 
       v2 = new Vertex(tds.finite_vertex()->point());
       V[tds.finite_vertex()]=v2;
       set_finite_vertex(v2);
       return;
    }

    // create the finite vertices
    {
      Vertex_iterator it=tds.vertices_begin();
      while (it != tds.vertices_end()) {
	V[&(*it)] = new Vertex( it->point() );
	++it;
      }
    }
  
     // create the finite faces
    {
      Face_iterator it = tds.faces_begin();
      while(it != tds.faces_end()){
	F[&(*it)]=  new Face( (Vertex*) V[it->vertex(0)],
			      (Vertex*) V[it->vertex(1)],
			      (Vertex*) V[it->vertex(2)] );
	++(it);
        }
    }

    //create the infinite faces
    {
      Face_circulator fc = tds.infinite_vertex()->incident_faces();
      Face_circulator done(fc);
      do{
	F[&(*fc)]=  new Face( (Vertex*) V[fc->vertex(0)],
			      (Vertex*) V[fc->vertex(1)],
			      (Vertex*) V[fc->vertex(2)] );
      }while(++fc != done);
    }

    // link the infinite vertex to a triangle
    infinite_vertex()->set_face( (Face*) F[tds.infinite_face()] );

    // link the finite vertices to a triangle
    {
      Vertex_iterator it = tds.vertices_begin();
      while(it != tds.vertices_end()) {
            v2 = (Vertex*) V[&(*it)];
            v2->set_face( (Face*) F[it->face()] );
            ++it;
        }
    }

    // hook neighbor pointers of the finite faces
    {
        Face_iterator  it = tds.faces_begin();
        while(it != tds.faces_end()){
          for(int j = 0; j < 3; j++){
            f2 = ((Face*) F[&(*it)]);
            f2->set_neighbor(j, (Face*) F[it->neighbor(j)] );
          }
          ++it;
        }
    }

    // hook neighbor pointers of the infinite faces
    {
        Face_circulator  fc = tds.infinite_vertex()->incident_faces(),
            done(fc);
        do{
          f2 = ((Face*) F[&(*fc)]);
          for(int j = 0; j < 3; j++){
            f2->set_neighbor(j, (Face*) F[fc->neighbor(j)] );
          }
        }while(++fc != done);
    }
    CGAL_triangulation_postcondition( is_valid() );
    return;

  }
 
  
  void swap(Tds &tds)
  {
    Geom_traits   t  = geom_traits();
    Vertex*  fv = finite_vertex();
    Vertex*  iv = infinite_vertex();
    int     nv = number_of_vertices();

    _geom_traits = tds.geom_traits(); 
    //the class Geom_traits is required to have a pertinent operator=
    set_finite_vertex(tds.finite_vertex());
    set_infinite_vertex(tds.infinite_vertex());
    set_number_of_vertices(tds.number_of_vertices());

    tds._geom_traits = t;
    tds.set_finite_vertex(fv);
    tds.set_infinite_vertex(iv);
    tds.set_number_of_vertices(nv);
  }

  void clear()
  {
   if(number_of_vertices() == 0) return;

   if(number_of_vertices()==1) {
     delete infinite_vertex();
     set_infinite_vertex(NULL);
     set_dimension(0);
     set_number_of_vertices(0);
     return;
   }
      
   list<Face*> Faces;
   list<Vertex*> Vertices;

    {
        Vertex_iterator it = vertices_begin(), done = vertices_end();
        do{
            Vertices.push_front(&(*it));
        }while(++it!=done);
    }
    {
        Face_iterator it = faces_begin(), done = faces_end();
        while(it!=done){
            Faces.push_front(&(*it));
            ++it;
        }
        // Face_circulator fc = infinite_vertex()->incident_faces(),
//             fcdone(fc);
//         do{
//             Faces.push_front(&(*fc));
//         }while(++fc != fcdone);
    }
    CGAL_triangulation_assertion( number_of_faces() == (int) Faces.size());
     
    {
        list<Face*>::iterator
          it=Faces.begin(),done=Faces.end();
        do{
            delete *it;
        }while(++it!=done);
    }
    {
        list<Vertex*>::iterator
          it=Vertices.begin(),done=Vertices.end();
        do{
            delete *it;
        }while(++it!=done);
    }

    set_infinite_vertex(NULL);
    set_number_of_vertices(0);
    set_dimension(0);
  }

};


template < class Gt , class Vb, class Fb>
istream&
operator>>(istream& is,  
	   CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb>& tds)
{
  typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb>::Vertex Vertex;
  typedef CGAL_Triangulation_ds_face_2<Vb,Fb>::Face Face;
  typedef Gt Geom_traits;
  typedef typename Vb::Point Point;

  if(tds.number_of_vertices() != 0){ //clear but keep infinite_vertex;
    Vertex* vtemp = tds.infinite_vertex();
    tds.clear();
    tds.set_infinite_vertex(vtemp);
    tds.infinite_vertex()->set_face(NULL);
    }

    int i = 0;
    int n, m, d;
    is >> n >> m >> d;

    tds.set_number_of_vertices(n);

    vector<Vertex* > V(n);
    vector<Face*> F(m);
    Vertex * v;
    Face * f;

    // the first point is at infinity

    Point p;
    is >> p;
    v = tds.infinite_vertex();
    v->set_point(p);
    V[i] = v;
    ++i;

    if(n == 1){
        return is;
    }
    // for the other points we create new vertices
    for(; i < n; i++) {
        Point p;
        is >> p;
        V[i] = new Vertex(p);
    }


    if(n == 2){
        return is;
    }
    // Creation of the faces
    for(i = 0; i < m; i++) {
        int i0, i1, i2;
        is >> i0 >> i1 >> i2;
        f = new Face(V[i0], V[i1], V[i2]);
        F[i] = f;
        // The face pointer of vertices is set too often,
        // but otherwise we had to use a further map
        V[i0]->set_face(f);
        V[i1]->set_face(f);
        V[i2]->set_face(f);
    }

    // Setting the neighbor pointers is the same for the
    // faces on the other side of the plane  and the other faces
    for(i = 0; i < m; i++) {
        int i0, i1, i2;
        is >> i0 >> i1 >> i2;
        f = F[i];
        f->set_neighbor(0, F[i0]);
        f->set_neighbor(1, F[i1]);
        f->set_neighbor(2, F[i2]);
    }
    tds.infinite_vertex()->set_face(f);
    return is;
}


template < class Gt, class Vb, class Fb>
ostream&
operator<<(ostream& os, 
	   const  CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb>  &tds)
{
  map< void*, int, less<void*> > V;
  map< void*, int, less<void*> > F;
  typedef CGAL_Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
  typedef  Tds::Vertex  Vertex;
  typedef  Tds::Face   Face;
  typedef  Tds::Edge   Edge;
  typedef  Tds::Face_iterator  Face_iterator;
  typedef  Tds::Vertex_iterator  Vertex_iterator;
  typedef  Tds::Edge_iterator  Edge_iterator;
  typedef  Tds::Face_circulator  Face_circulator;

  Vertex* v;

    int n = tds.number_of_vertices() + 1;
    int m = tds.number_of_faces();
    if(CGAL_is_ascii(os)){
        os << n << ' ' << m << ' ' << tds.dimension() << endl;
    } else {
        os << n << m << tds.dimension();
    }

    // write the vertex at infinity
    int i = 0;
    v = tds.infinite_vertex();
    V[v] = i;
    os << v->point();
    if(CGAL_is_ascii(os)){
        os << ' ';
    }
    if(n == 1){
        return os;
    }

    // write the finite vertices
    {
        Vertex_iterator
          it = tds.vertices_begin();

        while(it != tds.vertices_end()){
            V[&(*it)] = ++i;
            os << it->point();
            if(CGAL_is_ascii(os)){
                os << ' ';
            }
            ++it;
        }
    }
    CGAL_triangulation_assertion( (i+1) == n );
    if(CGAL_is_ascii(os)){ os << "\n";}

    if(n == 2){
        return os;
    }

    i = 0;
    // vertices of the finite faces
    {
        Face_iterator
          it = tds.faces_begin();

        while(it != tds.faces_end()){
            F[&(*it)] = i++;
            for(int j = 0; j < 3; j++){
                os << V[it->vertex(j)];
                if(CGAL_is_ascii(os)){
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
        Face_circulator
            fc = tds.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            F[&(*fc)] = i++;
            for(int j = 0; j < 3; j++){
                os << V[fc->vertex(j)];
                if(CGAL_is_ascii(os)){
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
        Face_iterator
            it = tds.faces_begin();
        while(it != tds.faces_end()){
            for(int j = 0; j < 3; j++){
                os << F[&(* it->neighbor(j))];
                if(CGAL_is_ascii(os)){
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
        Face_circulator
            fc = tds.infinite_vertex()->incident_faces(),
            done(fc);

        do{
            for(int j = 0; j < 3; j++){
                os << F[fc->neighbor(j)];
                if(CGAL_is_ascii(os)){
                    if(j==2) {
                        os << "\n";
                    } else {
                        os <<  ' ';
                    }
                }
            }
        }while(++fc != done);
    }

    
    return os;
}


#endif CGAL_TRIANGULATION_DEFAULT_DATA_STRUCTURE_2_H
