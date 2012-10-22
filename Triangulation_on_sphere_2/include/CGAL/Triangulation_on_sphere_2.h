#ifndef CGAL_TRIANGULATION_ON_SPHERE_2_H
#define CGAL_TRIANGULATION_ON_SPHERE_2_H

//#define HALF_SPHERE

#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_on_sphere_2.h>
//#include <CGAL/Triangulation_line_face_circulator_sphere_2.h>
#include <CGAL/Random.h>

#include <CGAL/spatial_sort.h>



//TO DO
//oriented_side
//insertion, location, removal dim 0,1
//March locate 2D improve code and cases Vertex and Edge of boundary
//Line face circulator & marchlocate LFC
//includes_edge
//constructions
//move->well_oriented, from_convex_hull, depends on how defining the convex hull
//compare walks
//?fill_hole_delaunay
//draw_triangulation, depends on constructions


namespace CGAL {

template < class Gt, class Tds > class Triangulation_on_sphere_2;
template < class Gt, class Tds > std::istream& operator>>
    (std::istream& is, Triangulation_on_sphere_2<Gt,Tds> &tr);
template < class Gt, class Tds >  std::ostream& operator<<
  (std::ostream& os, const Triangulation_on_sphere_2<Gt,Tds> &tr);
  

template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                             Triangulation_vertex_base_2<Gt>,
                             Triangulation_face_base_on_sphere_2<Gt> > >
class Triangulation_on_sphere_2
  : public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> <>
               (std::istream& is, Triangulation_on_sphere_2 &tr);
  typedef Triangulation_on_sphere_2<Gt,Tds>             Self;

public:
  typedef Tds                                     Triangulation_data_structure;
  typedef Gt                                      Geom_traits;
  typedef typename Geom_traits::Point_2           Point;
  // typedef typename Geom_traits::Segment_2      Segment;
  // typedef typename Geom_traits::Triangle_2     Triangle;
  typedef typename Geom_traits::Orientation_2     Orientation_2;
  typedef typename Geom_traits::Coradial_sphere_2 Coradial_sphere_2;
  typedef typename Geom_traits::Compare_x_2       Compare_x;
  typedef typename Geom_traits::Compare_y_2       Compare_y;
  typedef typename Geom_traits::FT                FT;
  typedef typename Tds::size_type                 size_type;
  typedef typename Tds::difference_type           difference_type;
 
  typedef typename Tds::Vertex                 Vertex;
  typedef typename Tds::Face                   Face;
  typedef typename Tds::Edge                   Edge;
  typedef typename Tds::Vertex_handle          Vertex_handle;
  typedef typename Tds::Face_handle            Face_handle;

  typedef typename Tds::Face_circulator        Face_circulator;
  typedef typename Tds::Vertex_circulator      Vertex_circulator;
  typedef typename Tds::Edge_circulator        Edge_circulator;

  typedef typename Tds::Face_iterator          Faces_iterator;
  typedef typename Tds::Edge_iterator          Edges_iterator;
  typedef typename Tds::Vertex_iterator        Vertices_iterator;

  enum Locate_type {VERTEX=0, 
		    EDGE, //1
		    FACE, //2
		    OUTSIDE_CONVEX_HULL, //3
		    OUTSIDE_AFFINE_HULL};//4 

protected:

  Gt  _gt;
  Tds _tds;
  Vertex_handle _pivot;
  bool _full_sphere;
  Face_handle _negative;

  mutable Random rng;  

public:

  // CONSTRUCTORS
  Triangulation_on_sphere_2(const Geom_traits& geom_traits=Geom_traits());
  Triangulation_on_sphere_2(const Point& sphere); 
  Triangulation_on_sphere_2(const Triangulation_on_sphere_2<Gt,Tds> &tr);        
  void clear();

  //Assignement
  Triangulation_on_sphere_2 &operator=(const Triangulation_on_sphere_2 &tr);

  //Helping
  void copy_triangulation(const Triangulation_on_sphere_2 &tr);
  void swap(Triangulation_on_sphere_2 &tr);
 
  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;  

  // TESTS
  bool is_edge(Vertex_handle va, Vertex_handle vb) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr,int & i) const;
  //bool includes_edge(Vertex_handle va, Vertex_handle vb,Vertex_handle& vbb, Face_handle& fr, int & i) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle &fr) const;
 
   //ACCESS FUNCTION
  const Geom_traits& geom_traits() const { return _gt;}
  const Tds & tds() const  { return _tds;}
  Tds & tds()   { return _tds;}
  const Vertex_handle& pivot() const { return _pivot; }
   bool& full_sphere()  { return _full_sphere; }
 
  int dimension() const { return _tds.dimension();}
  size_type number_of_vertices() const {return _tds.number_of_vertices();}
  size_type number_of_faces() const;         
  int number_of_negative_faces();
 

  
 //--------------------------------------------------------------TRAVERSING : ITERATORS AND CIRCULATORS------------------------------------------- 
  Faces_iterator faces_begin() const;
  Faces_iterator faces_end() const;
  Vertices_iterator vertices_begin() const;
  Vertices_iterator vertices_end() const;
  Edges_iterator edges_begin() const;
  Edges_iterator edges_end() const; 
  Face_circulator incident_faces( Vertex_handle v, Face_handle f = Face_handle()) const;
  Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f = Face_handle()) const;
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f = Face_handle()) const;

 
  size_type degree(Vertex_handle v) const;
  
  Vertex_handle mirror_vertex(Face_handle f, int i) const;
  int mirror_index(Face_handle v, int i) const;
  /*
  Line_face_circulator    line_walk(const Point& p,
				    const Point& q,
				    Face_handle f = Face_handle()) const;
  */

  //-----------------------------------------------------------------LOCATION-------------------------------------------------
  Face_handle march_locate_2D(Face_handle c, const Point& t,Locate_type& lt, int& li) const; 
  Face_handle march_locate_1D(const Point& t, Locate_type& lt, int& li) const ;
  Face_handle locate(const Point& p, Locate_type& lt, int& li, Face_handle start) const;
  Face_handle locate(const Point &p, Face_handle start) const;

  //------------------------------------------------------------------------PREDICATES----------------------------------------
  Orientation orientation(const Point& p, const Point& q, const Point& r) const;
  Orientation orientation_1(const Point& p, const Point& q) const;
  Orientation orientation(const Face_handle f) const;
  //Comparison_result compare_x(const Point& p, const Point& q) const;
  //Comparison_result compare_y(const Point& p, const Point& q) const;
	Oriented_side power_test(const Point &p,
							 const Point &q,
							 const Point &r,
							 const Point &s) const;
	Oriented_side power_test(const Point &p,
							 const Point &q,
							 const Point &r) const;
	Oriented_side power_test(const Point &p,
							 const Point &r) const;
	Oriented_side power_test(const Face_handle &f, 
							 const Point &p) const;
	Oriented_side power_test(const Face_handle& f, int i,
							 const Point &p) const;
	
  Oriented_side
  oriented_side(const Point &p0, const Point &p1,
	       const Point &p2, const Point &p) const;
    
  Bounded_side
  bounded_side(const Point &p0, const Point &p1,
	      const Point &p2, const Point &p) const;
    
  Oriented_side
  oriented_side(Face_handle f, const Point &p) const;

  bool xy_equal(const Point& p, const Point& q) const;

  bool 
  collinear_between(const Point& p, const Point& q, const Point& r) const;
 

  //------------------------------------------------------------------DEBUG---------------------------------------------------
  void show_all() const;
  void show_vertex(Vertex_handle vh) const;
  void show_face(Face_handle fh) const;

 
  //----------------------------------------------------------PUBLIC REMOVE---------------------------------------------------
 void make_hole(Vertex_handle v, std::list<Edge> & hole);
  
  Face_handle create_face(Face_handle f1, int i1,
			  Face_handle f2, int i2,
			  Face_handle f3, int i3);
  
  Face_handle create_face(Face_handle f1, int i1,
			  Face_handle f2, int i2);
  Face_handle create_face();
  
  Face_handle create_face(Face_handle f, int i, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3,
			  Face_handle f1, Face_handle f2, Face_handle f3);

  Face_handle create_face(Face_handle);
  void delete_face(Face_handle f);
  void delete_vertex(Vertex_handle v);
 

  //           IN/OUT
  Vertex_handle file_input(std::istream& is);
  void file_output(std::ostream& os) const;

   /*---------------------------------------------------------------------TEMPLATE MEMBERS--------------------------------------*/
 public: 
 
  template<class EdgeIt, class FaceIt>
  Vertex_handle star_hole( const Point& p, 
			   EdgeIt edge_begin,
			   EdgeIt edge_end,
			   FaceIt face_begin,
			   FaceIt face_end) {
    typedef typename Triangulation_data_structure::Edge  Tds_Edge;
    typedef typename Triangulation_data_structure::Face  Tds_Face;
    Vertex_handle v = _tds.star_hole( edge_begin, edge_end,
				      face_begin, face_end);
    v->set_point(p);
    return v;
  }
 
  template<class FaceIt>
  void delete_faces(FaceIt face_begin, FaceIt face_end)
  {
    FaceIt fit=face_begin;
    for(;fit!=face_end;++fit)
    {
      delete_face(*fit);
    }
  }

};



// CONSTRUCTORS
template <class Gt, class Tds >
Triangulation_on_sphere_2<Gt, Tds>::
  Triangulation_on_sphere_2(const Geom_traits& geom_traits) 
  : _gt(geom_traits), _tds(), _full_sphere(false)                           
{}

template <class Gt, class Tds >
Triangulation_on_sphere_2<Gt, Tds>::
Triangulation_on_sphere_2(const Point& sphere) 
  : _gt(sphere), _tds(), _full_sphere(false)                           
{}

// copy constructor duplicates vertices and faces
template <class Gt, class Tds >
Triangulation_on_sphere_2<Gt, Tds>::
Triangulation_on_sphere_2(const Triangulation_on_sphere_2 &tr) 
  : _gt(tr._gt), _tds(tr._tds)
{}

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>:: 
clear()
{
  _tds.clear(); 
  _full_sphere=false;
}
  // Helping functions
  
template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::   
copy_triangulation(const Triangulation_on_sphere_2 &tr)
{
  _tds.clear();
  _gt = tr._gt;
  _tds = tr._tds;
}

  //Assignement
template <class Gt, class Tds >
Triangulation_on_sphere_2<Gt, Tds> &
Triangulation_on_sphere_2<Gt, Tds>::
operator=(const Triangulation_on_sphere_2 &tr)
{
  copy_triangulation(tr);
  return *this;
}

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>:: 
swap(Triangulation_on_sphere_2 &tr)
{
  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t; 
}

//ACCESS FUNCTION

//CHECKING
template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::size_type
Triangulation_on_sphere_2<Gt, Tds>::
number_of_faces() const
{
  size_type count = _tds.number_of_faces();
  return count;
}

template <class Gt, class Tds >
bool
Triangulation_on_sphere_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  if (dimension() <= 0 ||
     (dimension()==1 && number_of_vertices() == 3 ) ) return result;
 
  if (dimension() == 1) {
    Vertices_iterator vit=vertices_begin();
    }    

  else { //dimension() == 2

    for(Faces_iterator it=faces_begin(); 
	it!=faces_end(); it++) {
      Orientation s = orientation(it->vertex(0)->point(),
				  it->vertex(1)->point(),
				  it->vertex(2)->point());
      CGAL_triangulation_assertion( s == LEFT_TURN || it->is_negative() );
      result = result && ( s == LEFT_TURN || it->is_negative() );
    }

    // check number of faces. This cannot be done by the Tds
    // which does not know the number of components nor the genus
    result = result && (number_of_faces() == (2*number_of_vertices()- 4) );
                                      
    CGAL_triangulation_assertion( result);
  }
  return result;
}

//TESTS

template <class Gt, class Tds >
inline bool
Triangulation_on_sphere_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb) const
{
  return _tds.is_edge( va, vb);
}

template <class Gt, class Tds >
inline bool
Triangulation_on_sphere_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i) const
{
  return _tds.is_edge(va, vb, fr, i);
}
/*
template <class Gt, class Tds >
bool 
Triangulation_on_sphere_2<Gt, Tds>::
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
  Edge_circulator ec = incident_edges(va), done(ec);
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
}*/

template <class Gt, class Tds >
inline bool 
Triangulation_on_sphere_2<Gt, Tds>::
is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
{
  return _tds.is_face(v1, v2, v3);
}

template <class Gt, class Tds >
inline bool 
Triangulation_on_sphere_2<Gt, Tds>::
is_face(Vertex_handle v1, 
	Vertex_handle v2, 
	Vertex_handle v3,
	Face_handle &fr) const
{
  return _tds.is_face(v1, v2, v3, fr);
}

//---------------------------------------------------------------------------/POINT LOCATION---------------------------------------//

template <class Gt, class Tds >    
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
march_locate_1D(const Point& t,
		Locate_type& lt,
		int& li) const
{
  Face_handle f = edges_begin()->first;

  //first check if the circle is coplanar with o
  Orientation pqr = orientation(f->vertex(0)->point(), 
				f->vertex(1)->point(),
				f->neighbor(0)->vertex(1)->point());
  if( pqr != ON_ORIENTED_BOUNDARY ){
    if(xy_equal(t,f->vertex(0)->point())){
      lt = VERTEX;
      li = 0;
      return f;
    }
    if(xy_equal(t,f->vertex(1)->point())){
      lt = VERTEX;
      li = 1;
      return f;
    }
    if(xy_equal(t,f->neighbor(0)->vertex(1)->point())){
      lt = VERTEX;
      li = 1;
      return f->neighbor(0);
    }
    lt = OUTSIDE_AFFINE_HULL;
    li = 4 ;// should not be used
    return f;
  }

  //if not ckeck wether t is in the plane of the circle or not
  Orientation pqt = orientation(f->vertex(0)->point(), 
				f->vertex(1)->point(),
				t);
  if(pqt != ON_ORIENTED_BOUNDARY) {
    lt = OUTSIDE_AFFINE_HULL;
    li = 4 ;// should not be used
    return f;
  }

  Edges_iterator eit=edges_begin();
  
  //find negative face if exists
  for(; eit!=edges_end(); ++eit){
    if(eit->first->is_negative())
    {f=eit->first; break;}
  }

  //if so check wether t is on the convex_hull or not
  if(f->is_negative()){
    //show_face(f);
    if(xy_equal(t,f->vertex(0)->point())){
      lt = VERTEX;
      li = 0;
      return f;
    }
    if(xy_equal(t,f->vertex(1)->point())){
      lt = VERTEX;
      li = 1;
      return f;
    }
    if (!collinear_between(f->vertex(0)->point(),f->vertex(1)->point(),t)) {
      lt = OUTSIDE_CONVEX_HULL;
      li = 4;
      return f;
    }
  }
  
  //then t is in the convex hull	
  eit=edges_begin();
  Vertex_handle u,v;
  for( ; eit!=edges_end() ; ++eit) {
    if(!eit->first->is_negative()){
      u = (*eit).first->vertex(0);
      v = (*eit).first->vertex(1);
      if(xy_equal(t,v->point())){
	lt = VERTEX;
	li = 1;
	return (*eit).first;
      }
      if(collinear_between(u->point(), v->point(), t)){
	lt = EDGE;
	li =  2;
	return (*eit).first;
      }
    }
  }
  //Provisoire pas satisfaisant du tout!
  lt=OUTSIDE_CONVEX_HULL;
  li=4;
  return f;
  CGAL_triangulation_assertion(false);
  return Face_handle();

}
  
template <class Gt, class Tds >   typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
march_locate_2D(Face_handle c,
		const Point& t,
		Locate_type& lt,
	        int& li) const
{
  CGAL_triangulation_assertion(!c->is_negative());

	
		
  Face_handle prev = Face_handle();
  bool first=true;

  while (1) {
	  
	  const Point & p0 = c->vertex( 0 )->point();
	  const Point & p1 = c->vertex( 1 )->point();
	  const Point & p2 = c->vertex( 2 )->point();
	  
	   
	  if( power_test(c,t)!=ON_NEGATIVE_SIDE){
		  if (c ->is_negative()){
			  lt = OUTSIDE_CONVEX_HULL;
			  li = 4;
		      return c;
		  }
		
		  else {
			  lt = FACE;
			  li = 4;		  
			  return c;
		  }
	  }

		

  
    
    // Instead of testing c's edges in a random order we do the following
    // until we find a neighbor to go further:
    // As we come from prev we do not have to check the edge leading to prev
    // Now we flip a coin in order to decide if we start checking the
    // edge before or the edge after the edge leading to prev
    // We do loop unrolling in order to find out if this is faster.
    // In the very beginning we do not have a prev, but for the first step 
    // we do not need randomness
    int left_first =rng.template get_bits<1>();
    Orientation o0, o1, o2;
  
    /************************FIRST*************************/
    if(first){
      prev = c;
      first = false;
      o0 = orientation(p0,p1,t);//classic march locate
      if ( o0 == NEGATIVE ) {
	c = c->neighbor( 2 );
	continue;
      }
      o1 = orientation(p1,p2,t);
      if ( o1 == NEGATIVE ) {
	c = c->neighbor( 0 );
	continue;
      }
      o2 = orientation(p2,p0,t);
      if ( o2 == NEGATIVE ) {
	c = c->neighbor( 1 );
	continue;
      }
    } /***********************END FIRST*****************/
  
    //****LEFT****//

    else if(left_first){
      if(c->neighbor(0) == prev){
	prev = c;
	o0 = orientation(p0,p1,t);
	if ( o0 == NEGATIVE ) {
	  c = c->neighbor( 2 );
	  continue;
	}
	o2 = orientation(p2,p0,t);
	if ( o2 == NEGATIVE ) {
	  c = c->neighbor( 1 );
	  continue;
	}
	o1 = orientation(p1,p2,t);
      } else if(c->neighbor(1) == prev){
	prev = c;
	o1 = orientation(p1,p2,t);
	if ( o1 == NEGATIVE ) {
	  c = c->neighbor( 0 );
	  continue;
	}
	o0 = orientation(p0,p1,t);
	if ( o0 == NEGATIVE ) {
	  c = c->neighbor( 2 );
	  continue;
	}
	o2 = orientation(p2,p0,t);
      } else {
	prev = c;
	o2 = orientation(p2,p0,t);
	if ( o2 == NEGATIVE ) {
	  c = c->neighbor( 1 );
	  continue;
	}
	o1 = orientation(p1,p2,t);
	if ( o1 == NEGATIVE ) {
	  c = c->neighbor( 0 );
	  continue;
	}
	o0 = orientation(p0,p1,t);
      }
    }
   
    else {
      if(c->neighbor(0) == prev){
	prev = c;
	o2 = orientation(p2,p0,t);
	if ( o2 == NEGATIVE ) {
	  c = c->neighbor( 1 );
	  continue;
	}
	o0 = orientation(p0,p1,t);
	if ( o0 == NEGATIVE ) {
	  c = c->neighbor( 2 );
	  continue;
	}
	o1 = orientation(p1,p2,t);
      } else if(c->neighbor(1) == prev){
	prev = c;
	o0 = orientation(p0,p1,t);
	if ( o0 == NEGATIVE ) {
	  c = c->neighbor( 2 );
	  continue;
	}
	o1 = orientation(p1,p2,t);
	if ( o1 == NEGATIVE ) {
	  c = c->neighbor( 0 );
	  continue;
	}
	o2 = orientation(p2,p0,t);
      } else {
	prev = c;
	o1 = orientation(p1,p2,t);
	if ( o1 == NEGATIVE ) {
	  c = c->neighbor( 0 );
	  continue;
	}
	o2 = orientation(p2,p0,t);
	if ( o2 == NEGATIVE ) {
	  c = c->neighbor( 1 );
	  continue;
	}
	o0 = orientation(p0,p1,t);
      }
    }
  
 
    //********FACE LOCATED************/

    int sum = ( o0 == COLLINEAR )
      + ( o1 == COLLINEAR )
      + ( o2 == COLLINEAR );
    
    switch (sum) {
    case 0:
    case -1:
    case -2:
      {
	lt = FACE;
	li = 4;
	break;
      }
    case 1:
      { 
	lt = EDGE;
	li = ( o0 == COLLINEAR ) ? 2 :
	  ( o1 == COLLINEAR ) ? 0 :
	  1;
	break;
      }
    case 2:
      {
	lt = VERTEX;
	li = ( o0 != COLLINEAR ) ? 2 :
	  ( o1 != COLLINEAR ) ? 0 :
	  1;
	break;
      }
    case 3:
      {
	lt=FACE;
	li=4;
      }
    }
  
   return c;
    
  }
} 

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt,Tds>::
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle start) const
{
  if (dimension() == -2) {         //empty triangulation
      lt = OUTSIDE_AFFINE_HULL;
      li = 4; // li should not be used in this case
      return Face_handle();
  }
  
  if(dimension() == -1){      //1 vertex
    Point q=vertices_begin()->point();
    if(xy_equal(q,p)){
      lt=VERTEX;
    }
    else{
      lt=OUTSIDE_AFFINE_HULL;
    }
    li=4;
    return Face_handle();
  }
 
  if( dimension() == 0) {
    Vertex_handle v=vertices_begin();
    Face_handle f=v->face();
    if(xy_equal(p,v->point())){
     lt=VERTEX;      
     li=0;  
     return f;
    }else if(xy_equal(p,f->neighbor(0)->vertex(0)->point())){
      lt=VERTEX;      
      li=0;  
      return f->neighbor(0);
     }
    else{
      lt=OUTSIDE_AFFINE_HULL;
      li=4;  
    }
    return Face_handle();
  }
  
  if(dimension() == 1){
    return march_locate_1D(p, lt, li);
  }
   
  if(start==Face_handle()){
      start=faces_begin();
  }
 
  if(start->is_negative()){
	  
   /* if(!start->neighbor(0)->is_negative())
      start=start->neighbor(0);
    else if(!start->neighbor(1)->is_negative())
      start=start->neighbor(1);
    else {
		bool test = start->is_negative();
	
      start=start->neighbor(2);
	}*/
	  for (Faces_iterator it = this->_tds.face_iterator_base_begin(); it != faces_end(); it++) {  
	  
		  if(!it->is_negative()){
		  start = it;
			  break;
		  }
	  
	  }

  }

#if ( ! defined(CGAL_ZIG_ZAG_WALK)) && ( ! defined(CGAL_LFC_WALK))
#define CGAL_ZIG_ZAG_WALK
#endif

#ifdef CGAL_ZIG_ZAG_WALK
  Face_handle res1;
  res1 = march_locate_2D(start, p, lt, li);  
 
    
 
#endif

#ifdef CGAL_LFC_WALK
  Locate_type lt2;
  int li2;
  Face_handle res2 = march_locate_2D_LFC(start, p, lt2, li2);
#endif

#if defined(CGAL_ZIG_ZAG_WALK) && defined(CGAL_LFC_WALK)
  compare_walks(p,
		res1, res2,
		lt, lt2,
		li, li2);
#endif

#ifdef CGAL_ZIG_ZAG_WALK
  return res1;
#endif

#ifdef CGAL_LFC_WALK
  lt = lt2;
  li = li2;
  return res2;
#endif

}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>:: Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
locate(const Point &p,
       Face_handle start) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, start);
}




 //--------------------------------------------------------------------ITERATORS AND CIRCULATORS--------------------------------------
template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Faces_iterator
Triangulation_on_sphere_2<Gt, Tds>::
faces_begin() const
{
  return _tds.faces_begin();
} 


template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Faces_iterator
Triangulation_on_sphere_2<Gt, Tds>::
faces_end() const
{
  return _tds.faces_end();;
}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Vertices_iterator
Triangulation_on_sphere_2<Gt, Tds>::
vertices_begin() const
{
  return _tds.vertices_begin();
}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Vertices_iterator
Triangulation_on_sphere_2<Gt, Tds>::
vertices_end() const
{
  return _tds.vertices_end();
}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Edges_iterator
Triangulation_on_sphere_2<Gt, Tds>::
edges_begin() const
{
  return _tds.edges_begin();
}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Edges_iterator
Triangulation_on_sphere_2<Gt, Tds>::
edges_end() const
{
  return _tds.edges_end();
}

template <class Gt, class Tds >
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_circulator
Triangulation_on_sphere_2<Gt, Tds>::
incident_faces(Vertex_handle v, Face_handle f) const
{
  return _tds.incident_faces(v,f);
} 


template <class Gt, class Tds >
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Vertex_circulator
Triangulation_on_sphere_2<Gt, Tds>::  
incident_vertices(Vertex_handle v, Face_handle f) const
{
  return _tds.incident_vertices(v,f);
}

template <class Gt, class Tds >
inline
typename Triangulation_on_sphere_2<Gt, Tds>::size_type
Triangulation_on_sphere_2<Gt, Tds>::    
degree(Vertex_handle v) const
{
  return _tds.degree(v);
}

template <class Gt, class Tds >
inline
int
Triangulation_on_sphere_2<Gt, Tds>::    
mirror_index(Face_handle f, int i) const
{
  return _tds.mirror_index(f,i);
}

template <class Gt, class Tds >
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Vertex_handle
Triangulation_on_sphere_2<Gt, Tds>::    
mirror_vertex(Face_handle f, int i) const
{
  return _tds.mirror_vertex(f,i);
}


//------------------------------------------------------------------------------PREDICATES-----------------------------------------------------------------
template <class Gt, class Tds >
inline
Orientation
Triangulation_on_sphere_2<Gt, Tds>::
orientation(const Point& p, const Point& q,const Point& r ) const
{
  return geom_traits().orientation_2_object()(p,q,r);
}

template <class Gt, class Tds >
inline
Orientation
Triangulation_on_sphere_2<Gt, Tds>::
orientation_1(const Point& p, const Point& q ) const
{
  return geom_traits().orientation_1_object()(p,q);
}

template <class Gt, class Tds >
Orientation
Triangulation_on_sphere_2<Gt, Tds>::
orientation(const Face_handle f) const
{
  return  orientation(f->vertex(0)->point(),f->vertex(1)->point(),f->vertex(2)->point());
}
/*
template <class Gt, class Tds >
inline
Comparison_result
Triangulation_on_sphere_2<Gt, Tds>::
compare_x(const Point& p, const Point& q) const
{
  return geom_traits().compare_x_2_object()(p,q);
}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_on_sphere_2<Gt, Tds>::
compare_y(const Point& p, const Point& q) const
{
  return geom_traits().compare_y_2_object()(p,q);
}
*/
	
template <class Gt, class Tds >
Oriented_side
Triangulation_on_sphere_2<Gt, Tds>::
oriented_side(const Point &p0, const Point &p1,
	      const Point &p2, const Point &p) const
{
  Orientation o1 = orientation(p0, p1, p),
              o2 = orientation(p1, p2, p),
              o3 = orientation(p2, p0, p);

  if(orientation(p0, p1, p2)==POSITIVE){
    if(o1==POSITIVE)
      if(o2==POSITIVE)
	if(o3==POSITIVE)
	  return ON_POSITIVE_SIDE;

 if (o1 == COLLINEAR){
      if (o2 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
      return ON_NEGATIVE_SIDE;
    }
    if (o2 == COLLINEAR){
      if (o1 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
      return ON_NEGATIVE_SIDE;
    }
    if (o3 == COLLINEAR){
      if (o2 == COLLINEAR ||  o1 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
      return ON_NEGATIVE_SIDE;
    }
	  
	  if(o1 == COLLINEAR || o2==COLLINEAR || o3==COLLINEAR)
		  return ON_ORIENTED_BOUNDARY;
    return ON_NEGATIVE_SIDE;
	  
  }else{
    if(o1==POSITIVE && o2==POSITIVE && o3==POSITIVE)
	  return ON_POSITIVE_SIDE;

    if(o1==NEGATIVE && o2==NEGATIVE && o3==NEGATIVE )
      return ON_NEGATIVE_SIDE;
    else
      return ON_NEGATIVE_SIDE;
  }
}
	
	
	
	template < class Gt, class Tds >
	Oriented_side
	Triangulation_on_sphere_2<Gt,Tds>::
	power_test(const Face_handle &f, const Point &p) const
	{
		return power_test(f->vertex(0)->point(),
						  f->vertex(1)->point(),
						  f->vertex(2)->point(),p);
	}
	
	template < class Gt, class Tds >
	Oriented_side
	Triangulation_on_sphere_2<Gt,Tds>::
	power_test(const Face_handle& f, int i,
			   const Point &p) const
	{
		CGAL_triangulation_precondition (
										 orientation(f->vertex(ccw(i))->point(),
													 f->vertex( cw(i))->point(),
													 p)
										 == COLLINEAR);
		return  power_test(f->vertex(ccw(i))->point(),
						   f->vertex( cw(i))->point(),
						   p);
	}
	
	template < class Gt, class Tds >
	inline
	Oriented_side
	Triangulation_on_sphere_2<Gt,Tds>::
	power_test(const Point &p,
			   const Point &q,
			   const Point &r,
			   const Point &s) const
	{
		return geom_traits().power_test_2_object()(p,q,r,s);
	}
	
	template < class Gt, class Tds >
	inline
	Oriented_side
	Triangulation_on_sphere_2<Gt,Tds>::
	power_test(const Point &p,
			   const Point &q,
			   const Point &r) const
	{
		if(number_of_vertices()==2)
			if(orientation_1(p,q)==COLLINEAR)
				return ON_POSITIVE_SIDE;
		return geom_traits().power_test_2_object()(p,q,r);
	}
	
	
	
	
	

	/*template <class Gt, class Tds >
	Oriented_side
	Triangulation_on_sphere_2<Gt, Tds>::
	oriented_side(const Point &p0, const Point &p1,
				  const Point &p2, const Point &p) const
	{
		Orientation o1 = orientation(p0, p1, p),
		o2 = orientation(p1, p2, p),
		o3 = orientation(p2, p0, p);
		
		//if(orientation(p0, p1, p2)==POSITIVE &&orientation(p0, p1, p2, p)==POSITIVE ){
			//return ON_POSITIVE_SIDE;
		//	//}
		if(orientation(p0, p1, p2)==NEGATIVE &&geom_traits().power_test_2_object()(p0, p1, p2, p)== NEGATIVE ){
			return ON_POSITIVE_SIDE;
		}
		
			
			if (o1 == COLLINEAR){
				if (o2 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			if (o2 == COLLINEAR){
				if (o1 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			if (o3 == COLLINEAR){
				if (o2 == COLLINEAR ||  o1 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			
			return ON_NEGATIVE_SIDE;
		
				//return ON_NEGATIVE_SIDE;
		
	}	*/
	
	/*template <class Gt, class Tds >
	Oriented_side
	Triangulation_on_sphere_2<Gt, Tds>::
	oriented_side(const Point &p0, const Point &p1,
				  const Point &p2, const Point &p) const
	{
		Orientation o1 = orientation(p0, p1, p),
		o2 = orientation(p1, p2, p),
		o3 = orientation(p2, p0, p);
		
		if(orientation(p0, p1, p2)==POSITIVE){
			if(o1==POSITIVE)
				if(o2==POSITIVE)
					if(o3==POSITIVE)
						return ON_POSITIVE_SIDE;
			
			if (o1 == COLLINEAR){
				if (o2 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			if (o2 == COLLINEAR){
				if (o1 == COLLINEAR ||  o3 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			if (o3 == COLLINEAR){
				if (o2 == COLLINEAR ||  o1 == COLLINEAR) return ON_ORIENTED_BOUNDARY;
				return ON_NEGATIVE_SIDE;
			}
			
			return ON_NEGATIVE_SIDE;
		}else{
			if(o1==POSITIVE)
				if(o2==POSITIVE)
					if(o3==POSITIVE)
						return ON_POSITIVE_SIDE;
			
			if(o1==NEGATIVE && o2==NEGATIVE && o3==NEGATIVE )
				return ON_NEGATIVE_SIDE;
			else
				return ON_POSITIVE_SIDE;
		}
	}
	
	
	
	
	*/
	
	

template <class Gt, class Tds >
Bounded_side
Triangulation_on_sphere_2<Gt, Tds>::
bounded_side(const Point &p0, const Point &p1,
	     const Point &p2, const Point &p) const
{
  // return position of point p with respect to triangle p0p1p2, depends on the orientation of the vertives
  CGAL_triangulation_precondition( orientation(p0, p1, p2) != COLLINEAR);
  Orientation o1 = orientation(p0, p1, p),
              o2 = orientation(p1, p2, p),
              o3 = orientation(p2, p0, p);
    
  if (o1 == COLLINEAR){
    if (o2 == COLLINEAR ||  o3 == COLLINEAR) return ON_BOUNDARY;
    if (o2 == ON_NEGATIVE_SIDE && o3 == ON_NEGATIVE_SIDE)        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o2 == COLLINEAR){
    if (o3 == COLLINEAR)                     return ON_BOUNDARY;
    if (o3 == ON_NEGATIVE_SIDE && o1 == ON_NEGATIVE_SIDE)        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o3 == COLLINEAR){
    if (o1 == ON_NEGATIVE_SIDE && o2 == ON_NEGATIVE_SIDE)        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  // from here none ot, o1, o2 and o3 are known to be non null
    if (o1 == o2 && o2 == o3)  return ON_BOUNDED_SIDE;
    return ON_UNBOUNDED_SIDE;
}


template <class Gt, class Tds >
Oriented_side
Triangulation_on_sphere_2<Gt, Tds>::
oriented_side(Face_handle f, const Point &p) const
{
  CGAL_triangulation_precondition ( dimension()==2); 
  return oriented_side(f->vertex(0)->point(),
		       f->vertex(1)->point(),
		       f->vertex(2)->point(),
		       p);
}
/*
template < class Gt, class Tds >
Oriented_side
Triangulation_on_sphere_2<Gt,Tds>::
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
*/

template <class Gt, class Tds >
bool
Triangulation_on_sphere_2<Gt, Tds>::
xy_equal(const Point& p, const Point& q) const
{
  return geom_traits().coradial_sphere_2_object()(p,q);
}

template <class Gt, class Tds >
bool
Triangulation_on_sphere_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{
  // return true if r lies inside the cone defined by trait.sphere, p and q
  return geom_traits().inside_cone_2_object()(p,q,r);
}
//------------------------------------------------------------------------------DEBUG-------------------------------------------------

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::
show_all() const
{
  std::cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<std::endl;
  std::cerr << std::endl<<"====> "<< this; 
  std::cerr <<  " dimension " <<  dimension() << std::endl;
  std::cerr << "nb of vertices " << number_of_vertices() << std::endl;
    
  if (dimension() < 1) return;
  if(dimension() == 1) {
    std::cerr<<" all edges dim 1 "<<std::endl; 
    Edges_iterator aeit;
      for(aeit = edges_begin(); aeit != edges_end(); aeit++){
       show_face(aeit->first);
    }
    return;
  }
  
  std::cerr<<" faces "<<std::endl;
  Faces_iterator fi;
  for(fi = faces_begin(); fi != faces_end(); fi++) {
    show_face(fi);
	  std::cerr<<"   ------------   " <<std::endl;
  }

  
  if (number_of_vertices()>1) {
    std::cerr << "affichage des sommets de la triangulation"
	      <<std::endl;
    Vertices_iterator vi;
    for( vi = vertices_begin(); vi != vertices_end(); vi++){
      show_vertex(vi);
      std::cerr << "  / face associee : "
	     << (void*)(&(*(vi->face())))<< std::endl;;
      }
      std::cerr<<std::endl;
  }
  return;
}
	

template <class Gt, class Tds >
int
Triangulation_on_sphere_2<Gt, Tds>::
number_of_negative_faces() 
{
  int nb=0;
  for(Faces_iterator it=faces_begin();it!=faces_end();++it)
    {if(it->is_negative())nb++;}
  return nb;
}

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::
show_vertex(Vertex_handle vh) const
{
 
  std::cerr << vh->point() << "\t";
  return;
}

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::
show_face(Face_handle fh) const
{
  std::cerr << "face : "<<(void*)&(*fh)<<" => "<<std::endl;
  if(fh->is_negative()) std::cerr << "negative "<<std::endl;
 
  int i = fh->dimension(); 
  switch(i){
  case 0:
    std::cerr <<"point :" ; show_vertex(fh->vertex(0));
    std::cerr <<" / voisin " << &(*(fh->neighbor(0)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
    std::cerr <<"]"  << std::endl;
    break;
  case 1:
     std::cerr <<"point :" ; show_vertex(fh->vertex(0));
     std::cerr <<" / voisin " << &(*(fh->neighbor(0)));
     std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
     std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(1));
     std::cerr <<"]" <<std::endl;

     std::cerr <<"point :" ; show_vertex(fh->vertex(1));
     std::cerr <<" / voisin " << &(*(fh->neighbor(1)));
     std::cerr <<"[" ; show_vertex(fh->neighbor(1)->vertex(0));
     std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(1));
     std::cerr <<"]" <<std::endl;
     break;
  case 2:
    std::cerr <<"point :" ; show_vertex(fh->vertex(0));
    std::cerr <<" / voisin " << &(*(fh->neighbor(0)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(2));
    std::cerr <<"]" <<std::endl;

    std::cerr <<"point :" ; show_vertex(fh->vertex(1));
    std::cerr <<" / voisin " << &(*(fh->neighbor(1)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(1)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(2));
    std::cerr <<"]" <<std::endl;

    std::cerr <<"point :" ; show_vertex(fh->vertex(2));
    std::cerr <<" / voisin " << &(*(fh->neighbor(2)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(2)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(2)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(2)->vertex(2));
    std::cerr <<"]" <<std::endl;
    break;
  }
  return;
	
	
		
	
	
}


 /*---------------------------------------------------------------------PUBLIC REMOVE-------------------------------------*/

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::
make_hole ( Vertex_handle v, std::list<Edge> & hole)
{                
  std::list<Face_handle> to_delete;

  Face_handle  f, fn;
  int i, in ;
  Vertex_handle  vv;
      
  Face_circulator fc = incident_faces(v);
  Face_circulator done(fc);
  do {
    f = fc; fc++;
    i = f->index(v);
    fn = f->neighbor(i);
    in = fn->index(f);
    vv = f->vertex(cw(i));
    if( vv->face()==  f) vv->set_face(fn);
    vv = f->vertex(ccw(i));
    if( vv->face()== f) vv->set_face(fn);
    fn->set_neighbor(in, Face_handle());
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
inline
void
Triangulation_on_sphere_2<Gt, Tds>::

delete_face(Face_handle f)
{
  _tds.delete_face(f);
}

template <class Gt, class Tds >    
inline
void
Triangulation_on_sphere_2<Gt, Tds>::
delete_vertex(Vertex_handle v)
{
  _tds.delete_vertex(v);
}

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2,
	 Face_handle f3, int i3)
{
  return _tds.create_face(f1, i1, f2, i2, f3, i3);
}


template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2)
{
  return _tds.create_face(f1, i1, f2, i2);
}  

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face()
{
  return _tds.create_face();
}

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Face_handle f, int i, Vertex_handle v)
{
  return _tds.create_face(f, i, v);
}

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  return _tds.create_face(v1, v2, v3);
}

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
	    Face_handle f1, Face_handle f2,  Face_handle f3)
{
  return _tds.create_face(v1, v2, v3, f1, f2, f3);
}

template <class Gt, class Tds >    
inline
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
create_face(Face_handle fh)
{
  return _tds.create_face(fh);
}

//----------------------------------------------------------------I/O----------------------------------------------------------------//

template <class Gt, class Tds >
void
Triangulation_on_sphere_2<Gt, Tds>::
file_output(std::ostream& os) const
{
  _tds.file_output(os, Vertex_handle(), true);
}

template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt, Tds>::Vertex_handle
Triangulation_on_sphere_2<Gt, Tds>::
file_input(std::istream& is)
{
  clear();
  Vertex_handle v= _tds.file_input(is, true);
  return v;
}

template <class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Triangulation_on_sphere_2<Gt, Tds> &tr)
{
  tr.file_output(os);
  return os ;
}


template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_on_sphere_2<Gt, Tds> &tr)
{
  tr.file_input(is);
  CGAL_triangulation_assertion(tr.is_valid());
  return is;
}


 

/*--------------------------------------------------------------------------------2D     PROTOTYPES--------------------------------------------------*/
  
    /*

 


 // GEOMETRIC FEATURES AND CONSTRUCTION
  Triangle triangle(Face_handle f) const;
  Segment segment(Face_handle f, int i) const;
  Segment segment(const Edge& e) const;
  Segment segment(const Edge_circulator& ec) const;
  Segment segment(const All_edges_iterator& ei) const;
  Segment segment(const Finite_edges_iterator& ei) const;
  Point circumcenter(Face_handle  f) const; 
  Point circumcenter(const Point& p0, 
		     const Point& p1, 
		     const Point& p2) const;
  

  //MOVE - INSERTION - DELETION - Flip
public:

  


//   template < class InputIterator >
//   int insert(InputIterator first, InputIterator last);

 
  void remove_first(Vertex_handle  v);
  void remove_second(Vertex_handle v);




  // POINT LOCATION

  Face_handle


  Face_handle
  march_locate_2D_LFC(Face_handle start,
		  const Point& t,
		  Locate_type& lt,
		  int& li) const;

  void
  compare_walks(const Point& p,
		Face_handle c1, Face_handle c2,
		Locate_type& lt1, Locate_type& lt2,
		int li1, int li2) const;

  Face_handle
  locate(const Point& p,
	 Locate_type& lt,
	 int& li,
	 Face_handle start = Face_handle()) const;

  Face_handle
  locate(const Point &p,
	 Face_handle start = Face_handle()) const;
    

  
 

  // IO
// template < class Stream >
// Stream&  draw_triangulation(Stream& os) const;

 //PREDICATES

 Oriented_side
 side_of_oriented_circle(Face_handle f, const Point & p) const; 
*/




  


  
  
} //namespace CGAL
    

#endif //CGAL_TRIANGULATION_OM_SPHERE_2_H

