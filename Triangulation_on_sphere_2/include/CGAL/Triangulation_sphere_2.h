#ifndef CGAL_TRIANGULATION_SPHERE_2_H
#define CGAL_TRIANGULATION_SPHERE_2_H


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
#include <CGAL/Triangulation_face_base_sphere_2.h>
#include <CGAL/Random.h>

#include <CGAL/spatial_sort.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/squared_distance_3.h>
#include <boost/type_traits.hpp>
#include <boost/type_traits/integral_constant.hpp>

// this class just provides some basic methods and can not be used independently. No insertion or remove implemented.

namespace CGAL {

template < class Gt, class Tds > class Triangulation_sphere_2;
template < class Gt, class Tds > std::istream& operator>>
    (std::istream& is, Triangulation_sphere_2<Gt,Tds> &tr);
template < class Gt, class Tds >  std::ostream& operator<<
  (std::ostream& os, const Triangulation_sphere_2<Gt,Tds> &tr);
  

template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                             Triangulation_vertex_base_2<Gt>,
                             Triangulation_face_base_sphere_2<Gt> > >
	

	
	
class Triangulation_sphere_2
  : public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> <>
               (std::istream& is, Triangulation_sphere_2 &tr);
  typedef Triangulation_sphere_2<Gt,Tds>             Self;
	

public:
  typedef Tds                                     Triangulation_data_structure;
  typedef Gt                                      Geom_traits;
	
  typedef typename Geom_traits::Point_2           Point;
  typedef typename Geom_traits::Orientation_2     Orientation_2;
  typedef typename Geom_traits::Coradial_sphere_2 Coradial_sphere_2;
  typedef typename Geom_traits::Compare_x_2       Compare_x;
  typedef typename Geom_traits::Compare_y_2       Compare_y;
  typedef typename Geom_traits::FT                FT;
  typedef typename Tds::size_type                 size_type;
  typedef typename Tds::difference_type           difference_type;
  typedef typename Tds::Vertex                    Vertex;
  typedef typename Tds::Face                      Face;
  typedef typename Tds::Edge                      Edge;
  typedef typename Tds::Vertex_handle            Vertex_handle;
  typedef typename Tds::Face_handle              Face_handle;
  typedef typename Tds::Face_circulator         Face_circulator;
  typedef typename Tds::Vertex_circulator       Vertex_circulator;
  typedef typename Tds::Edge_circulator         Edge_circulator;
  typedef typename Tds::Face_iterator           All_faces_iterator;
  typedef typename Tds::Edge_iterator           All_edges_iterator;
  typedef typename Tds::Vertex_iterator         All_vertices_iterator;

	
// This class is used to generate the Solid*_iterators.
class Ghost_tester
{
  const Triangulation_sphere_2 *t;
  public:
  Ghost_tester() {}
  Ghost_tester(const Triangulation_sphere_2 *tr)	  : t(tr) {}
  bool operator()(const All_faces_iterator & fit ) const {
   return fit->is_ghost();
  }
 bool operator()(const All_edges_iterator & eit) const {
	 //int dim = t->dimension();
   Face_handle f = eit->first;
	  bool edge1 = f->is_ghost();
	 Face_handle f2 = f->neighbor(eit->second);
	 //bool edge2b = f2->is_ghost();
   bool edge2 = (f->neighbor(eit->second))->is_ghost();
   bool result = edge1&&edge2;
   return !result;
 }
};
	
class Contour_tester
{
 const Triangulation_sphere_2 *t;
 public:
	Contour_tester() {}
	Contour_tester(const Triangulation_sphere_2 *tr)	  : t(tr) {}
	bool operator() (const All_edges_iterator & eit) const {
	 Face_handle f = eit->first;
	 bool edge1 = f->is_ghost();
	 bool edge2 = f->neighbor(eit->second)->is_ghost();
	return edge1 !=edge2;
	}
};
			
		
//We derive in order to add a conversion to handle.
class Solid_faces_iterator
  : public Filter_iterator<All_faces_iterator, Ghost_tester> 
{
 typedef Filter_iterator<All_faces_iterator, Ghost_tester> Base;
 typedef Solid_faces_iterator                           Self;
 public:
  Solid_faces_iterator() : Base() {}
 Solid_faces_iterator(const Base &b) : Base(b) {}
 Self & operator++() { Base::operator++(); return *this; }
 Self & operator--() { Base::operator--(); return *this; }
 Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
 Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
 operator const Face_handle() const { return Base::base(); }
};
	
	
typedef Filter_iterator<All_edges_iterator, Ghost_tester> 	Solid_edges_iterator;  //Solid edges : both adjacent faces are solid
typedef Filter_iterator<All_edges_iterator, Contour_tester> Contour_edges_iterator; //one solid and one ghost face adjacent to this face
	
  enum Locate_type {VERTEX=0, 
		    EDGE, //1
		    FACE, //2
		    OUTSIDE_CONVEX_HULL, //3
	        OUTSIDE_AFFINE_HULL,//4 
	        CONTOUR,//5
	        NOT_ON_SPHERE,
	        TOO_CLOSE};
	
	
protected:

  Gt  _gt;
  Tds _tds;
 Face_handle _ghost;		//stores an arbitary ghost-face
 double _minDistSquared;   //minimal distance of two points to each other
  double _minRadiusSquared;//minimal distance of a point from center of the sphere
 double _maxRadiusSquared; //maximal distance of a point from center of the sphere
 mutable Random rng;		//used to decide how tostart the march_locate
  


public:

  // CONSTRUCTORS
  Triangulation_sphere_2(const Geom_traits& geom_traits=Geom_traits());
  Triangulation_sphere_2(const Point& sphere); 
  Triangulation_sphere_2(const Triangulation_sphere_2<Gt,Tds> &tr);  
   void clear();
	
private:
	void init(double radius);
	
public:
	
//setting function for the radius of the sphere. Attention: the triangulation is cleared in this step!	
void set_radius(double radius){
		clear();
		init(radius);
}

  //Assignement
  Triangulation_sphere_2 &operator=(const Triangulation_sphere_2 &tr);

  //Helping
  void copy_triangulation(const Triangulation_sphere_2 &tr);
  void swap(Triangulation_sphere_2 &tr);
 
  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const; 
  double squared_distance(const Point& p, const Point& q)const;	
	
	

  // TESTS
  bool is_edge(Vertex_handle va, Vertex_handle vb) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr,int & i) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle &fr) const;
 
   //ACCESS FUNCTION
  const Geom_traits& geom_traits() const { return _gt;}
  const Tds & tds() const  { return _tds;}
  Tds & tds()   { return _tds;}
    
 
  int dimension() const { return _tds.dimension();}
  size_type number_of_vertices() const {return _tds.number_of_vertices();}
size_type number_of_faces() const{return _tds.number_of_faces();}//total number of faces (solid + ghost)
		
  int number_of_ghost_faces();
  
  
  //-----------------------------------------------------------------LOCATION-------------------------------------------------
  Face_handle march_locate_2D(Face_handle c, const Point& t,Locate_type& lt, int& li) const; 
  Face_handle march_locate_1D(const Point& t, Locate_type& lt, int& li) const ;
  Face_handle locate(const Point& p, Locate_type& lt, int& li, Face_handle start) const;
  Face_handle locate(const Point &p, Face_handle start) const;
  Face_handle locate_edge(const Point& p, Locate_type& lt, int& li, bool plane)const;
  void test_distance( const Point& p, Face_handle f, Locate_type &lt, int &li)const;
  bool is_too_close(const Point& p, const Point& q)const;
  //------------------------------------------------------------------------PREDICATES----------------------------------------
  Orientation orientation(const Point& p, const Point& q, const Point& r) const;
  Orientation orientation(const Face_handle f) const;
  Orientation orientation(const Face_handle f, const Point& p)const;
  Orientation orientation(const Point&p, const Point& q, const Point& r, const Point &s) const;
  Comparison_result compare_xyz(const Point& p, const Point& q) const;
  bool equal (const Point& p, const Point& q) const;
  Oriented_side oriented_side(Face_handle f, const Point &p) const;
  bool xy_equal(const Point& p, const Point& q) const;
  bool collinear_between(const Point& p, const Point& q, const Point& r) const;
  //Orientation coplanar_orientation(const Point& p, const Point& q,const Point& r ) const;
Orientation coplanar_orientation(const Point& p, const Point& q,const Point& r, const Point& s ) const;
  //------------------------------------------------------------------DEBUG---------------------------------------------------
  void show_all() const;
  void show_vertex(Vertex_handle vh) const;
  void show_face(Face_handle fh) const;

  //----------------------------------------------------------Creation---------------------------------------------------
 void make_hole(Vertex_handle v, std::list<Edge> & hole);
 Face_handle create_face(Face_handle f1, int i1,Face_handle f2, int i2,
			            Face_handle f3, int i3);
 Face_handle create_face(Face_handle f1, int i1, Face_handle f2, int i2);
 Face_handle create_face();
 Face_handle create_face(Face_handle f, int i, Vertex_handle v);
 Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3);
 Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3,
            			  Face_handle f1, Face_handle f2, Face_handle f3);
 Face_handle create_face(Face_handle);
 void delete_face(Face_handle f);
 void delete_vertex(Vertex_handle v);
	
//-----------------------GEOMETRIC FEATURES AND CONSTRUCTION---------------------------
Point circumcenter(Face_handle  f) const; 
Point circumcenter(const Point& p0, 
					   const Point& p1, 
					   const Point& p2) const;
	
	
	
  //           IN/OUT
  Vertex_handle file_input(std::istream& is);
  void file_output(std::ostream& os) const;

	
	//--------------------------------------------------------------TRAVERSING : ITERATORS AND CIRCULATORS------------------------------------------- 
All_faces_iterator all_faces_begin() const {
  return _tds.faces_begin();
}
	
All_faces_iterator all_faces_end() const {
  return _tds.faces_end();
}
	
Solid_faces_iterator solid_faces_begin() const {
 if (dimension() < 2)
	return solid_faces_end();
return CGAL::filter_iterator( all_faces_end(),
							 Ghost_tester(this), all_faces_begin() );
} 
		
	
Solid_faces_iterator solid_faces_end() const {
	return CGAL::filter_iterator(  all_faces_end(),
     							 Ghost_tester(this)   );
}
	
Solid_edges_iterator solid_edges_begin() const {
		if ( dimension() < 1 )
			return solid_edges_end();
		return CGAL::filter_iterator (all_edges_end(), Ghost_tester(this),
									  all_edges_begin());
	}
	
Solid_edges_iterator solid_edges_end() const {
	return CGAL::filter_iterator (all_edges_end(), Ghost_tester(this));
}
	
Contour_edges_iterator contour_edges_begin() const{
 if(dimension()<1)
	return contour_edges_begin();
 return CGAL::filter_iterator (all_edges_end(), Ghost_tester(this),
									  all_edges_begin());
}
	
Contour_edges_iterator contour_edges_end() const{
	return CGAL::filter_iterator (all_edges_end(), Contour_tester(this));
}

All_vertices_iterator vertices_begin() const{
	return _tds.vertices_begin();
}

All_vertices_iterator vertices_end() const {
	return _tds.vertices_end();
}
	
All_edges_iterator all_edges_begin() const{
	return _tds.edges_begin();
}
	
All_edges_iterator all_edges_end() const{
	return _tds.edges_end();
}
	
Face_circulator incident_faces( Vertex_handle v, Face_handle f = Face_handle()) const{
	return _tds.incident_faces(v,f);
}
	
Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f = Face_handle()) const{
	return _tds.incident_vertices(v,f);
}
	
	
Edge_circulator incident_edges(Vertex_handle v, Face_handle f = Face_handle()) const{
	return _tds.incident_edges(v,f);
}
	
	
size_type degree(Vertex_handle v) const {
	return _tds.degree(v);
}
	
Vertex_handle mirror_vertex(Face_handle f, int i) const{
	return _tds.mirror_vertex(f,i);
}
	
int mirror_index(Face_handle v, int i) const{
	return _tds.mirror_index(v,i);
}
	
Edge mirror_edge(const Edge e) const
{
 return _tds.mirror_edge(e);
}
	
/*---------------------------------------------------------------------TEMPLATE MEMBERS--------------------------------------*/
 public: 
 
 template<class FaceIt>
 void delete_faces(FaceIt face_begin, FaceIt face_end)
 {
   FaceIt fit=face_begin;
   for(;fit!=face_end;++fit)
       delete_face(*fit);    
 }
	
//is_on_sphere test whether a given point lies close enough to the sphere. Whether a test is necessary or not is defined in the traits (requires test)	
//additional test needed	
private:
void
is_on_sphere(boost::true_type, const Point &p, Locate_type &lt) const {
 double distance2 = pow(p.x(),2)+pow(p.y(),2)+pow(p.z(),2);
 bool test = _minRadiusSquared<distance2&& distance2<_maxRadiusSquared;
 if (!test)
			lt = NOT_ON_SPHERE;		
}
			
// no test needed -> empty function
void
is_on_sphere(boost::false_type, const Point &p, Locate_type &lt) const{
}
			

	
};
	
// CONSTRUCTORS
	
template <class Gt, class Tds >
Triangulation_sphere_2<Gt, Tds>::
Triangulation_sphere_2(const Geom_traits& geom_traits) 
: _gt(geom_traits), _tds()
 {init(1);}
	
template <class Gt, class Tds >
Triangulation_sphere_2<Gt, Tds>::
Triangulation_sphere_2(const Point& sphere) 
: _gt(sphere), _tds() 	
 { init(1);}	
	
	
//initializes the requiered data to proof the preconditons for the lemma about hidden vertices. By default radius ==1;
template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::
init(double radius){
	double minRadius = radius*(1-pow(2, -50));
	_minRadiusSquared = minRadius * minRadius;
	double maxRadius = radius*(1+pow(2, -50));
	_maxRadiusSquared = maxRadius * maxRadius;
	double minDist = radius*pow (2, -23);
	_minDistSquared = minDist *minDist;
	_gt.set_radius(radius);
	
}
	

// copy constructor duplicates vertices and faces
template <class Gt, class Tds >
Triangulation_sphere_2<Gt, Tds>::
Triangulation_sphere_2(const Triangulation_sphere_2 &tr) 
  : _gt(tr._gt), _tds(tr._tds)
	{init(_gt._radius);}

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>:: 
clear()
{
  _tds.clear(); 
 
}
  
  
template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::   
copy_triangulation(const Triangulation_sphere_2 &tr)
{
  _tds.clear();
  _gt = tr._gt;
  _tds = tr._tds;
   init(_gt._radius);
}

  //Assignement
template <class Gt, class Tds >
Triangulation_sphere_2<Gt, Tds> &
Triangulation_sphere_2<Gt, Tds>::
operator=(const Triangulation_sphere_2 &tr)
{
  copy_triangulation(tr);
  return *this;
}

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>:: 
swap(Triangulation_sphere_2 &tr)
{
  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t; 
}

//--------------------------CHECKING---------------------------

//is_valid is used by in the instream method, Delaunay_tri...
	//has its own is_valid

template <class Gt, class Tds >
bool
Triangulation_sphere_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  if (dimension() <= 0 ||
     (dimension()==1 && number_of_vertices() == 3 ) ) return result;
 
  if (dimension() == 1) {
    All_vertices_iterator vit=vertices_begin();
    }    

  else { //dimension() == 2

    for(All_faces_iterator it=all_faces_begin(); 
	it!=all_faces_end(); it++) {
      Orientation s = orientation(it->vertex(0)->point(),
				  it->vertex(1)->point(),
				  it->vertex(2)->point());
      CGAL_triangulation_assertion( s == LEFT_TURN || it->is_ghost() );
      result = result && ( s == LEFT_TURN || it->is_ghost() );
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
Triangulation_sphere_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb) const
{
  return _tds.is_edge( va, vb);
}

template <class Gt, class Tds >
inline bool
Triangulation_sphere_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i) const
{
  return _tds.is_edge(va, vb, fr, i);
}

template <class Gt, class Tds >
inline bool 
Triangulation_sphere_2<Gt, Tds>::
is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
{
  return _tds.is_face(v1, v2, v3);
}

template <class Gt, class Tds >
inline bool 
Triangulation_sphere_2<Gt, Tds>::
is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,Face_handle &fr) const
{
  return _tds.is_face(v1, v2, v3, fr);
}

//---------------------------------------------------------------------------/POINT LOCATION---------------------------------------//
// tests whether the two points p and q are too close according to the lemma about hidden vertices.	
template<class Gt, class Tds>
inline bool
Triangulation_sphere_2<Gt, Tds> ::
is_too_close(const Point& p, const Point& q)const{
	 return squared_distance(p,q)<=_minDistSquared;
}

/*location for degenerated cases: locates the conflicting edge in a 1 dimensional triangulation.
 This methode is used, when the new point is coplanar with the existing 	vertices.
	bool plane defines whether the points are also coplanar with the center of the sphere (true) or not (false).
 */
	
template <class Gt, class Tds>
typename Triangulation_sphere_2<Gt, Tds> ::Face_handle
Triangulation_sphere_2<Gt, Tds>::
locate_edge(const Point& p, Locate_type& lt, int& li, bool plane)const
{
  Face_handle loc;
  if(plane){
	All_edges_iterator eit;
	for(eit = all_edges_begin(); eit !=all_edges_end(); eit ++){
		if(!eit->first->is_ghost())
		  if(collinear_between(eit->first->vertex(0)->point(), eit->first->vertex(1)->point(),p)){
			test_distance( p, (*eit).first, lt, li);
			return (*eit).first;
		  }
			
		if (eit->first->is_ghost())
			loc = eit->first;
		}//end for
	  
	 test_distance(p, loc, lt, li);
	return loc;
	}
	
  else {//not plane
   All_edges_iterator eit;
   Face_handle f;
   for(eit = all_edges_begin(); eit!=all_edges_end(); eit ++){
	 f=eit->first;
	 Vertex_handle v1 = f->vertex(0);
	 Vertex_handle v2 = f -> vertex(1);
	  if(orientation(v1->point(), v2->point(), p)==RIGHT_TURN){
		lt=EDGE;
		li=2;
		test_distance( p, (*eit).first, lt, li);
		return (*eit).first;
	   }	
	}//end for
		
	test_distance(p, loc, lt, li);
	return loc;
  }//end else
}
	
/*
 calls too_close for possible conflicts.
 If the point p is too close to an existing vertex, this vertex is returned.
 should be replaced by a nearest neighbor search.
 */
template <class Gt, class Tds > 
inline
void Triangulation_sphere_2<Gt, Tds>::
test_distance( const Point& p, Face_handle f, Locate_type& lt, int& li)const{
CGAL_precondition(dimension()>=-1);
 switch (dimension()){
 case -1 : {	
	 if(is_too_close(p, vertices_begin()->point())){
		lt = TOO_CLOSE;
		li=0;
	   return;
	 }
 } break;
		
 case 0:
 case 1:{
	 
	 All_vertices_iterator vi=vertices_begin(); 
  	 for(;vi!=vertices_end();vi++)
	 if(is_too_close(vi->point(),p)){
		 lt = TOO_CLOSE;
		 li=1;
		 return;
	 }
		
 } break;
		
  case 2:{
	Vertex_handle v0 = f->vertex(0);
	Vertex_handle v1= f->vertex(1);
	Vertex_handle v2 = f->vertex(2);
	 
	if(is_too_close(v0->point(),p)){
	  lt = TOO_CLOSE;
	  li = 0;
		return;
	}
	if(is_too_close(v1->point(),p)){
		  lt = TOO_CLOSE;
		  li = 1;
		  return;
	  }
	if(is_too_close(v2->point(),p)){
		  lt = TOO_CLOSE;
		  li = 2;
		  return;
	  }
  }break;
  }
}	
	
		
template <class Gt, class Tds >    
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
march_locate_1D(const Point& p, Locate_type& lt, int& li) const
{	
  Face_handle f =all_edges_begin()->first;
 //check if p is coplanar with existing points
	
 //first three points of triangulation
 Vertex_handle v1=f->vertex(0);
 Vertex_handle v2=f->vertex(1);
 Vertex_handle v3=f->neighbor(0)->vertex(1);
 
 Orientation orient = orientation(v1->point(), v2->point(), v3->point(),p);
 if(orient !=ON_ORIENTED_BOUNDARY){
   lt = OUTSIDE_AFFINE_HULL;
   li = 4;
	test_distance(p, f, lt, li);
	return f;
 }
	
 //check if p is coradial with one existing point
 All_vertices_iterator vi;
 for( vi = vertices_begin(); vi != vertices_end(); vi++){
	if (xy_equal(vi->point(), p)){
	  lt = VERTEX;
	  li = 1;
	 return f;
	}
  }
	
//***find conflicting edge***
 lt=EDGE;
 li = 2;
 Orientation pqr = orientation(v1->point(),v2->point(),v3->point());
 if(pqr == ON_ORIENTED_BOUNDARY)
	return locate_edge(p, lt, li, true);
	else 
		return locate_edge(p, lt,li,false);
}
	

	  
template <class Gt, class Tds >   typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
march_locate_2D(Face_handle c,const Point& t,Locate_type& lt,int& li) const
{
	
 CGAL_triangulation_assertion(!c->is_ghost());
 Face_handle prev = Face_handle();
 bool first=true;
 while (1) {
  if ( c->is_ghost() ) {
 	if(orientation(c, t)==ON_POSITIVE_SIDE){ //conflict with the corresponding face
	  lt = OUTSIDE_CONVEX_HULL;
	  li = 4;
	  test_distance(t,c, lt, li); 
	  return c;
	 }
      else {//one of the neighbour face has to be in conflict with
		Face_handle next = Face_handle();
		for(int i=0; i<=2 ; i++){
			next=c->neighbor(i);
			//Orientation orient =orientation(next, t);
		   if(orientation(next, t)==ON_POSITIVE_SIDE){
		    	lt = CONTOUR;
				li = 4;
				test_distance(t, next, lt, li);
			    return next;
			}
		 }
CGAL_triangulation_precondition(false);
				
	   }
	}// end if ghost
	
		
	
	const Point & p0 = c->vertex( 0 )->point();
	const Point & p1 = c->vertex( 1 )->point();
	const Point & p2 = c->vertex( 2 )->point();
		
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
    case -2:{
	 lt = FACE;
	 li = 4;
	 break;}
    case 1:{ 
	lt = EDGE;
	li = ( o0 == COLLINEAR ) ? 2 :
	  ( o1 == COLLINEAR ) ? 0 : 1;
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
	 test_distance(t,c,lt,li);
   return c;
    
  }
} 
	
template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt,Tds>::
locate(const Point& p,Locate_type& lt,int& li, Face_handle start) const
{
	is_on_sphere( typename Gt::requires_test(), p, lt);
	
	if(lt == NOT_ON_SPHERE)
		return Face_handle();
	
  switch (dimension()){
	case-2 : {         //empty triangulation
      lt = OUTSIDE_AFFINE_HULL;
      li = 4; // li should not be used in this case
      return start;
    }
  
	case -1: {      //1 vertex
     Point q=vertices_begin()->point();
     if(xy_equal(q,p)){
      lt=VERTEX;
	  li=0;
     } else{
      lt=OUTSIDE_AFFINE_HULL;
	  li=4;
	  test_distance(p, all_faces_begin(), lt, li);
    }
  	//test_distance(p, all_faces_begin(), lt, li);
    return Face_handle();
		//return start;
   }
 
	case 0: {
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
	}else{
      lt=OUTSIDE_AFFINE_HULL;
      li=4;  
    }
	test_distance(p, f, lt, li);
    return Face_handle();
  }
  
	case 1:{
     return march_locate_1D(p, lt, li);
     }
  }
   
  if(start==Face_handle()){
      start=all_faces_begin();
  }
 
  if(start->is_ghost()){
	for (All_faces_iterator it = this->_tds.face_iterator_base_begin(); it !=all_faces_end(); it++) {  
	   if(!it->is_ghost()){
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
typename Triangulation_sphere_2<Gt, Tds>:: Face_handle
Triangulation_sphere_2<Gt, Tds>::
locate(const Point &p,
       Face_handle start) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, start);
}



//------------------------------------------------------------------------------PREDICATES-----------------------------------------------------------------
template <class Gt, class Tds >
Comparison_result
Triangulation_sphere_2< Gt,  Tds>::
compare_xyz(const Point &p, const Point &q) const
{
		return geom_traits().compare_xyz_3_object()(p, q);
}

	
template <class Gt, class Tds >	
bool
Triangulation_sphere_2<Gt, Tds>::
equal(const Point &p, const Point &q) const
{
		return compare_xyz(p, q) == EQUAL;
}
	

template <class Gt, class Tds >
inline
Orientation
Triangulation_sphere_2<Gt, Tds>::
coplanar_orientation(const Point& p, const Point& q,const Point& r, const  Point &s ) const
{
	return geom_traits().orientation_1_object()(p,q,r,s);
}
	
	
template <class Gt, class Tds >
inline
Orientation
Triangulation_sphere_2<Gt, Tds>::
orientation(const Point& p, const Point& q,const Point& r ) const
{
  return geom_traits().orientation_2_object()(p,q,r);
}

template <class Gt, class Tds >
inline
Orientation
Triangulation_sphere_2<Gt, Tds>::
orientation(const Face_handle fh,const Point& r ) const
{
	return orientation(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point(),r);
}	
	

template <class Gt, class Tds >
Orientation
Triangulation_sphere_2<Gt, Tds>::
orientation(const Face_handle f) const
{
  return  orientation(f->vertex(0)->point(),f->vertex(1)->point(),f->vertex(2)->point());
}


template <class Gt, class Tds>
Orientation
Triangulation_sphere_2<Gt, Tds>::
orientation(const Point&p, const Point &q, const Point &r, const Point & s)const
{
	return geom_traits().orientation_2_object()(p,q,r,s);
}
	
//returns true if p,q, and O (center of the sphere) are aligned and if p (q) is located in the segment Oq (Op)		
template <class Gt, class Tds >
bool
Triangulation_sphere_2<Gt, Tds>::
xy_equal(const Point& p, const Point& q) const
{
  return geom_traits().coradial_sphere_2_object()(p,q);
}

// return true if r lies inside the cone defined by trait.sphere, p and q
template <class Gt, class Tds >
bool
Triangulation_sphere_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{  
  return geom_traits().inside_cone_2_object()(p,q,r);
}
//------------------------------------------------------------------------------DEBUG-------------------------------------------------

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::
show_all() const
{
	//Triangulation_2::show_all();
  std::cerr<< "AFFICHE TOUTE LA TRIANGULATION :"<<std::endl;
  std::cerr << std::endl<<"====> "<< this; 
  std::cerr <<  " dimension " <<  dimension() << std::endl;
  std::cerr << "nb of vertices " << number_of_vertices() << std::endl;
    
  if (dimension() < 1) return;
  if(dimension() == 1) {
    std::cerr<<" all edges dim 1 "<<std::endl; 
    All_edges_iterator aeit;
      for(aeit =all_edges_begin(); aeit !=all_edges_end(); aeit++){
       show_face(aeit->first);
		std::cerr<<"   ------------   " <<std::endl;  
    }
    return;
  }
  
  std::cerr<<" faces "<<std::endl;
  All_faces_iterator fi;
  for(fi = all_faces_begin(); fi !=all_faces_end(); fi++) {
    show_face(fi);
	  std::cerr<<"   ------------   " <<std::endl;
  }

  
  if (number_of_vertices()>1) {
    std::cerr << "affichage des sommets de la triangulation"
	      <<std::endl;
    All_vertices_iterator vi;
    for( vi = vertices_begin(); vi != vertices_end(); vi++){
      show_vertex(vi);
      std::cerr << "  / face associee : "
	     << (void*)(&(*(vi->face())))<< std::endl;;
      }
      std::cerr<<std::endl;
  }
  return;
}
	
//returns the number of ghost_faces
template <class Gt, class Tds >
int
Triangulation_sphere_2<Gt, Tds>::
number_of_ghost_faces() 
{
  int nb=0;
  for(All_faces_iterator it=all_faces_begin();it!=all_faces_end();++it)
    if(it->is_ghost())
		nb++;
  return nb;
}

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::
show_vertex(Vertex_handle vh) const
{
 
  std::cerr << vh->point() << "\t";
  return;
}

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::
show_face(Face_handle fh) const
{
  std::cerr << "face : "<<(void*)&(*fh)<<" => "<<std::endl;
  if(fh->is_ghost()) std::cerr << "ghost "<<std::endl;
 
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
Triangulation_sphere_2<Gt, Tds>::
make_hole ( Vertex_handle v, std::list<Edge> & hole)
{     
	return this->_tds.make_hole(v, hole);
}
	
template <class Gt, class Tds >    
inline
void
Triangulation_sphere_2<Gt, Tds>::
delete_face(Face_handle f)
{
  _tds.delete_face(f);
}

template <class Gt, class Tds >    
inline
void
Triangulation_sphere_2<Gt, Tds>::
delete_vertex(Vertex_handle v)
{
  _tds.delete_vertex(v);
}

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2,
	 Face_handle f3, int i3)
{
  return _tds.create_face(f1, i1, f2, i2, f3, i3);
}


template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
	 Face_handle f2, int i2)
{
  return _tds.create_face(f1, i1, f2, i2);
}  

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face()
{
  return _tds.create_face();
}

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Face_handle f, int i, Vertex_handle v)
{
  return _tds.create_face(f, i, v);
}

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  return _tds.create_face(v1, v2, v3);
}

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
	    Face_handle f1, Face_handle f2,  Face_handle f3)
{
  return _tds.create_face(v1, v2, v3, f1, f2, f3);
}

template <class Gt, class Tds >    
inline
typename Triangulation_sphere_2<Gt, Tds>::Face_handle
Triangulation_sphere_2<Gt, Tds>::
create_face(Face_handle fh)
{
  return _tds.create_face(fh);
}
	
//----------------CONSTRUCTION-----------------------------------
	

template<class Gt, class Tds>
inline
typename Triangulation_sphere_2<Gt,Tds>::Point
Triangulation_sphere_2<Gt,Tds>::
circumcenter (const Point& p0, const Point& p1, const Point&  p2) const
{
	return 
	geom_traits().construct_circumcenter_2_object()(p0,p1,p2);
}
	
	
template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt, Tds>::Point
Triangulation_sphere_2<Gt, Tds>::
circumcenter(Face_handle  f) const
{
	return circumcenter((f->vertex(0))->point(), 
							(f->vertex(1))->point(), 
							(f->vertex(2))->point());
}
	
template<class Gt, class Tds>
double
Triangulation_sphere_2<Gt, Tds> ::
squared_distance(const Point& p, const Point& q)const
{return geom_traits().compute_squared_distance_3_object()(p,q);}
	

//----------------------------------------------------------------I/O----------------------------------------------------------------//

template <class Gt, class Tds >
void
Triangulation_sphere_2<Gt, Tds>::
file_output(std::ostream& os) const
{
  _tds.file_output(os, Vertex_handle(), true);
}

template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt, Tds>::Vertex_handle
Triangulation_sphere_2<Gt, Tds>::
file_input(std::istream& is)
{
  clear();
  Vertex_handle v= _tds.file_input(is, true);
  return v;
}

template <class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Triangulation_sphere_2<Gt, Tds> &tr)
{
  tr.file_output(os);
  return os ;
}


template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_sphere_2<Gt, Tds> &tr)
{
  tr.file_input(is);
  CGAL_triangulation_assertion(tr.is_valid());
  return is;
}
  
} //namespace CGAL
    

#endif //CGAL_TRIANGULATION_ON_SPHERE_2_H

