#ifndef CGAL_Delaunay_triangulation_sphere_2_H
#define CGAL_Delaunay_triangulation_sphere_2_H

#define HOLE_APPROACH

#include <CGAL/Triangulation_sphere_2.h>
#include <CGAL/Triangulation_face_base_sphere_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/utility.h>
#include <fstream>



namespace CGAL { 


template < class Gt,
           class Tds  = Triangulation_data_structure_2 <
             Triangulation_vertex_base_2<Gt>,
		     Triangulation_face_base_sphere_2<Gt> > >
class Delaunay_triangulation_sphere_2
  : public Triangulation_sphere_2<Gt,Tds>
{

   typedef Delaunay_triangulation_sphere_2<Gt, Tds>                         Self;
   typedef Triangulation_sphere_2<Gt,Tds>                                  Base;
	
   	
public:
  typedef Tds									  Triangulation_data_structure;
  typedef Gt                                      Geom_traits;
  typedef typename Gt::Point_2                    Point;
  typedef typename Gt::Segment_3				  Segment;
  typedef typename Base::size_type                size_type;
  typedef typename Base::Face_handle              Face_handle;
  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Vertex                   Vertex;
  typedef typename Base::Edge                     Edge;
  typedef typename Base::Locate_type              Locate_type;
  typedef typename Base::Face_circulator          Face_circulator;
  typedef typename Base::Edge_circulator		  Edge_circulator;
  typedef typename Base::Vertex_circulator        Vertex_circulator;
  typedef typename Base::All_vertices_iterator    All_vertices_iterator;		
  typedef typename Base::All_edges_iterator       All_edges_iterator;
  typedef typename Base::All_faces_iterator       All_faces_iterator;
  typedef typename Base::Solid_faces_iterator     Solid_faces_iterator;
  typedef typename Base::Solid_edges_iterator     Solid_edges_iterator;
  typedef typename Base::Contour_edges_iterator Contour_edges_iterator;
	

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;
  using Base::dimension;
  using Base::geom_traits;
  using Base::create_face;
  using Base::number_of_faces;
  using Base::number_of_vertices;
  using Base::all_faces_begin;
  using Base::all_faces_end;
  using Base::all_edges_begin;
  using Base::all_edges_end;
  using Base::solid_faces_begin;
  using Base::solid_faces_end;
  using Base::solid_edges_begin;
  using Base::solid_edges_end;
  using Base::contour_edges_begin;
  using Base::contour_edges_end;	
  using Base::vertices_begin;	
  using Base::vertices_end;		
  using Base::OUTSIDE_AFFINE_HULL;
  using Base::FACE;
  using Base::OUTSIDE_CONVEX_HULL;
  using Base ::NOT_ON_SPHERE;
  using Base :: VERTEX;
  using Base :: TOO_CLOSE;
  using Base::orientation;
  using Base ::show_all;
  using Base ::show_face;
  using Base ::show_vertex;
  using Base::clear;
  using Base::delete_faces;
  using Base::compare_xyz;
  using Base::coplanar_orientation;
  using Base::set_radius;
  using Base::circumcenter;
	
 #endif
	
//class for sorting points lexicographically.
//This sorting is used for the perturbation in power_test
class Perturbation_order {
  const Self *t;
  public:
	Perturbation_order(const Self *tr)
	: t(tr) {}
		
	bool operator()(const Point *p, const Point *q) const {
     return t->compare_xyz(*p, *q) == SMALLER;
     }
};

public: 
  //CONSTRUCTORS
  Delaunay_triangulation_sphere_2(const Gt& gt=Gt()) 
	:Base(Gt(gt))
  {}


//CHECK
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid_face(Face_handle fh,bool verbose = false, int level = 0 ) const;
  bool is_valid_vertex(Vertex_handle fh, bool verbose = false, int level = 0) const;
  bool is_plane() const;
//checks whether neighboring faces are linked correctly to each other.
  void check_neighboring()
  {
    All_faces_iterator eit;
    if(dimension()==1){
      for(eit=all_faces_begin();eit!=all_faces_end();++eit)
      {
	  Face_handle f1 = eit->neighbor(0);
	  Face_handle f2 = eit->neighbor(1);
	  CGAL_triangulation_assertion(f1->has_neighbor(eit));
	  CGAL_triangulation_assertion(f2->has_neighbor(eit));
      }
      
    }

    for(eit=all_faces_begin();eit!=all_faces_end();++eit)
      {
	  Face_handle f1 = eit->neighbor(0);
	  Face_handle f2 = eit->neighbor(1);
	  Face_handle f3 = eit->neighbor(2);
	  CGAL_triangulation_assertion(f1->has_neighbor(eit));
	  CGAL_triangulation_assertion(f2->has_neighbor(eit));
	  CGAL_triangulation_assertion(f3->has_neighbor(eit));
      }

  }
 

  //INSERTION
  Vertex_handle insert(const Point &p, Locate_type  lt, Face_handle loc, int li );
  Vertex_handle insert(const Point &p, Face_handle f = Face_handle() );
  Vertex_handle insert_first(const Point &p);
  Vertex_handle insert_second(const Point &p);
  Vertex_handle insert_outside_affine_hull_regular(const Point& p);
  Vertex_handle insert_cocircular(const Point &p, Locate_type lt, Face_handle loc);
	
  bool test_conflict(const Point  &p, Face_handle fh) const;
  bool update_ghost_faces(Vertex_handle v=Vertex_handle(), bool first = false);
	
  //REMOVAL
  void remove_degree_3(Vertex_handle v, Face_handle f = Face_handle());
  void remove(Vertex_handle v);
  void remove_1D(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  bool test_dim_down(Vertex_handle v);
  bool test_dim_up(const Point &p)const;
  void fill_hole_regular(std::list<Edge> & hole);

	
  Oriented_side power_test(const Point &p,const Point &q, const Point &r, const Point &s, bool perturb = false) const;
  Oriented_side power_test(const Point &p,const Point &q, const Point &r) const;
  Oriented_side power_test(const Point &p, const Point &r) const;
  Oriented_side power_test(const Face_handle &f,const Point &p, bool perturb = false) const;
  Oriented_side power_test(const Face_handle& f, int i, const Point &p) const;	
	
  //dual
	Point dual (Face_handle f) const;
	Object dual(const Edge &e) const ;
	Object dual(const Edge_circulator& ec) const;
	Object dual(const All_edges_iterator& ei) const;	

	
  //TEMPLATE MEMBERS
  //----------------------------------------------------------------------HOLE APPROACH
	
template <class Stream>
Stream &write_vertices(Stream &out,std::vector<Vertex_handle> &t)
  {
    for(typename std::vector<Vertex_handle>::iterator it= t.begin(); it!= t.end(); ++it)
      {if((*it)->face()==Face_handle()){
	  Point p=(*it)->point();
	  out << p.x() << " " 
	  << p.y() << " " 
	  << p.z() << std::endl;
	}
      }
    return out;
  }
	
template <class Stream>
 Stream &write_triangulation_to_off_2(Stream &out,Stream &out2)
{

    // Points of triangulation
 for (All_faces_iterator it = this->_tds.face_iterator_base_begin(); it !=all_faces_end(); it++) {  
    if(!it->is_ghost()){//assert(orientation(it)==POSITIVE);
	for (int i=0 ; i<3 ; i++) {
	  if(it->vertex(i)!=Vertex_handle())
	  {
	    Point p = it->vertex(i)->point();
	    out << p.x() << " " 
		<< p.y() << " " 
		<< p.z() << std::endl;
	  }
	}
      }
  else{
	for (int i=0 ; i<3 ; i++){
	  if(it->vertex(i)!=Vertex_handle())
	    {
	      Point p = it->vertex(i)->point();
	      out2 << p.x() << " " 
	           << p.y() << " " 
	           << p.z() << std::endl;
	    }	  
	  }
      }
    }
 

  return out;
}

template <class Stream>
Stream &write_triangulation_to_off(Stream &out) 
{
  std::vector<Face_handle> faces;

  // Points of triangulation
  for (All_faces_iterator it =all_faces_begin(); it !=all_faces_end(); it++) {
    for (int i=0 ; i<3 ; i++) {
      Point p = it->vertex(i)->point();
      out << p.x() << " " 
	  << p.y() << " " 
	  << p.z() << std::endl;
    }
  }

  return out;
}
 
template <class Stream >
Stream &write_face_to_off(Stream &out,Face_handle f)
{

      for (int i=0 ; i<3 ; i++) {
	Point p = f->vertex(i)->point();
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
      }
  

  return out;
}


template <class Stream,class FaceIt>
Stream &write_faces_to_off(Stream &out,FaceIt face_begin, FaceIt face_end)
{

    FaceIt fit=face_begin;
    for(;fit!=face_end;++fit)
    {
      for (int i=0 ; i<3 ; i++) {
	Point p = (*fit)->vertex(i)->point();
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
      }
  }

  return out;
}

template <class Stream,class FaceIt>
Stream &write_edges_to_off(Stream &out,FaceIt face_begin, FaceIt face_end)
{

    FaceIt fit=face_begin;
    for(;fit!=face_end;++fit)
    {
      Face_handle f=(*fit).first;
      int i=(*fit).second;
     
      Point p=f->vertex(cw(i))->point();
      Point q=f->vertex(ccw(i))->point();
     
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;

	out << q.x() << " " 
	    << q.y() << " " 
	    << q.z() << std::endl;
    }

  return out;
}
//insert points in a given range using spacial sorting
template < class InputIterator >
int insert(InputIterator first, InputIterator last)
{
	int n = number_of_vertices();
		
	std::vector<Point> points (first, last);
	std::random_shuffle (points.begin(), points.end());
	spatial_sort (points.begin(), points.end());
		
	Face_handle hint;
	Vertex_handle v;
	for (typename std::vector<Point>::const_iterator p = points.begin(),end = points.end(); p != end; ++p){
		v = insert (*p, hint);
		if( v != Vertex_handle())//could happen if the point is not on the sphere
			hint = v->face();
	}
		
	return number_of_vertices() - n;
}

template <class OutputItFaces, class OutputItBoundaryEdges> 
std::pair<OutputItFaces,OutputItBoundaryEdges>
get_conflicts_and_boundary(const Point  &p, OutputItFaces fit, OutputItBoundaryEdges eit, 
						   Face_handle fh ) const
{
	CGAL_triangulation_precondition( dimension() == 2);
	CGAL_triangulation_precondition(test_conflict(p,fh));  
		
	*fit++ = fh; //put fh in OutputItFaces
	fh->set_in_conflict_flag(1);
	std::pair<OutputItFaces,OutputItBoundaryEdges>
	pit = std::make_pair(fit,eit);
	pit = propagate_conflicts(p,fh,0,pit);
	pit = propagate_conflicts(p,fh,1,pit);
	pit = propagate_conflicts(p,fh,2,pit);
					
	return std::make_pair(fit,eit);
} 

	
private:
template <class OutputItFaces, class OutputItBoundaryEdges> 
std::pair<OutputItFaces,OutputItBoundaryEdges>
propagate_conflicts (const Point  &p, Face_handle fh, int i,
						 std::pair<OutputItFaces,OutputItBoundaryEdges> pit)  const 
{
	Face_handle fn = fh->neighbor(i);
	if (fn->get_in_conflict_flag() ==1)
		return pit;
	if (! test_conflict(p,fn)) 
			*(pit.second)++ = Edge(fn, fn->index(fh));
	else {
			*(pit.first)++ = fn;
			fn->set_in_conflict_flag(1);
			int j = fn->index(fh);
			pit = propagate_conflicts(p,fn,ccw(j),pit);
			pit = propagate_conflicts(p,fn,cw(j), pit);
		}
		return pit;
  }

	
};
	
//------------power-test

template < class Gt, class Tds >
Oriented_side
Delaunay_triangulation_sphere_2<Gt,Tds>::
power_test(const Face_handle &f, const Point &p, bool perturb) const
{
	return power_test(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(),p, perturb);
}
	
template < class Gt, class Tds >
Oriented_side
Delaunay_triangulation_sphere_2<Gt,Tds>::
power_test(const Face_handle& f, int i, const Point &p) const
{
	CGAL_triangulation_precondition ( orientation(f->vertex(ccw(i))->point(),
												  f->vertex( cw(i))->point(),  p)
            									 == COLLINEAR);
		
	return  power_test(f->vertex(ccw(i))->point(), f->vertex( cw(i))->point(), p);
}
//computes the power-test of 4 points. perturb defines whether a symbolic perturbation
//is used (by default : perturb == false)
	//in the perturbation the smalest vertex is in conflict with the others
template < class Gt, class Tds >
inline
Oriented_side
Delaunay_triangulation_sphere_2<Gt,Tds>::
power_test(const Point &p0, const Point &p1, const Point &p2, const Point &p, bool perturb) const
{
	Oriented_side os = geom_traits().power_test_2_object()(p0,p1,p2,p);
	if(os!=ON_ORIENTED_BOUNDARY || !perturb)
		return os;
   // We are now in a degenerate case => we do a symbolic perturbation.
		
	// We sort the points lexicographically.
	const Point * points[3] = {&p0, &p1, &p2};
	std::sort(points, points+3, Perturbation_order(this) );
		
		
	if(points[0] == &p0){
	   	if(compare_xyz(p, p0)==SMALLER)
			  return ON_POSITIVE_SIDE;
		if(coplanar_orientation(p0,p1,p,p2)==ON_NEGATIVE_SIDE)
			  return ON_NEGATIVE_SIDE;
		if(coplanar_orientation(p0,p2,p,p1)==ON_NEGATIVE_SIDE)
			   return ON_NEGATIVE_SIDE;
			
	  return ON_POSITIVE_SIDE;
	 }
				
	 if(points[0] == &p1){
		if(compare_xyz(p, p1)==SMALLER)
			return ON_POSITIVE_SIDE;
		if(coplanar_orientation(p1,p0,p,p2)==ON_NEGATIVE_SIDE)
			return ON_NEGATIVE_SIDE;
		if(coplanar_orientation(p1,p2,p,p0)==ON_NEGATIVE_SIDE)
			return ON_NEGATIVE_SIDE;
			
	  return ON_POSITIVE_SIDE;
	 }
		
	 if(points[0] == &p2){
		if(compare_xyz(p, p2)==SMALLER)
	       return ON_POSITIVE_SIDE;
		if(coplanar_orientation(p2,p1,p,p0)==ON_NEGATIVE_SIDE)
			return ON_NEGATIVE_SIDE;
		if(coplanar_orientation(p2,p0,p1,p)==ON_NEGATIVE_SIDE)
			return ON_NEGATIVE_SIDE;
			
	 return ON_POSITIVE_SIDE;		
	}
				
 CGAL_assertion(false);
 return ON_NEGATIVE_SIDE;
								
}
			
	
template < class Gt, class Tds >
inline
Oriented_side
Delaunay_triangulation_sphere_2<Gt,Tds>::
power_test(const Point &p, const Point &q, const Point &r) const
{
	if(number_of_vertices()==2)
		if(orientation_1(p,q)==COLLINEAR)
			return ON_POSITIVE_SIDE;
	return geom_traits().power_test_2_object()(p,q,r);
}
//----------------------------------------------------------------------CHECK---------------------------------------------------------------//
//checks whether a given triangulation is plane (all points are coplanar)
template <class Gt, class Tds>
bool
Delaunay_triangulation_sphere_2<Gt, Tds>::
is_plane()const{
  bool plane = true;
 if(dimension()==2)
	return false;
	
 if(number_of_vertices() > 3){
	All_vertices_iterator it1 = vertices_begin(),
	it2(it1), it3(it1), it4(it1);
	++it2;
	++it3; ++it3;
	++it4; ++it4; ++it4;
	 while( it4 != vertices_end()) {
		Orientation s =power_test(  it1->point(),it2->point(),it3->point(), it4->point()); 
		plane = plane && s == ON_ORIENTED_BOUNDARY ;
		++it1 ; ++it2; ++it3;++it4;
		 if (!plane)
			 return plane;
	 }return true;
    }
	
  if(number_of_vertices()==3)
		return true;
		
  return plane;
}
	
template < class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
is_valid(bool verbose, int level ) const 
{
	bool result = true;
	
	if ( !this-> _tds.is_valid(verbose,level) ) {
		if (verbose)
			std::cerr << "invalid data structure" << std::endl;
	 CGAL_triangulation_assertion(false);
	 return false;
	}
			
	for(All_faces_iterator fit =all_faces_begin(); fit !=all_faces_end(); ++fit) 
		is_valid_face(fit, verbose, level);
			  
   for(All_vertices_iterator vit = vertices_begin(); vit != vertices_end(); ++vit) 
        is_valid_vertex(vit, verbose, level);


   switch(dimension()) {
   case 0 :
     break;
   case 1:
	CGAL_triangulation_assertion(this->is_plane()); 
		   break;
   case 2 :
    for(All_faces_iterator it=all_faces_begin(); it!=all_faces_end(); it++) {
      Orientation s = orientation(it->vertex(0)->point(),  it->vertex(1)->point(), it->vertex(2)->point());
	  result = result && (s != NEGATIVE || it->is_ghost());
	  CGAL_triangulation_assertion(result);
	 }
   
     result = result && (number_of_faces() == 2*(number_of_vertices()) - 4 );
     CGAL_triangulation_assertion( result);
     break;
   }
  
   // in any dimension
   if(verbose) 
     std::cerr << " number of vertices " << number_of_vertices() << "\t" << std::endl;
   
     CGAL_triangulation_assertion( result);
   return result;
}

	
template < class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
is_valid_vertex(Vertex_handle vh, bool verbose, int level) const
{
	bool result = vh->face()->has_vertex(vh);
	if ( !result ) {
	 if(verbose) {
	   std::cerr << " from is_valid_vertex " << std::endl;
       std::cerr << "normal vertex " << &(*vh) << std::endl;
       std::cerr << vh->point() << " " << std::endl;
       std::cerr << "vh_>face " << &*(vh->face())  << " " << std::endl;
	  show_face(vh->face());
     }
  CGAL_triangulation_assertion(false);
  return false;
	}
return true;
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
is_valid_face(Face_handle fh, bool verbose, int level) const
{
   bool result = fh->get_in_conflict_flag()==0;
  	for (int i = 0; i<+2; i++)    {
	   Orientation test = power_test(fh, fh->vertex(i)->point());
	  result = result && test == ON_ORIENTED_BOUNDARY;
	  CGAL_triangulation_assertion(result); 
	}
	if (!result)   
		if (verbose){
			std::cerr<<" from is_valid_face "<<std::endl;
			std::cerr<<" face "<<std::endl;
			show_face(fh);
		}
    CGAL_triangulation_assertion(result); 
  return result;
}

	
	//tests whether there is a conflict between p and the face fh
template < class Gt, class Tds >
inline bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
test_conflict(const Point  &p, Face_handle fh) const
{
	return(power_test(fh,p, true) != ON_NEGATIVE_SIDE);
}  


//----------------------------------------------------------------------INSERTION-------------------------------------------------------------//

template < class Gt, class Tds >
typename Delaunay_triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt,Tds>::
insert(const Point &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = Base::locate(p, lt, li, start);
  switch (lt){
	case NOT_ON_SPHERE: 		
      return Vertex_handle();
	case TOO_CLOSE:{
		 if(dimension()==2)  
	        return loc->vertex(li);
		 return Vertex_handle();
	  }
	case VERTEX:{
	 if(number_of_vertices()==1)
		 return vertices_begin();
	 return (loc->vertex(li));
   }
	default: //point can be inserted
     return insert(p, lt, loc, li);
  }
}

//inserts a new point to a 1d triangulation, the new point is also coplanar with the existing points.	
template < class Gt, class Tds>
typename Delaunay_triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt, Tds>::
insert_cocircular(const Point &p, Locate_type lt, Face_handle loc)
{
	
	CGAL_triangulation_precondition(!test_dim_up(p));
	CGAL_triangulation_precondition(dimension()==1);
	
			
		Vertex_handle v0 = loc->vertex(0);
		Vertex_handle v1 = loc->vertex(1);
		Vertex_handle v=this->_tds.create_vertex();
		v->set_point(p);
		
		Face_handle f1 = this->_tds.create_face(v0, v, Vertex_handle());
		Face_handle f2 = this->_tds.create_face(v, v1, Vertex_handle());
		
		v->set_face(f1);
		v0->set_face(f1);
		v1->set_face(f2);
		
		this->_tds.set_adjacency(f1,0,f2,1);
		this->_tds.set_adjacency(f1,1,loc->neighbor(1),0);
		this->_tds.set_adjacency(f2,0,loc->neighbor(0),1);
		
		
		this->delete_face(loc);
		
	   update_ghost_faces(v);
		return v;
}
	

template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt,Tds>::
insert_first(const Point& p)
{
	CGAL_triangulation_precondition(number_of_vertices() == 0);
	Vertex_handle v =this->_tds.insert_first();
	v->set_point(p);
	return  v;
}
	
template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt,Tds>::
insert_second(const Point& p)
{
	CGAL_triangulation_precondition(number_of_vertices() == 1);
	Vertex_handle v =this->_tds.insert_second();
	v->set_point(p);
	return v;
}
	

//inserts a point which location is known. Calls the corresponding insert-function
	//e.g. insert_first
template < class Gt, class Tds >
typename Delaunay_triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt,Tds>::
insert(const Point &p, Locate_type lt, Face_handle loc, int li)
{
 Vertex_handle v;
 switch (dimension()){
	case -2 :
		return insert_first(p);
	case -1:
		return insert_second(p);
	case 0:
		return insert_outside_affine_hull_regular(p);
	case 1:{
  	  if(test_dim_up(p)){
		Face_handle f=all_edges_begin()->first;
		Vertex_handle v1=f->vertex(0);
		Vertex_handle v2=f->vertex(1);
		Vertex_handle v3=f->neighbor(0)->vertex(1);
		return insert_outside_affine_hull_regular(p);
       } 
	   else 
		return insert_cocircular(p,lt,loc);	
	}
	case 2:{
		std::vector<Face_handle> faces;
		std::vector<Edge> edges;
		faces.reserve(32);
		edges.reserve(32);
						
		get_conflicts_and_boundary(p, std::back_inserter(faces), std::back_inserter(edges), loc);
		v =this->_tds.star_hole(edges.begin(), edges.end());
		v->set_point(p);
		delete_faces(faces.begin(), faces.end());
							
		if( lt != FACE )
			update_ghost_faces(v);
												
		return v;		
	      }
	}
	CGAL_assertion(false);
	return v;
}	

	
//inserts a new point which lies outside the affine hull of the other points
template <class Gt, class Tds >
typename Triangulation_sphere_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_sphere_2<Gt,Tds>::
insert_outside_affine_hull_regular(const Point& p)
{
 if(dimension()==0){
	Vertex_handle v = vertices_begin();			
	Vertex_handle u=v->face()->neighbor(0)->vertex(0);
	Vertex_handle nv;
			
	//orientation is given by the 2 first points
	if( this->collinear_between(v->point(),u->point(),p) || orientation(u->point(),v->point(),p) == LEFT_TURN ) 
		nv=Base::tds().insert_dim_up(v,false);
	else 
		nv=Base::tds().insert_dim_up(v,true); 
	
	nv->set_point(p);
		
	Face_handle f=all_edges_begin()->first; 
	CGAL_triangulation_assertion( orientation(f->vertex(0)->point(), f->vertex(1)->point(),
											   f->neighbor(0)->vertex(1)->point()) != RIGHT_TURN );
			
	 update_ghost_faces(nv);
	 return nv;
   }
		
 else{ 	//dimension 1	
	bool conform = false;
	Face_handle f = (all_edges_begin())->first;
	Face_handle fn=f->neighbor(0);
	const Point p0=f->vertex(0)->point();
	const Point p1=f->vertex(1)->point();
	const Point p2=fn->vertex(1)->point();
	
	 CGAL_triangulation_assertion(orientation(p0, p1, p2)!=NEGATIVE);
	Orientation orient2 = power_test(p0, p1, p2, p);
	
	 
	if(orient2==POSITIVE)
	  conform =true;
				
	//find smalest vertex this step garanties a unique triangulation
	Vertex_handle w=vertices_begin();
	 All_vertices_iterator vi;
	for( vi = vertices_begin(); vi != vertices_end(); vi++){
		if(compare_xyz(vi->point(), w->point())==SMALLER){
				w=vi;
			}
	}
				
	Vertex_handle v = this->_tds.insert_dim_up( w, conform);
	v->set_point(p);
	this->_ghost=all_faces_begin();
	update_ghost_faces(v, true);
			
	return v; 
  }
}
	

	
/*method to mark faces incident to v as ghost-faces or solid-faces.
 bool first defines whether dimension of the triangulation increased
 from one to two by inserting v
 If first == true all faces are updated. 
 */
template < class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
update_ghost_faces(Vertex_handle v, bool first)
{
	if (number_of_vertices()<3)
		return false;
	
  bool neg_found=false;
  if(dimension()==1){
	 All_edges_iterator eit=all_edges_begin();
	 for (; eit!=all_edges_end(); eit++){
        Face_handle f=eit->first;
        Face_handle fn=f->neighbor(0);
        Point q=fn->vertex(1)->point();
		 if(this->collinear_between(f->vertex(0)->point(),f->vertex(1)->point(),q)){
	        f->ghost()=true;
		    neg_found = true;
		 } 
          else 
	        f->ghost()=false;
        }
	}
  else{//dimension==2
	if (first){ //first time dimension 2
	   All_faces_iterator fi = all_faces_begin();
		for(;fi!=all_faces_end();fi++){
			if(orientation(fi)!=POSITIVE){
				fi->ghost()=true;
				neg_found=true;
				this->_ghost=fi;
			  }
			  else
				  fi->ghost()=false;
		 }
	}
			  
	else { //not first
	  Face_circulator fc=this->incident_faces(v,v->face());
	  Face_circulator done(fc);
       do{
		if(orientation(fc)!=POSITIVE){
	      fc->ghost()=true;
	      neg_found=true;
	     this->_ghost=fc;
        }
        else
	      fc->ghost()=false;
	   }while(++fc!=done);
	 }
    
   }
	return neg_found;

}

//-------------------------------------------------------------------------------REMOVAL----------------------------------------------------//
template < class Gt, class Tds >
void
Delaunay_triangulation_sphere_2<Gt,Tds>::
remove_degree_3(Vertex_handle v, Face_handle f) 
{
  if (f == Face_handle())    f=v->face();
  this->_tds.remove_degree_3(v,f);
}

	
	
	
template < class Gt, class Tds >
void
Delaunay_triangulation_sphere_2<Gt,Tds>::
remove(Vertex_handle v )
{
CGAL_triangulation_precondition( v != Vertex_handle() );
	if(number_of_vertices()<=3)
		this->_tds.remove_dim_down(v);
	
	else if(dimension()==2) 
		remove_2D (v);
	
	else 
		remove_1D(v);
	
}

	
	
template <class Gt, class Tds >
void
Delaunay_triangulation_sphere_2<Gt, Tds>::
remove_1D(Vertex_handle v)
{
	this->_tds.remove_1D(v);
	update_ghost_faces();
}
	
	
	
	
template < class Gt, class Tds >
void
Delaunay_triangulation_sphere_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
	CGAL_triangulation_precondition(dimension()==2);
	
  if (test_dim_down(v)) { //resulting triangulation has dim 1
    this->_tds.remove_dim_down(v);
	update_ghost_faces(); //1d triangulation, no vertex needed to update ghost-faces
  }
  else {
    std::list<Edge> hole;
    this->make_hole(v, hole);
    fill_hole_regular(hole);
    }
  return;   
}

	
	
//tests whether the dimension will decrease when removing v from the triangulation.	
template <class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
test_dim_down(Vertex_handle v)
{
	CGAL_triangulation_precondition(dimension()==2);
	bool  dim1 = true; 
  if(number_of_vertices()==4){
    return dim1;
  }
 
	std::vector<Point> points;
	All_vertices_iterator it = vertices_begin();
	for(; it!=vertices_end();it ++)
		if (it != v)
			points.push_back(it->point());
	
	for(int i=0; i<(int)points.size()-4; i++){
		Orientation s = power_test(points.at(i), points.at(i+1),points.at(i+2),points.at(i+3));
		dim1 = dim1 && s == ON_ORIENTED_BOUNDARY ;
		if (!dim1)
			return dim1;
	}return true;
	
}

	
//tests whether the dimension will increase when adding p to the triangulation.	
template <class Gt, class Tds >
bool
Delaunay_triangulation_sphere_2<Gt,Tds>::
test_dim_up(const Point &p) const
	{ 
	CGAL_triangulation_precondition(dimension()!=2);
	Face_handle f=all_edges_begin()->first;
	Vertex_handle v1=f->vertex(0);
	Vertex_handle v2=f->vertex(1);
	Vertex_handle v3=f->neighbor(0)->vertex(1);
		
	v1->point();
	return (power_test(v1->point(), v2->point(), v3->point(),p)!=ON_ORIENTED_BOUNDARY);
		
}
	
//fill the hole in a triangulation after vertex removal.	
template < class Gt, class Tds >
void
Delaunay_triangulation_sphere_2<Gt,Tds>::
fill_hole_regular(std::list<Edge> & first_hole)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;
  
  Hole hole;
  Hole_list hole_list;
  Face_handle ff, fn;
  int i, ii, in;
	
  hole_list.push_front(first_hole);
  
  while (! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      typename Hole::iterator hit = hole.begin();
	    
      // if the hole has only three edges, create the triangle
      if (hole.size() == 3)
	{
	  Face_handle  newf = create_face();
	  hit = hole.begin();
	  for(int j=0; j<3; j++)
	    {
	      ff = (*hit).first;
	      ii = (*hit).second;
	      hit++;
	      ff->set_neighbor(ii,newf);
	      newf->set_neighbor(j,ff);
	      newf->set_vertex(newf->ccw(j),ff->vertex(ff->cw(ii)));
	    }
	  if(orientation(newf) != POSITIVE){
	    this->_ghost=newf;
	    newf->ghost()=true;
	  }
	  continue;
	}
  
      // else find an edge with two finite vertices
      // on the hole boundary
      // and the new triangle adjacent to that edge
      //  cut the hole and push it back
 
      // take the first neighboring face and pop it;
      ff = hole.front().first;
      ii = hole.front().second;
      hole.pop_front();
 
      Vertex_handle  v0 = ff->vertex(ff->cw(ii)); 
      const Point& p0 = v0->point();
      Vertex_handle  v1 = ff->vertex(ff->ccw(ii)); 
      const Point& p1 = v1->point();
      Vertex_handle  v2 = Vertex_handle(); 
      Point p2;
      Vertex_handle  vv;
      Point p;
 
      typename Hole::iterator hdone = hole.end();
      hit = hole.begin();
      typename Hole::iterator cut_after(hit);
 
      // if tested vertex is c with respect to the vertex opposite
      // to NULL neighbor,
      // stop at the before last face;
      hdone--;
      while (hit != hdone) 
	{
	  fn = (*hit).first;
	  in = (*hit).second;
	  vv = fn->vertex(ccw(in));
	  p = vv->point();
	  if ( /*orientation(p0,p1,p) == COUNTERCLOCKWISE*/ 1)
	    {
	      if (v2==Vertex_handle())
	      {
	         v2=vv;
	         p2=p;
	         cut_after=hit;
              }
	      else if (power_test(p0,p1,p2,p) == 
		       ON_POSITIVE_SIDE)
              {
	         v2=vv;
	         p2=p;
	         cut_after=hit;
	      }
            }
	    
	  ++hit;
	}
 
      // create new triangle and update adjacency relations
      Face_handle newf = create_face(v0,v1,v2);
      newf->set_neighbor(2,ff);
      ff->set_neighbor(ii, newf);
      if(orientation(newf) != POSITIVE){
	this->_ghost=newf;
	newf->ghost()=true;
      }
      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes
 
      fn = hole.front().first;
      in = hole.front().second;
      if (fn->has_vertex(v2, i) && i == (int)fn->ccw(in)) 
	{
	  newf->set_neighbor(0,fn);
	  fn->set_neighbor(in,newf);
	  hole.pop_front();
	  hole.push_front(Edge(newf,1));
	  hole_list.push_front(hole);
	}
      else
	{
	  fn = hole.back().first;
	  in = hole.back().second;
	  if (fn->has_vertex(v2, i) && i == (int)fn->cw(in)) 
	    {
	      newf->set_neighbor(1,fn);
	      fn->set_neighbor(in,newf);
	      hole.pop_back();
	      hole.push_back(Edge(newf,0));
	      hole_list.push_front(hole);
	    }
	  else
	    { // split the hole in two holes
	      Hole new_hole;
	      ++cut_after;
	      while (hole.begin() != cut_after)
		{
		  new_hole.push_back(hole.front());
		  hole.pop_front();
		}
 
	      hole.push_front(Edge(newf,1));
	      new_hole.push_front(Edge(newf,0));
	      hole_list.push_front(hole);
	      hole_list.push_front(new_hole);
	    }
	}
    }
}

//-----------------dual------------------------	
	
	//Following methods are used to compute the Voronoi-Diagram
template<class Gt, class Tds>
inline
typename Delaunay_triangulation_sphere_2<Gt, Tds>::	Point
Delaunay_triangulation_sphere_2<Gt, Tds>::
dual (Face_handle f)const 
{
	CGAL_triangulation_precondition(this->_tds.is_face(f));
	CGAL_triangulation_precondition (this->dimension()==2);
	return circumcenter(f);
}
	
template < class Gt, class Tds >
Object
Delaunay_triangulation_sphere_2<Gt,Tds>::
dual(const Edge &e) const
{
	CGAL_triangulation_precondition(this->_tds.is_edge(e.first,e.second));
	CGAL_triangulation_precondition(dimension()==2);
 	 Segment s = this->geom_traits().construct_segment_2_object()
	 (dual(e.first),dual(e.first->neighbor(e.second)));
	 return make_object(s);
}
	
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_sphere_2<Gt,Tds>::  
dual(const Edge_circulator& ec) const
{
	return dual(*ec);
}
	
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_sphere_2<Gt,Tds>::
dual(const All_edges_iterator& ei) const
{
	return dual(*ei);
}
	
} //namespace CGAL 

#endif // CGAL_Delaunay_triangulation_sphere_2_H
