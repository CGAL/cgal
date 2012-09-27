#ifndef CGAL_Regular_triangulation_on_sphere_2_H
#define CGAL_Regular_triangulation_on_sphere_2_H

#define HOLE_APPROACH

#include <CGAL/Triangulation_on_sphere_2.h>
#include <CGAL/Regular_triangulation_face_base_on_sphere_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/utility.h>
#include <fstream>


namespace CGAL { 

//TODO :
//power_test(const Face_handle& f, int i,const Point &p)
//copy constructor and assignation
//insertion, location, removal, flips and hidding dim 0 and 1
//Optimize Hole Approach and write it clearfully
//And more...
//test orientation in secure conflict for negative faces


template < class Gt,
           class Tds  = Triangulation_data_structure_2 <
                        Regular_triangulation_vertex_base_2<Gt>,
		        Regular_triangulation_face_base_on_sphere_2<Gt> > >
class Regular_triangulation_on_sphere_2 
  : public Triangulation_on_sphere_2<Gt,Tds>
{
  typedef Regular_triangulation_on_sphere_2<Gt, Tds>                         Self;
  typedef Triangulation_on_sphere_2<Gt,Tds>                                  Base;
public:
  typedef Tds                                  Triangulation_data_structure;
  typedef Gt                                   Geom_traits;
  typedef typename Gt::Point_2                 Point;

  typedef typename Base::size_type             size_type;
  typedef typename Base::Face_handle           Face_handle;
  typedef typename Base::Vertex_handle         Vertex_handle;
  typedef typename Base::Vertex                Vertex;
  typedef typename Base::Edge                  Edge;
  typedef typename Base::Locate_type           Locate_type;
  typedef typename Base::Face_circulator       Face_circulator;
  typedef typename Base::Edge_circulator       Edge_circulator;
  typedef typename Base::Vertex_circulator     Vertex_circulator;
  typedef typename Base::Edges_iterator        Edges_iterator;
  typedef typename Base::Faces_iterator        Faces_iterator;
  typedef typename Base::Face::Vertex_list     Vertex_list;
  typedef typename Vertex_list::iterator       Vertex_list_iterator;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;
  using Base::dimension;
  using Base::full_sphere;
  using Base::geom_traits;
  using Base::create_face;
  using Base::number_of_faces;
  using Base::faces_begin;
  using Base::faces_end;
  using Base::edges_begin;
  using Base::edges_end;
  using Base::OUTSIDE_AFFINE_HULL;
  using Base::VERTEX;
  using Base::FACE;
  using Base::EDGE;
  using Base::OUTSIDE_CONVEX_HULL;
  using Base::orientation;
  using Base::insert_four_init_vertices;
#endif

private:

  class Hidden_tester {
  public:
    bool operator()(const typename Base::Vertices_iterator&  it){
      return it->is_hidden();
     }
  };

  class Unhidden_tester {
  public:
    bool operator()(const typename Base::Vertices_iterator&  it){
      return ! it->is_hidden();
    }
  };

  typedef typename Base::Vertices_iterator     vertices_ib;

public:
  // We derive in order to add a conversion to handle.
  class Vertices_iterator :
    public Filter_iterator<vertices_ib, Hidden_tester> {
    typedef Filter_iterator<vertices_ib, Hidden_tester> Base;
    typedef Vertices_iterator                     Self;
     public:
    Vertices_iterator() : Base() {}
    Vertices_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator Vertex_handle() const { return Base::base(); } 
  };


  class Hidden_vertices_iterator :
    public Filter_iterator<vertices_ib, Unhidden_tester> {
    typedef Filter_iterator<vertices_ib, Unhidden_tester> Base; 
    typedef Hidden_vertices_iterator                     Self;
  public:
    Hidden_vertices_iterator() : Base() {}
    Hidden_vertices_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator Vertex_handle() const { return Base::base(); }
 };


 private:
  typedef std::list<Face_handle>      Faces_around_stack; 
  size_type _hidden_vertices;
   

  //PREDICATES

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


 public: 
  //CONSTRUCTORS
  Regular_triangulation_on_sphere_2(const Gt& gt=Gt()) 
    : Base(Gt(gt)), _hidden_vertices(0)
  {}

  void clear();
  //CHECK
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid_face(Face_handle fh) const;
  bool is_valid_vertex(Vertex_handle fh) const;

  void show_face(Face_handle fh) const;

  void show_all() const;

  bool test_conflict(const Point  &p, Face_handle fh) const;
  bool test_conflict_1(const Point  &p, Face_handle fh) const;

  bool update_negative_faces(Vertex_handle v=Vertex_handle());

  bool check_neighboring()
  {
    Faces_iterator eit;
    if(dimension()==1){
      for(eit=faces_begin();eit!=faces_end();++eit)
      {
	  Face_handle f1 = eit->neighbor(0);
	  Face_handle f2 = eit->neighbor(1);
	  CGAL_triangulation_assertion(f1->has_neighbor(eit));
	  CGAL_triangulation_assertion(f2->has_neighbor(eit));
      }
      
    }

    for(eit=faces_begin();eit!=faces_end();++eit)
      {
	  Face_handle f1 = eit->neighbor(0);
	  Face_handle f2 = eit->neighbor(1);
	  Face_handle f3 = eit->neighbor(2);
	  CGAL_triangulation_assertion(f1->has_neighbor(eit));
	  CGAL_triangulation_assertion(f2->has_neighbor(eit));
	  CGAL_triangulation_assertion(f3->has_neighbor(eit));
      }

  }

  //ACCES FUNCTION
  size_type number_of_vertices() const {
    return Base::number_of_vertices() - _hidden_vertices;
  }
 
  size_type number_of_hidden_vertices() const {
    return _hidden_vertices;
  }
  Face_handle create_star_2(const Vertex_handle& v, const Face_handle& c, int li);

  //INSERTION
  Vertex_handle insert(const Point &p, Locate_type  lt, Face_handle loc, int li );
  Vertex_handle insert_in_face(const Point &p, Face_handle f);
  Vertex_handle insert_in_edge(const Point &p, Face_handle f, int i);
  Vertex_handle insert(const Point &p, Face_handle f = Face_handle() );
  Vertex_handle reinsert(Vertex_handle v, Face_handle start);
  Vertex_handle insert_outside_affine_hull_regular(const Point& p,bool plane=false);

  Vertex_handle insert_hole_approach_2(const Point &p, Locate_type lt, Face_handle loc, int li) ;
  //REMOVAL
  void remove_degree_3(Vertex_handle v, Face_handle f = Face_handle());
  void remove(Vertex_handle v);
  void remove_hidden(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  bool test_dim_down(Vertex_handle v);
  void fill_hole_regular(std::list<Edge> & hole);

  //FLIP
  void flip(Face_handle f, int i);

  //HIDDING
  Vertex_handle hide_new_vertex(Face_handle f, const Point& p);
  void hide_vertex(Face_handle f, Vertex_handle v);
  void update_hidden_points_2_2(const Face_handle& f1, const Face_handle& f2);
  void update_hidden_points_3_1(const Face_handle& f1, const Face_handle& f2, const Face_handle& f3);
  void update_hidden_points_1_3(const Face_handle& f1, const Face_handle& f2,const Face_handle& f3);
  void hide_remove_degree_3(Face_handle fh, Vertex_handle vh);
  void exchange_incidences(Vertex_handle va, Vertex_handle vb);
  void set_face(Vertex_list& vl, const Face_handle& fh);
  void exchange_hidden(Vertex_handle va, Vertex_handle vb);



  //REGULAR 
  void regularize(Vertex_handle v);
  void stack_flip(Vertex_handle v, Faces_around_stack &faces_around);
  void stack_flip_4_2(Face_handle f, int i, int j, 
		      Faces_around_stack &faces_around);
  void stack_flip_3_1(Face_handle f, int i, int j,
		      Faces_around_stack &faces_around);
  void stack_flip_2_2(Face_handle f, int i, 
		      Faces_around_stack &faces_around);
  void stack_flip_dim1(Face_handle f, int i,
		       Faces_around_stack &faces_around);

  //ITERATORS
  Vertices_iterator vertices_begin () const;
  Vertices_iterator vertices_end () const;

  Hidden_vertices_iterator hidden_vertices_begin () const;
  Hidden_vertices_iterator hidden_vertices_end () const;

  //---------------------------------------------------------------------HOLE APPROACH



  //TEMPLATE MEMBERS
  //----------------------------------------------------------------------HOLE APPROACH
    template <class InputIterator, class OutputIterator>
    void process_faces_in_conflict(InputIterator start, InputIterator end, OutputIterator vertices) const
    {
      int dim = dimension();
      while (start != end) {
	for( typename Vertex_list::iterator it = (*start)->vertex_list().begin();
	     it != (*start)->vertex_list().end();
	     ++it){
	  (*it)->set_face(Face_handle());
	  *vertices++ = (*it);
	}
	(*start)->vertex_list().clear();
	
	for (int i=0; i<=dim; i++) {
	  Vertex_handle v = (*start)->vertex(i);
	  if (v->face() != Face_handle()) {
	    *vertices++= v;
	    v->set_face(Face_handle());
	  }
	}
	start ++;
      }
    }

	template <class InputIterator>
    void reinsert_vertices(Vertex_handle v, InputIterator start, InputIterator end ) {
	  Face_handle hf = v->face();
	  while(start != end){
	    if ((*start)->face() == Face_handle()){
	    hf = locate ((*start)->point(), hf);
	    hide_vertex(hf, *start);
	    }
	    start++;
	  }
	}

template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
      int n = number_of_vertices();

      std::vector<Point> points (first, last);
      std::random_shuffle (points.begin(), points.end());
      //~ spatial_sort (points.begin(), points.end(), geom_traits());
      spatial_sort (points.begin(), points.end());

      Face_handle hint;
      for (typename std::vector<Point>::const_iterator p = points.begin(),
		      end = points.end();
              p != end; ++p)
          hint = insert (*p, hint)->face();

      return number_of_vertices() - n;
  }

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
  Stream &write_triangulation_to_off_2(Stream &out,Stream &out2){

    // Points of triangulation
    for (Faces_iterator it = this->_tds.face_iterator_base_begin(); it != faces_end(); it++) {  
    if(!it->is_negative()
       /*(t.orientation(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point())==1)*/
       ){//assert(orientation(it)==POSITIVE);
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
Stream &write_triangulation_to_off(Stream &out) {
  std::vector<Face_handle> faces;

  // Points of triangulation
  for (Faces_iterator it = faces_begin(); it != faces_end(); it++) {
    for (int i=0 ; i<3 ; i++) {
      Point p = it->vertex(i)->point();
      out << p.x() << " " 
	  << p.y() << " " 
	  << p.z() << std::endl;
    }
  }

  return out;
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

template <class Stream >
Stream &write_face_to_off(Stream &out,Face_handle f){

      for (int i=0 ; i<3 ; i++) {
	Point p = f->vertex(i)->point();
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
      }
  

  return out;
}


template <class Stream,class FaceIt>
Stream &write_faces_to_off(Stream &out,FaceIt face_begin, FaceIt face_end){

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
Stream &write_edges_to_off(Stream &out,FaceIt face_begin, FaceIt face_end){

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



  template <
            class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces,
            class OutputIteratorInternalEdges>
  Triple<OutputIteratorBoundaryEdges,
         OutputIteratorFaces,
         OutputIteratorInternalEdges>
  find_conflicts(Point p,Face_handle f, 
		 Triple<OutputIteratorBoundaryEdges,
                        OutputIteratorFaces,
		        OutputIteratorInternalEdges> it) const
  {
 
    std::stack<Face_handle> Face_stack;
    Face_stack.push(f);
    f->set_in_conflict_flag(1);
    *it.second++ = f;

    do {
        Face_handle c = Face_stack.top();
		Face_stack.pop();
		
        for (int i=0; i<3; ++i) {
          Face_handle test = c->neighbor(i);
          if (test->get_in_conflict_flag() == 1) {
	    continue; // test was already in conflict.
          }

          if (test->get_in_conflict_flag() == 0) {
	    if (test_conflict(p,test)) {
		  Face_stack.push(test);
		  test->set_in_conflict_flag(1);
		  *it.second++ = test;
	          continue;
	        }
     	    test->set_in_conflict_flag(2); // test is on the boundary.
          }
          *it.first++ = Edge(c, i);
        }
    } while(!Face_stack.empty());
    return it;
  }
  
};
//------------------------------------------------------------------PREDICATES---------------------------------------------------------------------//
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
clear()
{
  Base::clear();
  _hidden_vertices = 0;
}

template < class Gt, class Tds >
Oriented_side
Regular_triangulation_on_sphere_2<Gt,Tds>::
power_test(const Face_handle &f, const Point &p) const
{
  return power_test(f->vertex(0)->point(),
		    f->vertex(1)->point(),
		    f->vertex(2)->point(),p);
}

template < class Gt, class Tds >
Oriented_side
Regular_triangulation_on_sphere_2<Gt,Tds>::
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
Regular_triangulation_on_sphere_2<Gt,Tds>::
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
Regular_triangulation_on_sphere_2<Gt,Tds>::
power_test(const Point &p,
	   const Point &q,
	   const Point &r) const
{
  if(number_of_vertices()==2)
    if(orientation_1(p,q)==COLLINEAR)
      return ON_POSITIVE_SIDE;
  return geom_traits().power_test_2_object()(p,q,r);
}

template < class Gt, class Tds >
inline
Oriented_side
Regular_triangulation_on_sphere_2<Gt,Tds>::
power_test(const Point &p,
	   const Point &r) const
{
  return geom_traits().power_test_2_object()(p,r);
}

//----------------------------------------------------------------------CHECK---------------------------------------------------------------//
template < class Gt, class Tds >
bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
  is_valid(bool verbose, int ) const //int level
{
  // cannot call for is_valid() of Base Triangulation class
  // because 1) number of vertices of base class does not match
  // tds.is_valid calls is_valid for each vertex
  // and the test is not fullfilled by  hidden vertices ...
  // result = result && Triangulation_on_sphere_2<Gt,Tds>::is_valid(verbose, level);
  bool result = true;
  for(Faces_iterator fit = faces_begin(); 
      fit != faces_end(); ++fit) {
    result = result && is_valid_face(fit);
  }CGAL_triangulation_assertion(result);
  
  for(Vertices_iterator vit = vertices_begin(); 
                            vit != vertices_end(); ++vit) {
    result = result && is_valid_vertex(vit);
  }CGAL_triangulation_assertion(result);

   for(Hidden_vertices_iterator hvit = hidden_vertices_begin(); 
                                hvit != hidden_vertices_end(); ++hvit) {
    result = result && is_valid_vertex(hvit);
  }CGAL_triangulation_assertion(result);

   switch(dimension()) {
   case 0 :
     break;
   case 1:
      if (number_of_vertices() > 3 ) {
       Vertices_iterator it1 = vertices_begin(),
	 it2(it1), it3(it1);
       ++it2;
       ++it3; ++it3;
       while( it3 != vertices_end()) {
	 Orientation s = orientation(it1->point(),
				    it2->point(),
				    it3->point()); 
	 result = result && s == COLLINEAR ;
	 CGAL_triangulation_assertion(result);
	 ++it1 ; ++it2; ++it3;
       }
      }
     break;
   case 2 :
    for(Faces_iterator it=faces_begin(); 
	 it!=faces_end(); it++) {
      Orientation s = orientation(it->vertex(0)->point(),
				  it->vertex(1)->point(),
				  it->vertex(2)->point());
      CGAL_triangulation_assertion( s == LEFT_TURN || it->is_negative());
      result = result && ( s == LEFT_TURN || it->is_negative());

     
    }

 
     // check number of faces. This cannot be done by the Tds
     // which does not know the number of components nor the genus
     result = result && (number_of_faces() == 2*(number_of_vertices()) - 4 );
     CGAL_triangulation_assertion( result);
     break;
   }
  
   // in any dimension
   if(verbose) {
     std::cerr << " nombres de sommets " << number_of_vertices() << "\t"
	       << "nombres de sommets  caches " << number_of_hidden_vertices()
	       << std::endl;
   }
   result = result && ( Base::number_of_vertices() ==
			number_of_vertices() + number_of_hidden_vertices());
   CGAL_triangulation_assertion( result);
   return result;
}

template < class Gt, class Tds >
bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
is_valid_vertex(Vertex_handle vh) const
{
  bool result = true;
  if (vh->is_hidden()) {
    Locate_type lt; 
    int li;
    Face_handle loc  = locate(vh->point(), lt, li, vh->face());
    if (dimension() == 0) {
      result = result && lt == Base::VERTEX;
        result = result && power_test (vh->face()->vertex(0)->point(), vh->point()) <= 0;
    } else {
        result = result && (loc == vh->face() ||
                (lt == Base::VERTEX && vh->face()->has_vertex(loc->vertex(li))) ||
                (lt == Base::EDGE && vh->face() ==
                 loc->neighbor(li)) );
	if(dimension()>1)
	  {result = result && 
	      power_test(vh->face(),vh->point()) == ON_NEGATIVE_SIDE;}
	if ( !result) {
	  std::cerr << " from is_valid_vertex " << std::endl;
	  std::cerr << "sommet cache " << &*(vh) 
		    << "vh_point " <<vh->point() << " " << std::endl;
               std::cerr << "vh_>face " << &*(vh->face())  << " " << std::endl;
               std::cerr <<  "loc      " <<  &*(loc )
        	        << " lt " << lt  << " li " << li << std::endl;
               show_face(vh->face());
               show_face(loc);
	}
    }
  }
  else { // normal vertex
    result = result && vh->face()->has_vertex(vh);
     if ( !result) {
       std::cerr << " from is_valid_vertex " << std::endl;
       std::cerr << "normal vertex " << &(*vh) << std::endl;
       std::cerr << vh->point() << " " << std::endl;
       std::cerr << "vh_>face " << &*(vh->face())  << " " << std::endl;
       show_face(vh->face());
     }
  }
  CGAL_triangulation_assertion(result);
  return result;
}

template < class Gt, class Tds >
bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
is_valid_face(Face_handle fh) const
{
  bool result = true;
  CGAL_triangulation_assertion(result);
  result = fh->get_in_conflict_flag()==0;

  typename Vertex_list::iterator vlit = fh->vertex_list().begin(),
	                       vldone = fh->vertex_list().end();
  for (; vlit != vldone; vlit++)    {
    result = result && power_test(fh, (*vlit)->point()) == ON_NEGATIVE_SIDE;
    result = result && ((*vlit)->face() == fh);
    if (!result)     show_face(fh);
    CGAL_triangulation_assertion(result); 
  }
  return result;
}


template <class Gt, class Tds >
void
 Regular_triangulation_on_sphere_2<Gt, Tds>::
show_face(Face_handle fh) const
{
  Base::show_face(fh);

  typename Vertex_list::iterator current;
  std::cerr << "  +++++>>>    ";
  for (current= fh->vertex_list().begin(); 
       current!= fh->vertex_list().end() ; current++ ) {
    std::cerr <<"[ "<< ((*current)->point()) << " ] ,  ";
  }
  std::cerr <<std::endl;
}


template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
show_all() const
{
  std::cerr<< "AFFICHE TOUTE LA TRIANGULATION :" << std::endl;
  std::cerr << std::endl<<"====> "<< this ;
  std::cerr <<  " dimension " <<  dimension() << std::endl;
  std::cerr << "nb of vertices " << number_of_vertices() 
	    << " nb of hidden vertices " << number_of_hidden_vertices() 
	    <<   std::endl;

  if (dimension() == 0){
    show_face(vertices_begin()->face());
    show_face(vertices_begin()->face()->neighbor(0));
  }

  if(dimension() == 1) {
    std::cerr<<" all edges "<<std::endl; 
    Edges_iterator aeit;
    for(aeit = edges_begin(); aeit != edges_end(); aeit++){
      show_face(aeit->first);
    }
   }
  
  else{ //dimension ==2
    std::cerr<<" faces "<<std::endl;
    Faces_iterator fi;
    for(fi = faces_begin(); fi != faces_end(); fi++) {
      show_face(fi);
    }
   }
  
  if (number_of_vertices()>1) {
    std::cerr << "affichage des sommets de la triangulation reguliere"
	      <<std::endl;
    Vertices_iterator vi;
    for( vi = vertices_begin(); vi != vertices_end(); vi++){
      show_vertex(vi);
      std::cerr << "  / face associee : "
	     <<  &*(vi->face()) << std::endl;
      }
      std::cerr<<std::endl;
  }
  
   std::cerr << "sommets caches "  << std::endl;
   Hidden_vertices_iterator hvi = hidden_vertices_begin();
   for( ; hvi != hidden_vertices_end(); hvi++) {
     show_vertex(hvi);
      std::cerr << "  / face associee : "
	     << &*(hvi->face()) << std::endl;
   }
  return;
}

template < class Gt, class Tds >
inline bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
test_conflict(const Point  &p, Face_handle fh) const
{
    return(power_test(fh,p) != ON_NEGATIVE_SIDE);
}  

template < class Gt, class Tds >
inline bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
test_conflict_1(const Point  &p, Face_handle fh) const
{
  return( !fh->is_negative() && power_test(fh->vertex(0)->point(),fh->vertex(1)->point(),p) == ON_POSITIVE_SIDE );
}

//----------------------------------------------------------------------INSERTION-------------------------------------------------------------//

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
insert(const Point &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(p, lt, li, start);
  #ifdef HOLE_APPROACH
  return insert_hole_approach_2(p, lt, loc, li);
  #else
  return insert(p, lt, loc, li);
  #endif
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
insert(const Point &p, Locate_type lt, Face_handle loc, int li) 
{
    Vertex_handle v;
    switch (lt) {
    case Base::VERTEX:
      {
	CGAL_precondition (dimension() >= -1);

	if (dimension() < 1 ) {
	  // in this case locate() oddly returns loc = NULL and li = 4,
	  // so we work around it.
	  loc = faces_begin();
	  li = 0;
	}

	Vertex_handle vv = loc->vertex(li);
	CGAL::Oriented_side side = power_test (vv->point(), p);
	switch(side) {
	      
	case ON_NEGATIVE_SIDE:
	  return hide_new_vertex (loc, p);
	      
	case ON_POSITIVE_SIDE:
	  v = this->_tds.create_vertex(); 
	  v->set_point(p);
	  exchange_incidences(v,vv);
	  hide_vertex(loc, vv);
	  regularize (v);
	      return v;
	      
	    case ON_ORIENTED_BOUNDARY:
	      return vv;
	      }
      }
    case Base::EDGE:
        {
            CGAL_precondition (dimension() >= 1);
	    /* Oriented_side os = dimension() == 1 ? power_test (loc, li, p) :
	       power_test (loc, p);*/
	    Oriented_side os = power_test(loc,p);

            if (os < 0) {
                return hide_new_vertex (loc, p);
            }
            v = insert_in_edge (p, loc, li);
            break;
        }	
    case Base::FACE:
      {
            CGAL_precondition (dimension() >= 2);
            if (power_test (loc, p) < 0) {
                return hide_new_vertex(loc,p);
            }
            v = insert_in_face (p, loc);
            break;
      }
    case Base::OUTSIDE_AFFINE_HULL:
      {
	if(dimension()==0){
	  Vertex_handle v=vertices_begin();
	  Vertex_handle u=v->face()->neighbor(0)->vertex(0);
	  loc=v->face();

	  if( collinear_between(u->point(),v->point(),p) && power_test(u->point(),v->point(),p) == ON_NEGATIVE_SIDE )
	    return hide_new_vertex (loc, p);

	  v = Base::insert_outside_affine_hull(p);    
	}
	  
      }
    default:
        v = Base::insert (p, lt, loc, li);
    }
    /*
    if (lt == OUTSIDE_AFFINE_HULL) {
        //clear vertex list of infinite faces which have been copied
        for ( All_faces_iterator afi = all_faces_begin();
                afi != all_faces_end(); afi++)
            if (is_infinite (afi))
                afi->vertex_list().clear();
		}*/


    regularize (v);


    return v;
}


template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Face_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
create_star_2(const Vertex_handle& v, const Face_handle& c, int li )
{
  CGAL_triangulation_assertion( dimension() == 2 );
  Face_handle cnew;

  // i1 i2 such that v,i1,i2 positive
  int i1=ccw(li);
  // traversal of the boundary of region in ccw order to create all
  // the new facets
  Face_handle bound = c;
  Vertex_handle v1 = c->vertex(i1);
  int ind = c->neighbor(li)->index(c); // to be able to find the
                                       // first cell that will be created 
  Face_handle cur;
  Face_handle pnew = Face_handle();
  do {
    cur = bound;
    // turn around v2 until we reach the boundary of region
    while ( cur->neighbor(cw(i1))->get_in_conflict_flag() == 1 ) {
      // neighbor in conflict
      cur = cur->neighbor(cw(i1));
      i1 = cur->index( v1 );
    }
    cur->neighbor(cw(i1))->set_in_conflict_flag(0);
    // here cur has an edge on the boundary of region
    cnew = create_face( v, v1, cur->vertex( ccw(i1) ) );
    this->_tds.set_adjacency(cnew, 0, cur->neighbor(cw(i1)),
	                   cur->neighbor(cw(i1))->index(cur));
    cnew->set_neighbor(1, Face_handle());
    cnew->set_neighbor(2, pnew);
    // pnew is null at the first iteration
    v1->set_face(cnew);
    //pnew->set_neighbor( cw(pnew->index(v1)), cnew );
    if (pnew != Face_handle()) { pnew->set_neighbor( 1, cnew );}

    bound = cur;
    i1 = ccw(i1);
    v1 = bound->vertex(i1);
    pnew = cnew;
    //} while ( ( bound != c ) || ( li != cw(i1) ) );
  } while ( v1 != c->vertex(ccw(li)) );
  // missing neighbors between the first and the last created cells
  cur = c->neighbor(li)->neighbor(ind); // first created cell
  this->_tds.set_adjacency(cnew, 1, cur, 2);
  return cnew;
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
insert_hole_approach_2(const Point &p, Locate_type lt, Face_handle loc, int li) 
{
  Vertex_handle v;
    switch (lt) {
    case Base::VERTEX:
      {
	Vertex_handle vv;

	if(dimension() == -1) {
	  vv=vertices_begin();
	  loc=vv->face();
	}else
	  vv = loc->vertex(li);
	
	//bug here because of removing faces_begin_2
	//CGAL_triangulation_assertion(xy_equal(vv->point(), p));

	CGAL::Oriented_side side = power_test (vv->point(), p);

	switch(side) {
	case ON_NEGATIVE_SIDE:	    
	  return hide_new_vertex (loc, p);
	  
	case ON_POSITIVE_SIDE:
	  if(dimension()<1){
	    v = this->_tds.create_vertex(); 
	    v->set_point(p);
	    exchange_incidences(v,vv);
	    hide_vertex(loc, vv);
	    return v;
	  }break;	
	  
	case ON_ORIENTED_BOUNDARY:
	  return vv;
	}
      }break;

    case Base::EDGE:
        {
            CGAL_precondition (dimension() >= 1);
	    Oriented_side os = dimension() == 1 ? power_test (loc, li, p) : power_test (loc, p);

            if (os < 0) {
                return hide_new_vertex (loc, p);
            }
        }break;	
    case Base::FACE:
	{
            CGAL_precondition (dimension() >= 2);

	    if (power_test (loc, p) < 0) {
                return hide_new_vertex(loc,p);
            }
	}break;
    case Base::OUTSIDE_AFFINE_HULL:
      {
	if(dimension()==0){
	  loc=vertices_begin()->face();
	  v=loc->vertex(0);
	  Vertex_handle u=loc->neighbor(0)->vertex(0);

	  bool between = collinear_between(u->point(),v->point(),p);
	  Oriented_side over = power_test(v->point(),u->point(),p);
	  //p under the edge
	  if( between  &&  over == ON_NEGATIVE_SIDE )
	    return hide_new_vertex (loc, p); 
	  
	  //p over edge and not inside the cone of edge -> hidding one existing vertex
	  if(  !between && over == ON_POSITIVE_SIDE && orientation(v->point(),u->point(),p)==COLLINEAR ) {
	    Vertex_handle nv = this->_tds.create_vertex(); 
	    nv->set_point(p);
	    if( collinear_between(u->point(),p,v->point())){
	      exchange_incidences(nv,v);
	      hide_vertex(nv->face(),v);
	    }else{
	      exchange_incidences(nv,u);
	      hide_vertex(nv->face(),u);
	    }
	    
	    return nv;
	  }
	  
	  return insert_outside_affine_hull_regular(p);

	}else 

	  if(dimension()==1){
	    Face_handle f=edges_begin()->first;
	    Vertex_handle v1=f->vertex(0);
	    Vertex_handle v2=f->vertex(1);
	    Vertex_handle v3=f->neighbor(0)->vertex(1);

	    Orientation orient=orientation(v1->point(),v2->point(),v3->point()) ;

	    if (orient != COLLINEAR 	&&
		power_test (v1->point(),v2->point(),v3->point(),p) == ON_NEGATIVE_SIDE    && 
		oriented_side (v1->point(),v2->point(),v3->point(),p) == ON_POSITIVE_SIDE  
		){
	      return hide_new_vertex(f,p);}
	    else{
	      return insert_outside_affine_hull_regular(p,orient==COLLINEAR);//boolean to know if the triangulation is plane
	    }
	  }       
	  }
    default:
      {
		if(dimension()<1){
			v = Base::insert (p, lt, loc, li);
		}
      }
    }
   if (dimension()==2){
      std::vector<Face_handle> faces;
      std::vector<Vertex_handle> hid_verts;
      Edge edge;
      faces.reserve(32);
      hid_verts.reserve(96);

      find_conflicts
	(p,loc, make_triple(Oneset_iterator<Edge>(edge),
			std::back_inserter(faces),
				Emptyset_iterator()));

      process_faces_in_conflict(faces.begin(),faces.end(),std::back_inserter(hid_verts));
      
      Vertex_handle newv = this->_tds.create_vertex();
      Face_handle fnew;  
      fnew = create_star_2(newv, edge.first, edge.second);
      newv->set_face(fnew);
  
      newv->set_point(p);

      delete_faces(faces.begin(),faces.end());
      if( lt != FACE )
	update_negative_faces(newv);

      //TODO: in the minimal example of five points in
      //test/bugs_minimal_examples.cpp one of the vertices stored in
      //hid_verts has an empty face and therefore is hidden. This is
      //an incorrect behavior, figure out why its face is empty.
      reinsert_vertices(newv,hid_verts.begin(),hid_verts.end());

       return newv;
   }
   //find the conflicting edges
   if( dimension()==1 ){
	std::vector<Face_handle> faces;
	std::vector<Vertex_handle> hid_verts;
	Edge bound[2];
	faces.push_back(loc);

	for (int j = 0; j<2; j++) {
	  Face_handle f = loc->neighbor(j);
	  while ( test_conflict_1(p,f) || (lt==OUTSIDE_CONVEX_HULL && f->is_negative() ) ) {
		faces.push_back(f);
	    f = f->neighbor(j);
	  }
	  Edge e = Edge(f,2);
	  bound[j]=e;
	}
	process_faces_in_conflict(faces.begin(),faces.end(),std::back_inserter(hid_verts));

	delete_faces(faces.begin(),faces.end());
	
		//"star" the hole with v 
	v = this->_tds.create_vertex();
	v->set_point (p);
	Face_handle fb0 = bound[0].first;
	Face_handle fb1 = bound[1].first;
	Face_handle f0 = this->_tds.create_face(v,fb0 ->vertex(0), Vertex_handle());
	Face_handle f1 = this->_tds.create_face(fb1->vertex(1), v, Vertex_handle());
	this->_tds.set_adjacency(fb0,1,f0,0);
	this->_tds.set_adjacency(f0,1,f1,0);
	this->_tds.set_adjacency(f1,1,fb1,0);
	fb0->vertex(0)->set_face(fb0);
	fb1->vertex(1)->set_face(fb1);
	v->set_face(f0);
	v->set_point (p);


	//update the negative face if needed
	if(lt == OUTSIDE_CONVEX_HULL){
 	  Face_handle ff=v->face();
	  int i=ff->index(v);
	  int j=(i+1)%2;
	  Vertex_handle u=ff->vertex(j);
	  Face_handle fn= ff->neighbor(j);
	  Vertex_handle w= fn->vertex(i);
	  if( collinear_between(p,u->point(),w->point()) )
	    ff->negative()=true;
	  else if( collinear_between(p,w->point(),u->point()) )
	    fn->negative()=true;
	}

	reinsert_vertices(v,hid_verts.begin(),hid_verts.end());
	return v;
   }

   return v;
         
}

template < class Gt, class Tds >
bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
update_negative_faces(Vertex_handle v)
{
  if(dimension()==1){
    Edges_iterator eit=edges_begin();
    do{
      Face_handle f=eit->first;
      Face_handle fn=f->neighbor(0);
      Point q=fn->vertex(1)->point();
      if(collinear_between(f->vertex(0)->point(),f->vertex(1)->point(),q))
	 f->negative()=true;
      else 
	f->negative()=false;
      ++eit;
    }while( eit!=edges_end());
  }
  else{//dimension==2

    Face_circulator fc=incident_faces(v,v->face());
    Face_circulator done(fc);
    bool neg_found=false;
	
    do{
      if(orientation(fc)!=POSITIVE){
	fc->negative()=true;
	neg_found=true;
	this->_negative=fc;
      }
      else{
	fc->negative()=false;
      }
    }while(++fc!=done);

    return neg_found;
  }
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
insert_in_face(const Point &p, Face_handle f)
{
  Vertex_handle v = Base::insert_in_face(p,f);
  update_hidden_points_1_3(f, 
			   f->neighbor(ccw(f->index(v))), 
			   f->neighbor( cw(f->index(v))) );
  return v;
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
insert_in_edge(const Point &p, Face_handle f, int i)
{
  Vertex_handle v;
  if (dimension() == 1) {
    v = Base::insert_in_edge(p,f,i);
    Face_handle g = f->neighbor(1 - f->index(v));
    update_hidden_points_2_2(f,g);
  }
  else { //dimension()==2
    // don't use update_hidden_points_2_2 any more to split
    // hidden vertices list because new affectation of f and n
    // around new vertex is unknown
    Face_handle n = f->neighbor(i);
    Vertex_list p_list;
    p_list.splice(p_list.begin(),f->vertex_list());
    p_list.splice(p_list.begin(),n->vertex_list());
    v = Base::insert_in_edge(p,f,i);
    Face_handle loc;
    while ( ! p_list.empty() ){
      loc = locate(p_list.front()->point(), n);
      hide_vertex(loc, p_list.front());
      p_list.pop_front();
    }
  }
  return v;
} 

//The reinsert function  insert a weighted point which was in a hidden
//vertex.
//The new and old vertices are then exchanged ; this is required
//if the regular triangulation is used with a hierarchy because
//the old vertex has its up and down pointers set and other vertices
//pointing on him

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
reinsert(Vertex_handle v, Face_handle start)
{
  CGAL_triangulation_assertion(v->is_hidden());
  v->set_hidden(false);
  _hidden_vertices--;
 
//   //to debug 
//   std::cerr << "from reinsert " << std::endl;
//   show_vertex(v);
//   Locate_type lt;
//   int li;
//   Face_handle loc = locate(v->point(), lt, li, start);
//   std::cerr << "locate " << &(*loc) << "\t" << lt << "\t" << li <<
//     std::endl;
//   show_face(loc);
//    std::cerr << std::endl;

  Vertex_handle vh = insert(v->point(), start);
  if(vh->is_hidden()) exchange_hidden(v,vh);
  else  exchange_incidences(v,vh);
  this->_tds.delete_vertex(vh);
  return v;
}
//-------------------------------------------------------------------------------REMOVAL----------------------------------------------------//
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
remove_degree_3(Vertex_handle v, Face_handle f) 
{
  if (f == Face_handle())    f=v->face();
  update_hidden_points_3_1(f, f->neighbor( cw(f->index(v))),
			   f->neighbor(ccw(f->index(v))));
  Base::remove_degree_3(v,f);
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
remove(Vertex_handle v )
{
    CGAL_triangulation_precondition( v != Vertex_handle() );

    if (v->is_hidden())
        return remove_hidden (v);

    Face_handle hint;
    int ihint = 0;

    Vertex_list to_reinsert;
    switch (dimension()) {
    case -1:
      {
	while(!v->face()->vertex_list().empty()){
	  show_vertex( v->face()->vertex_list().front());
	  remove_hidden( v->face()->vertex_list().front());
	}
      	break;
      }
	
    case 0:
      {
	to_reinsert.splice (to_reinsert.begin(), v->face()->vertex_list());
	to_reinsert.splice (to_reinsert.begin(), v->face()->neighbor(0)->vertex_list());
	break;
      }
    case 1:
        {
	  Face_handle f1 = v->face();
	  ihint = f1->index(v);
	  hint = f1->neighbor(ihint);
	  Face_handle f2 = f1->neighbor(1 - ihint);
	  ihint = mirror_index (f1, ihint);

	  to_reinsert.splice (to_reinsert.begin(), f1->vertex_list());
	  to_reinsert.splice (to_reinsert.begin(), f2->vertex_list());
	  break;
	}
    case 2:
      {
	Face_circulator f = incident_faces (v), end = f;
	ihint = f->index(v);
	hint = f->neighbor(ihint);
	ihint = mirror_index (f, ihint);
	do  to_reinsert.splice (to_reinsert.begin(), f->vertex_list());
	while (++f != end);
	break;
      }
    }
 
    if (number_of_vertices() <= 3) {
       this->_tds.remove_dim_down(v);
    } else if ( dimension() < 2) {
        Base::remove (v);
	} else {
      remove_2D (v);
    }
    if (hint != Face_handle()) hint = hint->neighbor(ihint);
    for (typename Vertex_list::iterator i = to_reinsert.begin();
         i != to_reinsert.end(); ++i) {
      hint = reinsert (*i, hint)->face();
    }
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) { 
    this->_tds.remove_dim_down(v);
    update_negative_faces();
  }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole_regular(hole);
    delete_vertex(v);
  }
  return;   
}

template <class Gt, class Tds >
bool
Regular_triangulation_on_sphere_2<Gt,Tds>::
test_dim_down(Vertex_handle v)
{
  //test the dimensionality of the resulting triangulation
  //upon removing of vertex v
  //it goes down to 1 iff
  // 1) There is only 4 vertices
  //and/or
  // 2) every vertices appart from v are coplanar with _sphere
  CGAL_triangulation_precondition(dimension() == 2);

  bool  dim1 = true; 
  if(number_of_vertices()==4){
    return dim1;
  }
 
  Face_circulator fc=incident_faces(v,v->face());
  Face_circulator done(fc);
  
  do{					     
    int i=fc->index(v);
    Face_handle f=fc->neighbor(i);
    if(orientation(f->vertex(0)->point(),
		   f->vertex(1)->point(),
		   f->vertex(2)->point())
       !=COLLINEAR)
      dim1=false;
  }while(++fc!=done && dim1);
  
  return dim1;
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
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
	    this->_negative=newf;
	    newf->negative()=true;
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
	this->_negative=newf;
	newf->negative()=true;
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

//------------------------------------------------------------------------------FLIP---------------------------------------------------------//

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
flip(Face_handle f, int i)
{
  Face_handle n = f->neighbor(i);
  Base::flip(f,i);
  update_hidden_points_2_2(f,n);
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
remove_hidden(Vertex_handle v )
{
  _hidden_vertices--;
  v->face()->vertex_list().remove(v);
  delete_vertex(v);
  return;
}
//-----------------------------------------------------------------------------HIDDING------------------------------------------------------//
// create a vertex and hide it
template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
hide_new_vertex(Face_handle f, const Point& p)
{
  Vertex_handle v = this->_tds.create_vertex(); 
  v->set_point(p);
  hide_vertex(f, v);
  return v;
}

// insert the vertex to the hidden vertex list
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
hide_vertex(Face_handle f, Vertex_handle vh)
{ 
  if(! vh->is_hidden()) {
    vh->set_hidden(true);
    _hidden_vertices++;
  }
  vh->set_face(f);
  f->vertex_list().push_back(vh);
}

// the points of the lists of 2 faces are sorted
// because of a 2-2 flip
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
update_hidden_points_2_2(const Face_handle& f1, const Face_handle& f2)
{	
  CGAL_triangulation_assertion(f1->has_neighbor(f2));
    
  Vertex_list p_list;
  p_list.splice(p_list.begin(),f1->vertex_list());
  p_list.splice(p_list.begin(),f2->vertex_list());

  /* if (dimension() == 1) {
    const Point& a1 = f1->vertex(f1->index(f2))->point();
    const Point& a  = f1->vertex(1-f1->index(f2))->point();
    while ( ! p_list.empty() ) {
      if ( compare_x(a, p_list.front()->point()) == 
	   compare_x(a, a1)  &&
	   compare_y(a, p_list.front()->point()) == 
	   compare_y(a, a1))
	{
	  hide_vertex(f1, p_list.front());
	} else {
	hide_vertex(f2, p_list.front());
	}
      p_list.pop_front();
    }
    return;
  }*/

  // from here f1 and f2 are 2-dimensional faces
  int idx2 = f1->index(f2);
  Vertex_handle v0=f1->vertex(ccw(idx2));
  Vertex_handle v1=f1->vertex(cw(idx2)); 

  while ( ! p_list.empty() )
    {
      if (orientation(v0->point(), v1->point(), p_list.front()->point()) ==
	  COUNTERCLOCKWISE)
	hide_vertex(f1, p_list.front());
	else
	hide_vertex(f2, p_list.front());
      p_list.pop_front();
    }
}

// add the vertex_list of f2 and f3 to the point list of f1
// for the 3-1 flip
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
update_hidden_points_3_1(const Face_handle& f1, const Face_handle& f2, 
			 const Face_handle& f3)
{
  set_face(f2->vertex_list(), f1);
  set_face(f3->vertex_list(), f1);
  (f1->vertex_list()).splice(f1->vertex_list().begin(),f2->vertex_list());
  (f1->vertex_list()).splice(f1->vertex_list().begin(),f3->vertex_list());
  return;				  
}

// The point list of f1 is separated into 3 lists
// for a 1-3 flip
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
update_hidden_points_1_3(const Face_handle& f1, const Face_handle& f2, 
			 const Face_handle& f3)
{
    CGAL_triangulation_assertion(f1->has_neighbor(f2) &&
			      f2->has_neighbor(f3) &&
			      f3->has_neighbor(f1));


    Vertex_list p_list;
    p_list.splice(p_list.begin(),f1->vertex_list());
    if (p_list.empty())
	return;

    // the following does not work if 
    // two of f1,f2 and f3 are twice neighbors
    // but this cannot appear taking the assertion into account;
    int idx2 = f1->index(f2),
        idx3 = f1->index(f3);
    Vertex_handle v2 = f1->vertex(idx2),
                  v3 = f1->vertex(idx3),
                  v0 = f1->vertex(3-(idx2+idx3)),
                  v1 = f2->vertex(f2->index(f1));

    CGAL_triangulation_assertion(f2->has_vertex(v0) && f1->has_vertex(v0));
    CGAL_triangulation_assertion(f3->has_vertex(v1));

    
    while(! p_list.empty())
    {
      Vertex_handle v(p_list.front());
      if(orientation(v2->point(),v0->point(), v->point()) ==
 	 orientation(v2->point(),v0->point(),v3->point())  &&
	 orientation(v3->point(),v0->point(), v->point()) ==
	 orientation(v3->point(),v0->point(), v2->point()))
	hide_vertex(f1, v);
      else if (orientation(v1->point(), v0->point(), v->point()) ==
	       orientation(v1->point(), v0->point(), v3->point()) )
	hide_vertex(f2,v);
      else hide_vertex(f3,v);
      p_list.pop_front();
    }
}

// the vertex is a degree three vertex which has to be removed
// and hidden
// create first  a new hidden vertex and exchange with the vertex
// to be removed by the tds : 
// this is required to keep up and down pointers right when using a hierarchy
template < class Gt, class Tds >
void 
Regular_triangulation_on_sphere_2<Gt,Tds>::
hide_remove_degree_3(Face_handle fh, Vertex_handle vh)
{
 Vertex_handle vnew= this->_tds.create_vertex();
 exchange_incidences(vnew,vh);
 remove_degree_3(vnew, fh);
 hide_vertex(fh,vh);
}

// set to va the incidences of vb 
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
exchange_incidences(Vertex_handle va, Vertex_handle vb)
{
  CGAL_triangulation_assertion ( !vb->is_hidden());
  std::list<Face_handle> faces;
  if (dimension() < 1) {
    faces.push_back (vb->face());
  } else if (dimension() == 1) {
    faces.push_back(vb->face());
    int i = vb->face()->index(vb);
    faces.push_back(vb->face()->neighbor(1-i));
  }
  else {
    CGAL_triangulation_assertion (dimension() == 2);
    Face_circulator fc = incident_faces(vb), done(fc);
    do {
      faces.push_back(fc);
      fc++;
    }while(fc != done);
  }

  va->set_face(*(faces.begin()));
  for(typename std::list<Face_handle>::iterator it = faces.begin();
      it != faces.end(); it++){
    Face_handle fh = *it;
    fh->set_vertex(fh->index(vb), va);
  }
  return;
}


template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
set_face(Vertex_list& vl, const Face_handle& fh)
{
  for(typename Vertex_list::iterator it = vl.begin(); it != vl.end(); it++)
    (*it)->set_face(fh);
}

//push va instead of vb in the list of the face fb hiding vb
// vb must be the last inserted vertex in the list of fb
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
exchange_hidden(Vertex_handle va, Vertex_handle vb)
{ 
  CGAL_triangulation_assertion (vb->is_hidden());
  CGAL_triangulation_assertion (vb == vb->face()->vertex_list().back());
 
//   //to debug 
//   std::cerr << "from exchange hidden 1" << std::endl;
//   show_vertex(vb);
//   std::cerr << "  / face associee : "
// 	     << &*(vb->face()) << std::endl;
  
  vb->face()->vertex_list().pop_back();
  _hidden_vertices--;
  hide_vertex(vb->face(), va);

//  //to debug 
//   std::cerr << "from exchange hidden 1" << std::endl;
//   show_vertex(va);
//   std::cerr << "  / face associee : "
// 	     << &*(va->face()) << std::endl << std::endl; 
}
//-----------------------------------------------------------------------------REGULAR------------------------------------------------------//
template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
regularize(Vertex_handle v)
{
  Faces_around_stack faces_around;

  if (dimension() < 1) return;

  //initialise faces_around
  if (dimension() == 1) {
    faces_around.push_back(v->face());
    faces_around.push_back(v->face()->neighbor(1- v->face()->index(v)));
  }
  else{ //dimension==2
    Face_circulator fit = incident_faces(v), done(fit);
    do {
      faces_around.push_back(fit++);
    } while(fit != done);
  }

  while( ! faces_around.empty() )
    stack_flip(v, faces_around);
  return;
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
stack_flip(Vertex_handle v, Faces_around_stack &faces_around)
{
  Face_handle f=faces_around.front();
  faces_around.pop_front();
  int i = f->index(v);
  Face_handle n = f->neighbor(i);
 
  if (dimension() == 1 ) {
    if ( f->is_negative()  || power_test( v->point(),
		     n->vertex(n->index(f))->point(),
		     f->vertex(1-i)->point()) ==  ON_NEGATIVE_SIDE)
      stack_flip_dim1(f,i,faces_around);
    return;
  }  

  // now dimension() == 2
  //test the regularity of edge (f,i)
  //if( power_test(n, v->point()) == ON_NEGATIVE_SIDE)
  if( power_test(n, v->point()) != ON_POSITIVE_SIDE)
    return;
    
    int ni = n->index(f);
    Orientation occw = orientation(f->vertex(i)->point(),
				   f->vertex(ccw(i))->point(),
				   n->vertex(ni)->point());
    Orientation ocw  = orientation(f->vertex(i)->point(),
				   f->vertex(cw(i))->point(),
				   n->vertex(ni)->point());
    if (occw == LEFT_TURN && ocw == RIGHT_TURN) {
      // quadrilater (f,n) is convex
      stack_flip_2_2(f,i, faces_around);
      return;
    }
    if (occw == RIGHT_TURN && degree(f->vertex(ccw(i))) == 3) {
      stack_flip_3_1(f,i,ccw(i),faces_around);
      return;
    }
    if (ocw == LEFT_TURN && degree(f->vertex(cw(i))) == 3) {
      stack_flip_3_1(f,i,cw(i),faces_around);
      return;
    }
    if (occw == COLLINEAR && degree(f->vertex(ccw(i))) == 4) {
      stack_flip_4_2(f,i,ccw(i),faces_around);
      return;
    }
    if (ocw == COLLINEAR && degree(f->vertex(cw(i))) == 4)
      stack_flip_4_2(f,i,cw(i),faces_around);
    
    return;
}
template <class Gt, class Tds >
typename Triangulation_on_sphere_2<Gt,Tds>::Vertex_handle
Regular_triangulation_on_sphere_2<Gt,Tds>::
  insert_outside_affine_hull_regular(const Point& p,bool plane)
{

  if(dimension()==0){
    Vertex_handle v=vertices_begin();
    Vertex_handle u=v->face()->neighbor(0)->vertex(0);
    Vertex_handle nv;

    //orientation is given by the 2 first points
    if( collinear_between(v->point(),u->point(),p) || orientation(u->point(),v->point(),p) == LEFT_TURN ) 
      nv=Base::tds().insert_dim_up(v,false);
    else 
      nv=Base::tds().insert_dim_up(v,true);
                                
    nv->set_point(p);

    CGAL_triangulation_assertion( orientation(edges_begin()->first->vertex(0)->point(),
					      edges_begin()->first->vertex(1)->point(),
					      edges_begin()->first->neighbor(0)->vertex(1)->point())
				  != RIGHT_TURN ); 

    //seting negative edge if needed
    bool done=false;
    Edges_iterator eit=edges_begin();
    do{
      Face_handle f=eit->first;
      Face_handle fn=f->neighbor(0);
      Point q=fn->vertex(1)->point();
      if(collinear_between(f->vertex(0)->point(),f->vertex(1)->point(),q)){
	f->negative()=true;
      }
      else{
	f->negative()=false;
      }
      ++eit;
    }while( eit!=edges_end() && !done );
    
    return nv;
  }
  else{ //dimension=1
    bool conform;
    Face_handle f = (edges_begin())->first;
 
    if(plane){//points coplanar with geom_traits->sphere
      Orientation orient = orientation( f->vertex(0)->point(),f->vertex(1)->point(),p);
      conform = ( orient == COUNTERCLOCKWISE);
    }
    else{//three vertices non-coplanar with geom_traits->sphere
      Face_handle fn=f->neighbor(0);

      const Point p0=f->vertex(0)->point();
      const Point p1=f->vertex(1)->point();
      const Point p2=fn->vertex(1)->point();

      Oriented_side side = oriented_side(p0,p1,p2,p);
 
      conform=( side == ON_POSITIVE_SIDE );
    }

    Vertex_handle v = this->_tds.insert_dim_up( f->vertex(0), conform);
    v->set_point(p);
    this->_negative=faces_begin();

    //seting negative faces if needed
    //TODO: Probably there is a bug in the context of neg_found
    bool neg_found=false;
    Face_circulator fc = incident_faces(vertices_begin(),vertices_begin()->face()), done(fc);
    do {
      if(orientation(fc->vertex(0)->point(),fc->vertex(1)->point(),fc->vertex(2)->point())!=POSITIVE){
	fc->negative()=true;
	this->_negative=fc;
      }
      else
	fc->negative()=false;
    } 
    while(++fc != done);

    fc = incident_faces(v,v->face()),  done=fc;
    do {
      if(orientation(fc->vertex(0)->point(),fc->vertex(1)->point(),fc->vertex(2)->point())!=POSITIVE){
	fc->negative()=true;
	neg_found=true;
	this->_negative=fc;
      }
      else
	fc->negative()=false;
    } 
    while(++fc != done);
 
     
    //in some case, 1 existing vertex has to be hidden 
    if(!plane && neg_found){
      int i=this->_negative->index(v);
      Vertex_handle vv=mirror_vertex(this->_negative,i);
      if( power_test(this->_negative,vv->point())==ON_POSITIVE_SIDE ){
	Point q=vv->point();
	remove(vv);
	Face_handle loc = locate (q, fc);
	return hide_new_vertex(loc,q);
      }
    }
    
    return v; 
  }
 
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
stack_flip_4_2(Face_handle f, int i, int j, Faces_around_stack & faces_around)
{
    int k = 3-(i+j);
    Face_handle g=f->neighbor(k);
    if (!faces_around.empty())
    {
      if (faces_around.front() == g)
	  faces_around.pop_front();
      else if (faces_around.back() == g) 
	  faces_around.pop_back();
    }
    
    //union f with  g and f->neihgbor(i) with g->f->neihgbor(i)
    Face_handle fn = f->neighbor(i);
    //Face_handle gn = g->neighbor(g->index(f->vertex(i)));
    Vertex_handle vq = f->vertex(j);
    
    this->_tds.flip( f, i); //not using flip because the vertex j is flat.
    update_hidden_points_2_2(f,fn);
    Face_handle h1 = ( j == ccw(i) ? fn : f);
    //hide_vertex(h1, vq);
    hide_remove_degree_3(g,vq);
    if(j == ccw(i)) {
      faces_around.push_front(h1); 
      faces_around.push_front(g);
    }
    else {
      faces_around.push_front(g);
      faces_around.push_front(h1); 
    }
}


template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
stack_flip_3_1(Face_handle f, int i, int j, Faces_around_stack & faces_around)
{
  int k = 3-(i+j);
  Face_handle g=f->neighbor(k);
  if (!faces_around.empty())
  {
    if (faces_around.front()== g)
	faces_around.pop_front();
    else if ( faces_around.back() == g)
	faces_around.pop_back();
  }

  Vertex_handle vq= f->vertex(j);
  //hide_vertex(f,vq);
  hide_remove_degree_3(f,vq);
  faces_around.push_front(f);
}


template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
stack_flip_2_2(Face_handle f, int i, Faces_around_stack & faces_around)
{
    Vertex_handle vq = f->vertex(ccw(i));
    flip(f,i);
    if(f->has_vertex(vq)) {
      faces_around.push_front(f->neighbor(ccw(i)));
      faces_around.push_front(f);
    }
    else { 
      faces_around.push_front(f);
      faces_around.push_front(f->neighbor(cw(i)));
    }
}

template < class Gt, class Tds >
void
Regular_triangulation_on_sphere_2<Gt,Tds>::
stack_flip_dim1(Face_handle f, int i, Faces_around_stack &faces_around)
{
  Vertex_handle va = f->vertex(1-i);
  Face_handle n= f->neighbor(i);
  int in = n->index(f);
  Vertex_handle vb = n->vertex(in);
  f->set_vertex(1-i, n->vertex(in));
  vb->set_face(f);
  f->set_neighbor(i, n->neighbor(1-in));
  n->neighbor(1-in)->set_neighbor(n->neighbor(1-in)->index(n), f);
  (f->vertex_list()).splice(f->vertex_list().begin(),n->vertex_list());
  set_face(f->vertex_list(),f);
  delete_face(n);
  hide_vertex(f,va);
  faces_around.push_front(f);
  return;
}

 //-----------------------------------------------------------------------------------ITERATORS-----------------------------------------------------------//
template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertices_iterator 
Regular_triangulation_on_sphere_2<Gt,Tds>::
vertices_begin () const
{
  return CGAL::filter_iterator(Base::vertices_end(), 
			 Hidden_tester(),
			 Base::vertices_begin());
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Vertices_iterator 
Regular_triangulation_on_sphere_2<Gt,Tds>::
vertices_end () const
{
  return CGAL::filter_iterator(Base::vertices_end(), 
			 Hidden_tester() ); 
}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Hidden_vertices_iterator 
Regular_triangulation_on_sphere_2<Gt,Tds>::
hidden_vertices_begin () const
{
  return CGAL::filter_iterator(Base::vertices_end(), 
			 Unhidden_tester(), 
			 Base::vertices_begin() );

}

template < class Gt, class Tds >
typename Regular_triangulation_on_sphere_2<Gt,Tds>::Hidden_vertices_iterator 
Regular_triangulation_on_sphere_2<Gt,Tds>::
hidden_vertices_end () const
{
  return CGAL::filter_iterator(Base::vertices_end(), 
			 Unhidden_tester() );
}




//-------------------------------------------------------------------CLASS DEFINITION--------------------------------------------------------//
/*

 //Tag to distinguish Delaunay from Regular triangulations
  typedef Tag_true  Weighted_tag;


public:

  Regular_triangulation_on_sphere_2(const Regular_triangulation_on_sphere_2 &rt);
  
  Regular_triangulation_on_sphere_2 & operator=(const Regular_triangulation_on_sphere_2 &tr);

	
  
   //  //template member functions, declared and defined at the end 
  //  template <class OutputItFaces, class OutputItBoundaryEdges, 
  //                                       class OutputItHiddenVertices> 
  //   Triple<OutputItFaces,OutputItBoundaryEdges, OutputItHiddenVertices>
  //   get_conflicts_and_boundary_and_hidden_vertices (const
  //   Point  &p, 
  // 						  OutputItFaces fit, 
  // 						  OutputItBoundaryEdges eit,
  // 						  OutputItHiddenVertices vit,  
  // 						  Face_handle start = 
  //                                                 Face_handle()) const;
  // template <class OutputItFaces, class OutputItBoundaryEdges> 
  // std::pair<OutputItFaces,OutputItBoundaryEdges>
  // get_conflicts_and_boundary(const Point  &p, 
  // 		                OutputItFaces fit, 
  // 		                OutputItBoundaryEdges eit,
  // 		                Face_handle start) const;
  // template <class OutputItFaces>
  // OutputItFaces
  // get_conflicts (const Point  &p, 
  //                OutputItFaces fit, 
  // 		    Face_handle start ) const;
  // template <class OutputItBoundaryEdges>
  // OutputItBoundaryEdges
  // get_boundary_of_conflicts(const Point  &p, 
  // 			       OutputItBoundaryEdges eit, 
  // 			       Face_handle start ) const;
  //   template <class OutputItBoundaryEdges, class OutputItHiddenVertices> 
  //   std::pair<OutputItBoundaryEdges, OutputItHiddenVertices> 
  //   get_boundary_of_conflicts_and_hidden_vertices(const Point  &p, 
  // 						OutputItBoundaryEdges eit, 
  // 						OutputItHiddenVertices vit,
  // 						Face_handle start=
  //                                                Face_handle()) const;
  //   template <class OutputItHiddenVertices> 
  //   OutputItHiddenVertices
  //   get_hidden_vertices(const Point  &p, 
  // 			   OutputItHiddenVertices vit,
  // 			   Face_handle start= 
  //                       Face_handle()) const;
  
  // DUAL
  Bare_point dual (Face_handle f) const;
  Object dual(const Edge &e) const ;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;
  Bare_point weighted_circumcenter(Face_handle f) const; 
  Bare_point weighted_circumcenter(const Point& p0, 
			      const Point& p1, 
			      const Point& p2) const;

  // Insertion, Deletion and Flip
  Vertex_handle push_back(const Point &p);


  //  Vertex_handle file_input(std::istream& is);
  // void file_output(std::ostream& os) const;

public:

  void copy_triangulation(const Self& tr);
private:
  void copy_triangulation_();




public:


  template < class Stream>
  Stream& draw_dual(Stream & ps) const
    {
      Finite_edges_iterator eit = finite_edges_begin();
      for (; eit != finite_edges_end(); ++eit) {
	Object o = dual(eit);
	typename Geom_traits::Line_2  l;
	typename Geom_traits::Ray_2   r;
	typename Geom_traits::Segment_2 s;
	if (CGAL::assign(s,o)) ps << s;
	if (CGAL::assign(r,o)) ps << r;
	if (CGAL::assign(l,o)) ps << l;
      }
      return ps;
    }

  
  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  get_conflicts_and_boundary (const Point  &p, 
			      OutputItFaces fit, 
			      OutputItBoundaryEdges eit,
			      Face_handle start = Face_handle()) const
    {
      Triple<OutputItFaces,OutputItBoundaryEdges,Emptyset_iterator>
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p, fit, eit,
       						       Emptyset_iterator(), 
       						       start);
      return std::make_pair(pp.first, pp.second);
    }
  template <class OutputItFaces, class OutputItHiddenVertices> 
  std::pair<OutputItFaces, OutputItHiddenVertices> 
  get_conflicts_and_hidden_vertices(const Point  &p, 
				    OutputItFaces fit, 
				    OutputItHiddenVertices vit,
				    Face_handle start = 
				    Face_handle()) const
    {
      Triple<OutputItFaces, Emptyset_iterator,OutputItHiddenVertices> 
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p,fit,
						       Emptyset_iterator(), 
						       vit,
						       start);
      return std::make_pair(pp.first,pp.third);
    }


   template <class OutputItBoundaryEdges, class OutputItHiddenVertices> 
  std::pair<OutputItBoundaryEdges, OutputItHiddenVertices> 
  get_boundary_of_conflicts_and_hidden_vertices(const Point  &p, 
						OutputItBoundaryEdges eit, 
						OutputItHiddenVertices vit,
						Face_handle start = 
						Face_handle()) const
    {
      Triple<Emptyset_iterator,OutputItBoundaryEdges,
	OutputItHiddenVertices> 
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p,
						       Emptyset_iterator(), 
						       eit,vit,
						       start);
      return std::make_pair(pp.second,pp.third);
    }

  template <class OutputItFaces> 
  OutputItFaces
  get_conflicts (const Point  &p, 
		 OutputItFaces fit, 
		 Face_handle start= Face_handle()) const
    {
      Triple<OutputItFaces,Emptyset_iterator,Emptyset_iterator>
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p, fit, 
						       Emptyset_iterator(),
						       Emptyset_iterator(), 
						       start);
      return pp.first;
    }
  
  template <class OutputItBoundaryEdges> 
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Point  &p, 
			    OutputItBoundaryEdges eit, 
			    Face_handle start= Face_handle()) const
    {    
      Triple<Emptyset_iterator, OutputItBoundaryEdges,Emptyset_iterator>
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p,
						       Emptyset_iterator(),
						       eit,
						       Emptyset_iterator(), 
						       start);
      return pp.second;
    }
  template <class OutputItHiddenVertices> 
  OutputItHiddenVertices 
  get_hidden_vertices(const Point  &p, OutputItHiddenVertices vit,
		      Face_handle start= Face_handle()) const
    {
      Triple<Emptyset_iterator,Emptyset_iterator,
	OutputItHiddenVertices> 
	pp = 
	get_conflicts_and_boundary_and_hidden_vertices(p,Emptyset_iterator(), 
						       Emptyset_iterator(),vit,
						       start);
      return pp.third;
    }

  // nearest power vertex query
  Vertex_handle nearest_power_vertex(const Bare_point& p) const;
*/
//---------------------------------------------------------------------END CLASS DEFINITION--------------------------------------------------//






} //namespace CGAL 

#endif // CGAL_Regular_triangulation_on_sphere_2_H
