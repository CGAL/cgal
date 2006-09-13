#ifndef CGAL_SLICE_ARRANGEMENT_H
#define CGAL_SLICE_ARRANGEMENT_H

//#include "types.h"
#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Algebraic_kernel_2_2.h>
#include <CGAL/Circular_kernel.h>
#include <CGAL/Exact_circular_kernel.h>
#include <CGAL/Arr_circular_line_arc_traits.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <CGAL/Cartesian.h>
#include <sstream>
//#include <CGAL/Arr_landmarks_point_location.h>

class Qt_examiner_viewer_2;

class Slice_arrangement {
public:
  typedef Slice_data_structure::Curve Curve;
  typedef Slice_data_structure::Point Point;

  typedef CGAL::Cartesian<CGAL::Gmpq> K;
  //typedef K::FT                                          NT;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<CGAL::Gmpq>     Algebraic_k;
  typedef CGAL::Circular_kernel_2<K, Algebraic_k>        Circular_k;
  //typedef CGAL::Exact_circular_kernel_2                  Circular_k;
  typedef Circular_k::FT                                 NT;
  typedef Circular_k::Root_of_2                          Exact_NT;
  //typedef Circular_k::Root_for_circles_2_2               RFC;
  typedef Circular_k::Circular_arc_2                     Circular_arc_2;
  typedef Circular_k::Line_arc_2                         Line_arc_2;
  typedef Circular_k::Point_2                            Point_2;
  typedef Circular_k::Vector_2                           Vector_2;
  typedef Circular_k::Line_2                             Line_2;
  typedef Circular_k::Circle_2                           Circle_2;
  typedef Circular_k::Point_3                            Point_3;
  typedef Circular_k::Vector_3                           Vector_3;
  typedef Circular_k::Line_3                             Line_3;
  typedef Circular_k::Sphere_3                           Sphere_3;

  typedef CGAL::Arr_circular_line_arc_traits
  <Circular_k, Circular_arc_2, Line_arc_2>  Arr_traits;
  /*typedef CGAL::Arr_extended_dcel<Arr_traits,
    int, bool, int>      Dcel;*/
  typedef CGAL::Arrangement_with_history_2<Arr_traits>                     CArr;
  typedef CGAL::Arr_naive_point_location<CArr>             Point_location;
  
  
  template <class Point_3, class NT>
  inline  NT squared_depth(Point_3 p, const NT &z){
    return CGAL::square(NT(p[sweep_coordinate().index()])-z);
  }
  
  template <class Sphere, class NT>
  bool has_overlap(const Sphere &s, const NT &z) {
    NT dz2= squared_depth(s.center(), z);
    if (dz2 <= NT(s.squared_radius())) return true;
    else return false;
  }
  
  template <class Sphere, class NT>
  Circular_k::Circle_2 intersect(Sphere s, const NT &z){
    NT r2= NT(s.squared_radius()) - squared_depth(s.center(), z);
    CGAL_assertion(r2>=0);
    Circular_k::Circle_2  c(Circular_k::Point_2(s.center()[plane_coordinate(0).index()],
						s.center()[plane_coordinate(1).index()]), r2);
    return c;
  }

public:

  template <class It>
  Slice_arrangement(It b, It e, NT z, NT inf): pl_(arr_),  inf_(inf){
    std::vector<Circular_k::Circle_2> circles;
    std::vector<int> names;
    int num=0;
    nums_= std::distance(b,e);
    for (; b != e; ++b){
      if (has_overlap(*b, z)){
	//std::cout << circles.size() << " is " << i << std::endl;
	circles.push_back(intersect(*b, z));
	names.push_back(num);
      }
      ++num;
    }
    build_arrangement(circles, names);
   
  }
  

  
  std::size_t number_of_faces() const;
  std::size_t number_of_vertices() const;
  std::size_t number_of_halfedges() const;

  class Face_iterator;
  
  Face_iterator faces_begin() const;
  Face_iterator faces_end() const;

 

  
  void draw(Qt_examiner_viewer_2 *qtv) const;

protected:
  struct Get_segment;
  struct Get_circle;

  template <int C> class Rule;

  void build_arrangement(const std::vector<Circular_k::Circle_2> &circles, 
			 const std::vector<int> &names);
  typedef CArr::Vertex_iterator Vertex_iterator;
  Vertex_iterator vertices_begin();
  Vertex_iterator vertices_end();
  typedef boost::variant< Circular_arc_2, Line_arc_2>       Arc;
  void insert(const Circular_k::Line_arc_2 &la, Curve f);
  void insert(const Circular_k::Circular_arc_2 &la, Curve f);
  void audit() const;
  Curve curve(CArr::Halfedge_const_handle) const;
  friend class Face_iterator;
  //typedef Point Vertex;
  Point point(CArr::Vertex_const_handle h) const;
  //typedef CArr::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  struct Line_arc_less {
    bool operator()(Circular_k::Line_arc_2 a, 
		    Circular_k::Line_arc_2 b) const;
  };
  struct Circular_arc_less {
    bool operator()(Circular_k::Circular_arc_2 a, 
		    Circular_k::Circular_arc_2 b) const;
  };
  
  void write(std::ostream& out, const Circular_k::Circular_arc_2 &k) const;
  void write(std::ostream& out, const Circular_k::Line_arc_2 &k) const;

  /*static Sphere_3 unproject(Circle_2 c);
  static Line_3 unproject(Line_2 c);
  static Point_3 unproject(Point_2 c);*/

  struct Curve_handle_less {
    bool operator()(CArr::Curve_const_handle a, CArr::Curve_const_handle b) const {
      CGAL_precondition(a != CArr::Curve_const_handle());
      CGAL_precondition(b != CArr::Curve_const_handle());
      return &*a < &*b;
    }
  };
  struct Vertex_handle_less {
    bool operator()(CArr::Vertex_const_handle a, CArr::Vertex_const_handle b) const {
      CGAL_precondition(a != CArr::Vertex_const_handle());
      CGAL_precondition(b != CArr::Vertex_const_handle());
      return &*a < &*b;
    }
  };

  /*struct Curve_lookup: boost::static_visitor<>{
    Curve_lookup(Arr *m): a_(m){}
    typedef Curve result_type;
    result_type operator()(const Circular_k::Circular_arc_2 &k) const {
    if (a_->circle_map_.find(k) == a_->circle_map_.end()) {
    a_->write(std::cout, k);
    }
    CGAL_precondition(a_->circle_map_.find(k) != a_->circle_map_.end());
    return a_->circle_map_.find(k)->second;
    }
    result_type operator()(const Circular_k::Line_arc_2 &k) const {
    if (a_->line_map_.find(k) == a_->line_map_.end()) {
    a_->write(std::cout, k);
    }
    CGAL_precondition(a_->line_map_.find(k) != a_->line_map_.end());
    return a_->line_map_.find(k)->second;
    }
    Arr *a_;
    };
    friend struct Curve_lookup;*/

  struct VHandle: public CArr::Vertex_const_handle {
    typedef CArr::Vertex_const_handle P;
    VHandle(CArr::Vertex_const_handle v): CArr::Vertex_const_handle(v){
      //std::cout << v->point() << std::endl;
    }
    bool operator<( VHandle b) const {
      //CGAL_precondition(a != CArr::Vertex_const_handle());
      //CGAL_precondition(b != CArr::Vertex_const_handle());
      return & P::operator*() < &b.operator*();
    }
  private:
    VHandle(){}
  };
  
  CArr arr_;
  Point_location pl_;
  //std::map<Circular_k::Line_arc_2, Curve, Line_arc_less> line_map_;
  //std::map<Circular_k::Circular_arc_2, Curve, Circular_arc_less> circle_map_;
  std::map<CArr::Curve_handle, Curve, Curve_handle_less> map_;
  //mutable std::set<Point> vset_;
  mutable std::map<VHandle, Point> vmap_;
  NT inf_;
  int nums_;
  //Curve_lookup lookup_curve_data_;
};


class Slice_arrangement::Face_iterator {
public:
  Face_iterator(CArr::Face_const_iterator c, 
		CArr::Face_const_iterator e,
		const Slice_arrangement *a);
  
  typedef std::pair<Curve,  Point> FP;
  
  typedef std::vector<std::pair<Curve,  Point> > value_type;
  typedef const value_type& reference_type;
  typedef const value_type* pointer_type;
  
  bool operator!=(const Face_iterator &o) const {
    return c_!= o.c_;
  }

  Face_iterator operator++() {
    increment();
    return *this;
  }
  Face_iterator operator++(int i) {
    Face_iterator o=*this;
    increment();
    return o;
  }

  reference_type operator*() const{
    return cache_;
  }
  pointer_type operator->() const{
    return &cache_;
  }
private:
  void fill_cache();
  void increment();
  value_type cache_;
  CArr::Face_const_iterator c_, e_;
  const Slice_arrangement* a_;
};
#endif
