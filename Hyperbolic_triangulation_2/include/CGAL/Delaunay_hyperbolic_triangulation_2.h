// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
#define CGAL_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <stack>
#include <set>

namespace CGAL {
    
class Hyperbolic_face_info_2
{
public:
  Hyperbolic_face_info_2() : _is_finite_invisible(false), _invisible_edge(UCHAR_MAX)
  {
  }
  
  bool is_finite_invisible() const
  {
    return _is_finite_invisible;
  }
  
  void set_finite_invisible(bool is_finite_invisible)
  {
    _is_finite_invisible = is_finite_invisible;
  }
  
  // Supposed to be called before "get_invisible_edge"
  bool has_invisible_edge() const
  {
    return _invisible_edge <= 2;
  }
  
  // Higly recommended to call "has_invisible_edge" before 
  unsigned char get_invisible_edge() const
  {
    assert(_is_finite_invisible);
    assert(_invisible_edge <= 2);
    
    return _invisible_edge;
  }
  
  void set_invisible_edge(unsigned char invisible_edge)
  {
    assert(_is_finite_invisible);
    assert(invisible_edge <= 2); 
    
    _invisible_edge = invisible_edge;
  }
  
private:
  // a face is invisible if its circumscribing circle intersects the circle at infinity
  bool _is_finite_invisible;
  
  // defined only if the face is finite and invisible
  unsigned char _invisible_edge;
};
    
template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                           Triangulation_vertex_base_2<Gt>, 
                           Triangulation_face_base_with_info_2<Hyperbolic_face_info_2, Gt> > >
class Delaunay_hyperbolic_triangulation_2 : public Delaunay_triangulation_2<Gt,Tds>
{
public:
  typedef Delaunay_hyperbolic_triangulation_2<Gt, Tds> Self;
  typedef Delaunay_triangulation_2<Gt,Tds> Base;
  
  typedef Triangulation_face_base_with_info_2<Hyperbolic_face_info_2, Gt> Face_base;
  typedef typename Face_base::Info Face_info;
  
  typedef typename Base::size_type             size_type;
  
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Face_handle   Face_handle;
  typedef typename Base::Edge          Edge;
  
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
  using Base::cw;
  using Base::ccw;
  using Base::geom_traits;
#endif
  
  typedef typename Base::Edge_circulator       Edge_circulator;
  typedef typename Base::Face_circulator       Face_circulator;
  typedef typename Base::Vertex_circulator     Vertex_circulator;
  
  typedef typename Base::All_vertices_iterator    All_vertices_iterator;
  typedef typename Base::All_edges_iterator       All_edges_iterator;
  typedef typename Base::All_faces_iterator       All_faces_iterator;
 
  typedef Gt Geom_traits;
  typedef typename Geom_traits::FT            FT;
  typedef typename Geom_traits::Point_2       Point;
  typedef typename Geom_traits::Segment_2     Segment;
  
  typedef typename Geom_traits::Triangulation_euclidean_traits_2 Euclidean_geom_traits;
  
  /*
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::side_of_oriented_circle;
  using Triangulation::circumcenter;
  using Triangulation::collinear_between;
  using Triangulation::test_dim_down;
  using Triangulation::make_hole;
  using Triangulation::fill_hole_delaunay;
  using Triangulation::delete_vertex;
#endif
*/
  
  Delaunay_hyperbolic_triangulation_2(const Gt& gt = Gt())
  : Delaunay_triangulation_2<Gt,Tds>(gt) {}
  
  Delaunay_hyperbolic_triangulation_2(
	     const Delaunay_hyperbolic_triangulation_2<Gt,Tds> &tr)
       : Delaunay_triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition( this->is_valid() );}

  
  void mark_star(Vertex_handle v) const
  {
    if(!is_star_bounded(v)) {
      mark_star_faces(v);
    }
  }
  
  Vertex_handle insert(const Point  &p, 
                       Face_handle start = Face_handle() )
  {
    Vertex_handle v = Base::insert(p, start);
    mark_star(v);
    
    return v;
  }
  
  Vertex_handle insert(const Point& p,
                       typename Base::Locate_type lt,
                       Face_handle loc, int li )
  {
    Vertex_handle v = Base::insert(p, lt, loc, li);
    mark_star(v);    
    
    return v;
  }
    
#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first, InputIterator last,
         typename boost::enable_if<
         boost::is_base_of<
         Point,
         typename std::iterator_traits<InputIterator>::value_type
         >
         >::type* = NULL
         )
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO 
  {
    size_type n = Base::insert(first, last);
    
    mark_faces();
    
    return n;
  }
  
  //test version of insert function
  
#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert2( InputIterator first, InputIterator last,
         typename boost::enable_if<
         boost::is_base_of<
         Point,
         typename std::iterator_traits<InputIterator>::value_type
         >
         >::type* = NULL
         )
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert_2(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO 
  {
    size_type n = this->number_of_vertices();
    
    spatial_sort(first, last, geom_traits());
    Face_handle f;
    while(first != last) {
      f = insert (*first++, f)->face();
    }
  
    return this->number_of_vertices() - n;
  }
  
  bool is_infinite(Vertex_handle v) const
  {
    return Base::is_infinite(v);
  }
  
  bool is_infinite(Face_handle f) const
  {
    return has_infinite_vertex(f) || is_finite_invisible(f);
  }
  
  bool is_infinite(Face_handle f, int i) const 
  {
    return has_infinite_vertex(f, i) || is_finite_invisible(f, i);
  }
  
  bool is_infinite(const Edge& e) const 
  {
    return is_infinite(e.first, e.second);
  }
  
  bool is_infinite(const Edge_circulator& ec) const 
  {
    return is_infinite(*ec);
  }
  
  bool is_infinite(const All_edges_iterator& ei) const 
  {
    return is_infinite(*ei);
  }
  
private:
  
  bool has_infinite_vertex(Face_handle f) const
  {
    return Base::is_infinite(f);
  }
  
  bool has_infinite_vertex(Face_handle f, int i) const
  {
    return Base::is_infinite(f, i);
  }
  
  bool has_infinite_vertex(const Edge& e) const
  {
    return Base::is_infinite(e);
  }
  
  int get_finite_invisible_edge(Face_handle f) const
  {
    assert(is_finite_invisible(f));
    
    return f->info().get_invisible_edge(); 
  }
  
  bool is_finite_invisible(Face_handle f) const
  {
    return f->info().is_finite_invisible();
  }
  
  bool is_finite_invisible(Face_handle f, int i) const
  {
    if(this->dimension() <= 1) {
      return false;
    }
    
    if(is_finite_invisible(f) && get_finite_invisible_edge(f) == i) {
      return true;
    }
    
    // another incident face and corresponding index
    Face_handle f2 = f->neighbor(i);
    int i2 =  f2->index(f);
    
    if(is_finite_invisible(f2) && get_finite_invisible_edge(f2) == i2) {
      return true;
    }
    
    return false;
  }
  
  bool is_finite_invisible(const Edge& e) const
  {
    return is_finite_invisible(e.first, e.second);
  }
  
  // Depth-first search (dfs) and marking the finite invisible faces.
  void mark_faces() const
  {
    if(this->dimension() <= 1) return;
      
    std::set<Face_handle> visited_faces;
    
    // maintain a stack to be able to backtrack
    // to the most recent faces which neighbors are not visited
    std::stack<Face_handle> backtrack;
    
    // start from a face with infinite vertex
    Face_handle current = Base::infinite_face();
    
    // mark it as visited
    visited_faces.insert(current);
    
    // put the element whose neighbors we are going to explore.
    backtrack.push(current);
    
    // test whether a face is finite invisible or not
    Mark_face test(*this);
    
    Face_handle next;
    Face_info face_info;
    
    while(!backtrack.empty()) {
      // take a face
      current = backtrack.top();
      
      // start visiting the neighbors
      int i = 0;
      for(; i < 3; i++) {
        next = current->neighbor(i);
        
        // if a neighbor is already visited, then stop going deeper
        if(visited_faces.find(next) != visited_faces.end()) {
          continue;
        }
        
        visited_faces.insert(next);
        mark_face(next, test);
        
        // go deeper if the neighbor is infinite
        if(is_infinite(next)) {
          backtrack.push(next);
          break;
        }
      }
      
      // if all the neighbors are already visited, then remove "current" face.
      if(i == 3) {
        backtrack.pop();
      }
    }
    
  }
  
  // check if a star is bounded by finite faces
  // TODO: rename this function name
  bool is_star_bounded(Vertex_handle v) const
  {
    if(this->dimension() <= 1) {
      return true;
    }
    
    Face_handle f = v->face();
    Face_handle next;
    int i;
    Face_handle start(f);
    
    Face_handle opposite_face;
    
    do {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      
      opposite_face = f->neighbor(i);
      if(this->is_infinite(opposite_face)) {
        return false;
      }
      
      f = next;
    } while(next != start);
    
    return true;
  }
  
  //use the function: insert_and_give_new_faces?
  
  void mark_star_faces(Vertex_handle v) const
  {
    // TODO: think of it
    if(this->dimension() <= 1) return;
    
    Mark_face test(*this);
    
    Face_handle f = v->face();
    Face_handle start(f), next;
    int i;
    do {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      
      mark_face(f, test);
      
      f = next;
    } while(next != start);
    return;
  }
  
  template<class Mark_face_test>
  void mark_face(const Face_handle& f, const Mark_face_test& test) const
  {
    f->info() = test(f);
  }
  
  void mark_face(const Face_handle& f) const
  {
    Mark_face test(*this);
    mark_face(f, test);
  }
  
  class Is_hyperbolic
  {
  public:
    Is_hyperbolic()
    {
    }
    
    bool operator() (const Face_handle& f) const
    {
      typedef typename Gt::Vector_3 Vector_3;
      
      Point p0 = f->vertex(0)->point();
      Point p1 = f->vertex(1)->point();
      Point p2 = f->vertex(2)->point();
      
      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(), 
                                                 p1.x()*p1.x() + p1.y()*p1.y(), 
                                                 p2.x()*p2.x() + p2.y()*p2.y());
      
      Vector_3 v1 = Vector_3(p0.x(), p1.x(), p2.x());
      Vector_3 v2 = Vector_3(p0.y(), p1.y(), p2.y());
      Vector_3 v3 = Vector_3(FT(1), FT(1), FT(1));
      
      FT dt0 = CGAL::determinant(v0, v1, v3);
      FT dt1 = CGAL::determinant(v0, v2, v3);
      FT dt2 = CGAL::determinant(v0 - v3, v1, v2);
      
      return dt0*dt0 + dt1*dt1 - dt2*dt2 < 0;
    }
    
    // assume f is "non-hyperbolic"
    bool operator() (const Face_handle& f, int i) const
    {
      typedef typename Gt::Vector_2 Vector_2;
      typedef typename Gt::Vector_3 Vector_3;
      typedef typename Gt::Point_3 Point_3;
      
      Point p0 = f->vertex(cw(i))->point();
      Point p1 = f->vertex(ccw(i))->point();
      Point q0 = f->vertex(i)->point();
      
      FT p0p0 = p0.x()*p0.x() + p0.y()*p0.y();
      FT p1p1 = p1.x()*p1.x() + p1.y()*p1.y();
      FT q0q0 = q0.x()*q0.x() + q0.y()*q0.y();
      
      Vector_3 v0 = Vector_3(2*p0.x(), 2*p1.x(), 2*q0.x());
      Vector_3 v1 = Vector_3(2*p0.y(), 2*p1.y(), 2*q0.y());
      Vector_3 v2 = Vector_3(FT(-1), FT(-1), FT(-1));
      Vector_3 v = Vector_3(p0p0, p1p1, q0q0);
      
      FT x = CGAL::determinant(v, v1, v2);
      FT y = CGAL::determinant(v0, v, v2);
      FT z = CGAL::determinant(v0, v1, v);
      
      FT t = CGAL::determinant(v0, v1, v2);
      
      Point_3 A = Point_3(x/t, y/t, z/t);
      
      Vector_3 norm_p0(2*p0.x(), 2*p0.y(), FT(-1));
      Vector_3 norm_p1(2*p1.x(), 2*p1.y(), FT(-1));
      Vector_3 norm_q0(2*q0.x(), 2*q0.y(), FT(-1));
      
      Vector_3 main_direction = CGAL::cross_product(norm_p0, norm_p1);
      
      if(main_direction * norm_q0 > 0) {
        main_direction = -main_direction;
      }
      //std::cout << "direction " << main_direction << std::endl;
      
      Vector_3 temp = Vector_3(-q0.x(), -q0.y(), 1 - q0q0);
      assert(temp * norm_q0 < 0);
      
      FT k = (1 + A.z())/(2*main_direction.z());
      if (k < 0) {
        k = -k;
      }
      
      Point_3 B = Point_3();
      if(A.z() + main_direction.z() < -1) {
        B = A + k*main_direction;
      } else {
        B = A + main_direction;
      }
      
      Point phiA = Point((2.0*A.x())/(1.0 + A.z()), (2.0*A.y())/(1.0 + A.z()));
      Point phiB = Point((2.0*B.x())/(1.0 + B.z()), (2.0*B.y())/(1.0 + B.z()));
      
      Vector_2 phiAB = Vector_2(phiA, phiB);
      Vector_2 phiOA = Vector_2(Point(0, 0), phiA);
      
      return phiAB*phiOA < 0;
    }
    
    // assume the incident faces are "non-hyperbolic"
    bool operator() (const Edge& e) const
    {
      typedef typename Gt::Vector_2 Vector_2;
      typedef typename Gt::Vector_3 Vector_3;
      
      // endpoints of the edge
      Point p0 = e.first->vertex(cw(e.second))->point();
      Point p1 = e.first->vertex(ccw(e.second))->point();
      
      Point q0 = e.first->vertex(e.second)->point();
      // neighbor of e.first
      Face_handle f = e.first->neighbor(e.second);
      int ind = f->index(e.first);
      Point q1 = f->vertex(ind);
      
      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(), 
                             p1.x()*p1.x() + p1.y()*p1.y(), 
                             q0.x()*q0.x() + q0.y()*q0.y());
      
      Vector_3 v1 = Vector_3(2*p0.y(), 2*p1.y(), 2*q0.y());
      Vector_3 v2 = Vector_3(2*p0.x(), 2*p1.x(), 2*q0.x());
      Vector_3 v3 = Vector_3(FT(-1), FT(-1), FT(-1));
      
      FT dt0 = CGAL::determinant(v0, v1, v3);
      FT dt1 = CGAL::determinant(v2, v0, v3);
      
      FT m11 = p0.x() * (p1.x()*p1.x() + p1.y()*p1.y() - 2) - p1.x() * (p0.x()*p0.x() + p0.y()*p0.y() - 2); 
      FT m12 = p0.y() * (p1.x()*p1.x() + p1.y()*p1.y() - 2) - p1.y() * (p0.x()*p0.x() + p0.y()*p0.y() - 2); 
      
      Vector_3 v4 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(), 
                             p1.x()*p1.x() + p1.y()*p1.y(), 
                             q1.x()*q1.x() + q1.y()*q1.y());
      
      Vector_3 v5 = Vector_3(2*p0.y(), 2*p1.y(), 2*q1.y());
      Vector_3 v6 = Vector_3(2*p0.x(), 2*p1.x(), 2*q1.x());
      Vector_3 v7 = Vector_3(FT(-1), FT(-1), FT(-1));
      
      FT dt3 = CGAL::determinant(v4, v5, v7);
      FT dt4 = CGAL::determinant(v6, v4, v7);
      
      FT big_dt0 = CGAL::determinant(Vector_2(dt0, m11), Vector_2(dt1, m12));
      FT big_dt1 = CGAL::determinant(Vector_2(dt3, m11), Vector_2(dt4, m12));
      
      assert(big_dt0 == 0 || big_dt1 == 0);
      return (big_dt0 > 0 && big_dt1 > 0) || (big_dt0 < 0 && big_dt1 < 0);
    }
    
  private:
    Is_hyperbolic(const Is_hyperbolic&);
    Is_hyperbolic& operator= (const Is_hyperbolic&);
  };
  
  class Mark_face
  {
  public:
    typedef typename Gt::Circle_2 Circle_2;
    
    Mark_face(const Self& tr) :
      _tr(tr)
    {}
    
    Face_info operator ()(const Face_handle& f) const
    {
      Face_info info;
      if(_tr.has_infinite_vertex(f)) {
        return info;
      }
      
      if(Is_hyperbolic()(f) == false) {
        
        info.set_finite_invisible(true);
        info.set_invisible_edge(find_invisible_edge_2(f));
        
        return info;
      }
      
      return info;
    }
    
  private:
    
    // assume that the circumscribing circle of a given face intersects
    // the circle at infinity
    unsigned char find_invisible_edge(const Face_handle& f) const
    {
      typedef Euclidean_geom_traits Egt;
      typedef typename Egt::Construct_circumcenter_2 Construct_circumcenter_2;
      typedef typename Egt::Direction_2 Direction_2;
      
      assert(!_tr.has_infinite_vertex(f));
      
      Point p0 = f->vertex(0)->point();
      Point p1 = f->vertex(1)->point();
      Point p2 = f->vertex(2)->point();
      
      _circle = Circle_2(p0, p1, p2);
      Point c = _circle.center();
      
      Direction_2 d0(p0.x()-c.x(), p0.y()-c.y());
      Direction_2 d1(p1.x()-c.x(), p1.y()-c.y());
      Direction_2 d2(p2.x()-c.x(), p2.y()-c.y());
      
      Direction_2 d(c.x(), c.y());
      
      if(d.counterclockwise_in_between(d0, d1)) {
        return 2;
      }
      
      if(d.counterclockwise_in_between(d1, d2)) {
        return 0;
      }
      
      return 1;
    }
    
    unsigned char find_invisible_edge_2(const Face_handle& f) const
    {
      
      if( Is_hyperbolic()(f, 0) == false ) {
        return 0;
      }
      if( Is_hyperbolic()(f, 1) == false ) {
        return 1;
      };
      return 2;
    }
    
    Mark_face(const Mark_face&);
    Mark_face& operator= (const Mark_face&);
    
    mutable Circle_2 _circle;
    
    const Self& _tr;
  }; 
  
public:
  // This class is used to generate the Finite_*_iterators.
  class Infinite_hyperbolic_tester
  {
    const Self *t;
  public:
    Infinite_hyperbolic_tester() {}
    Infinite_hyperbolic_tester(const Self *tr)	  : t(tr) {}
    
    bool operator()(const All_vertices_iterator & vit) const  {
      return t->is_infinite(vit);
    }
    bool operator()(const All_faces_iterator & fit) const {
      return t->is_infinite(fit);
    }
    bool operator()(const All_edges_iterator & eit ) const {
      return t->is_infinite(eit);
    }
  };
  
  Infinite_hyperbolic_tester
  infinite_hyperbolic_tester() const
  {
    return Infinite_hyperbolic_tester(this);
  }
  
  //Finite faces iterator
  
  class Finite_faces_iterator
  : public Filter_iterator<All_faces_iterator, Infinite_hyperbolic_tester> 
  {
    typedef Filter_iterator<All_faces_iterator, Infinite_hyperbolic_tester> Base;
    typedef Finite_faces_iterator                           Self;
  public:
    Finite_faces_iterator() : Base() {}
    Finite_faces_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator const Face_handle() const { return Base::base(); }
  };
  
  Finite_faces_iterator
  finite_faces_begin() const
  {
    if ( this->dimension() < 2 )
      return finite_faces_end();
    return CGAL::filter_iterator(this->all_faces_end(),
                                 Infinite_hyperbolic_tester(this),
                                 this->all_faces_begin() );
  } 
  
  Finite_faces_iterator
  finite_faces_end() const
  {
    return CGAL::filter_iterator(this->all_faces_end(),
                                 Infinite_hyperbolic_tester(this)   );
  }
  
  //Finite edges iterator
  
  typedef Filter_iterator<All_edges_iterator, 
  Infinite_hyperbolic_tester>
  Finite_edges_iterator;
  
  Finite_edges_iterator
  finite_edges_begin() const
  {
    if ( this->dimension() < 1 )
      return finite_edges_end();
    return CGAL::filter_iterator(this->all_edges_end(),
                                 infinite_hyperbolic_tester(),
                                 this->all_edges_begin());
  }
  
  Finite_edges_iterator
  finite_edges_end() const
  {
    return CGAL::filter_iterator(this->all_edges_end(),
                                 infinite_hyperbolic_tester() );
  }
  
  using Base::dual;
  
  Object
  dual(const Finite_edges_iterator& ei) const
  {
    return this->dual(*ei);
  }
  
  Object
  dual(const Edge &e) const
  {
    CGAL_triangulation_precondition (!this->is_infinite(e));
    
    if(this->dimension() == 1) {
      const Point& p = (e.first)->vertex(cw(e.second))->point();
      const Point& q = (e.first)->vertex(ccw(e.second))->point();
      
      // hyperbolic line
      Segment line = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
      return make_object(line);
    }
    
    // incident faces
    Face_handle f1 = e.first;
    int i1 = e.second;
    
    Face_handle f2 = f1->neighbor(i1);
    int i2 = f2->index(f1);
    
    // boths faces are infinite, but the incident edge is finite
    if(is_infinite(f1) && is_infinite(f2)){
      const Point& p = (f1)->vertex(cw(i1))->point();
      const Point& q = (f1)->vertex(ccw(i1))->point();
      
      // hyperbolic line
      Segment line = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
      return make_object(line);
    }
    
    // both faces are finite
    if(!is_infinite(f1) && !is_infinite(f2)) {
      
      Segment s = this->geom_traits().construct_segment_2_object()
      (dual(f1),dual(f2));
      
      return make_object(s);
    }
    
    // one of the incident faces is infinite
    Face_handle finite_face = f1;
    int i = i1;
    
    if(is_infinite(f1)) {
      finite_face = f2;
      i = i2;
    }
    
    const Point& p = finite_face->vertex(cw(i))->point();
    const Point& q = finite_face->vertex(ccw(i))->point();
    
    // ToDo: Line or Segment?
    // hyperbolic line and ray
    Segment line = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
    Segment ray = this->geom_traits().construct_ray_2_object()(dual(finite_face), line);
    return make_object(ray);
  }
};
  
} //namespace CGAL

#endif // CGAL_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
