// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Mariette Yvinec, Jean Daniel Boissonnat
 
#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/triangulation_assertions.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_2/insert_constraints.h>

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/internal/info_check.h>
#include <CGAL/is_iterator.h>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>

namespace CGAL {

namespace internal{

template <class T,bool has_info=is_iterator<T>::value>
struct Get_iterator_value_type{
 struct type{};
};

template <class T>
struct Get_iterator_value_type<T,true>{
 typedef typename std::iterator_traits<T>::value_type type;
};

} } //namespace CGAL::internal
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO



namespace CGAL {


template <class Gt,
          class Tds_ = Default ,
	  class Itag_ = Default>
class Constrained_Delaunay_triangulation_2
  : public  Constrained_triangulation_2<Gt, Tds_, Itag_>
{
public:
  typedef Constrained_triangulation_2<Gt,Tds_,Itag_>            Ctr;
  typedef typename Ctr::Tds Tds;
  typedef typename Ctr::Itag Itag;

  typedef Constrained_Delaunay_triangulation_2<Gt,Tds_,Itag_>    CDt;
  typedef typename Ctr::Geom_traits      Geom_traits;
  typedef typename Ctr::Intersection_tag Intersection_tag;

  typedef typename Ctr::Constraint    Constraint;
  typedef typename Ctr::Vertex_handle Vertex_handle;
  typedef typename Ctr::Face_handle   Face_handle;
  typedef typename Ctr::Edge          Edge;
  typedef typename Ctr::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Ctr::Constrained_edges_iterator Constrained_edges_iterator;
  typedef typename Ctr::Face_circulator       Face_circulator;
  typedef typename Ctr::size_type             size_type;
  typedef typename Ctr::Locate_type           Locate_type;
 
  typedef typename Ctr::List_edges List_edges;  
  typedef typename Ctr::List_faces List_faces;
  typedef typename Ctr::List_vertices  List_vertices;
  typedef typename Ctr::List_constraints List_constraints;
  typedef typename Ctr::Less_edge less_edge;
  typedef typename Ctr::Edge_set Edge_set;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_false Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Ctr::geom_traits;
  using Ctr::number_of_vertices;
  using Ctr::finite_faces_begin;
  using Ctr::finite_faces_end;
  using Ctr::dimension;
  using Ctr::cw;
  using Ctr::ccw;
  using Ctr::infinite_vertex;
  using Ctr::side_of_oriented_circle;
  using Ctr::is_infinite;
  using Ctr::collinear_between;
  using Ctr::are_there_incident_constraints;
  using Ctr::make_hole;
  // The next using statement makes trouble for VC++ with the version
  // of the method that is templated with PointIterator
  // VC++ cannot disambiguate between this method and the method of the base class
  // using Ctr::insert_constraint;
  using Ctr::locate;
  using Ctr::test_dim_down;
  using Ctr::fill_hole_delaunay;
  using Ctr::update_constraints;
  using Ctr::delete_vertex;
  using Ctr::push_back;
#endif

  typedef typename Geom_traits::Point_2  Point;

  Constrained_Delaunay_triangulation_2(const Geom_traits& gt=Geom_traits()) 
    : Ctr(gt) { }

  Constrained_Delaunay_triangulation_2(const CDt& cdt)
    : Ctr(cdt) {}

  Constrained_Delaunay_triangulation_2(const List_constraints& lc,
				       const Geom_traits& gt=Geom_traits())
    : Ctr(gt) 
    {
      insert_constraints(lc.begin(), lc.end());
      CGAL_triangulation_postcondition( is_valid() );
    }

  template<class InputIterator>
  Constrained_Delaunay_triangulation_2(InputIterator it,
				       InputIterator last,
				       const Geom_traits& gt=Geom_traits() )
    : Ctr(gt) 
    {
      insert_constraints(it, last);
      CGAL_triangulation_postcondition( is_valid() );
    }

  virtual ~Constrained_Delaunay_triangulation_2() {}
  

  // FLIPS
  bool is_flipable(Face_handle f, int i, bool perturb = true) const;
  void flip(Face_handle& f, int i);
  void flip_around(Vertex_handle va);
  void flip_around(List_vertices & new_vertices);
#ifndef CGAL_CDT2_USE_RECURSIVE_PROPAGATING_FLIP
  void non_recursive_propagating_flip(Face_handle f,int i);
  void propagating_flip(Face_handle f,int i, int depth=0);
#else
  void propagating_flip(Face_handle f,int i);
#endif
  void propagating_flip(List_edges & edges);

  // CONFLICTS
  bool test_conflict(Face_handle fh, const Point& p) const; //deprecated
  bool test_conflict(const Point& p, Face_handle fh) const;
  void find_conflicts(const Point& p, std::list<Edge>& le,  //deprecated
		      Face_handle hint= Face_handle()) const;
  //  //template member functions, declared and defined at the end 
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

  // INSERTION-REMOVAL
  Vertex_handle insert(const Point & a, Face_handle start = Face_handle());
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  Vertex_handle push_back(const Point& a);
//   template < class InputIterator >
//   std::ptrdiff_t insert(InputIterator first, InputIterator last);

  void remove(Vertex_handle v);
  void remove_incident_constraints(Vertex_handle v);
  void remove_constrained_edge(Face_handle f, int i);
//  template <class OutputItFaces>
//  OutputItFaces
//  remove_constrained_edge(Face_handle f, int i, OutputItFaces out)
 
  //for backward compatibility
  void insert(Point a, Point b) { insert_constraint(a, b);}

  void insert(Vertex_handle va, Vertex_handle  vb) {insert_constraint(va,vb);}

  void insert_constraint(Vertex_handle va, Vertex_handle  vb)
  {
    ((Ctr*)this)->insert_constraint(va,vb);
  }

  void
  insert_constraint(const Point& a, const Point& b)
  {
    ((Ctr*)this)->insert_constraint(a,b);
  }

  template <class PointIterator>
  void insert_constraint(PointIterator first, PointIterator last, bool close=false)
  {
    if(first == last){
      return;
    }
    const Point& p0 = *first;
    Point p = p0;
    Vertex_handle v0 = insert(p0), v(v0), w(v0);
    ++first;
    for(; first!=last; ++first){
      const Point& q = *first;
      if(p != q){
        w = insert(q);
        insert_constraint(v,w);
        v = w;
        p = q;
      }
    }
    if(close && (p != p0)){
      insert(w,v0);
    }
  }


  void remove_constraint(Face_handle f, int i){remove_constrained_edge(f,i);}

  // CHECK
  bool is_valid(bool verbose = false, int level = 0) const;
 
protected:
  virtual Vertex_handle virtual_insert(const Point & a, 
				       Face_handle start = Face_handle());
  virtual Vertex_handle virtual_insert(const Point& a,
				       Locate_type lt,
				       Face_handle loc, 
				       int li );
//Vertex_handle special_insert_in_edge(const Point & a, Face_handle f, int i);
  void remove_2D(Vertex_handle v );
  virtual void triangulate_hole(List_faces& intersected_faces,
				List_edges& conflict_boundary_ab,
				List_edges& conflict_boundary_ba);

public:
  // MESHING 
  // suppressed meshing functions from here

  //template member functions
public:
  template < class InputIterator >
#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  std::ptrdiff_t
  insert( InputIterator first, InputIterator last,
          typename boost::enable_if<
            boost::is_convertible<
                typename internal::Get_iterator_value_type< InputIterator >::type,
                Point
            >
          >::type* = NULL
  )
#else
#if defined(_MSC_VER)
  std::ptrdiff_t insert(InputIterator first, InputIterator last, int i = 0)
#else
    std::ptrdiff_t insert(InputIterator first, InputIterator last) 
#endif
#endif
    {
      size_type n = number_of_vertices();

      std::vector<Point> points (first, last);
      spatial_sort (points.begin(), points.end(), geom_traits());
      Face_handle f;
      for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
              p != end; ++p)
          f = insert (*p, f)->face();

      return number_of_vertices() - n;
    }

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
private:
  //top stands for tuple-or-pair
  template <class Info>
  const Point& top_get_first(const std::pair<Point,Info>& pair) const { return pair.first; }
  template <class Info>
  const Info& top_get_second(const std::pair<Point,Info>& pair) const { return pair.second; }
  template <class Info>
  const Point& top_get_first(const boost::tuple<Point,Info>& tuple) const { return boost::get<0>(tuple); }
  template <class Info>
  const Info& top_get_second(const boost::tuple<Point,Info>& tuple) const { return boost::get<1>(tuple); }

  template <class Tuple_or_pair,class InputIterator>
  std::ptrdiff_t insert_with_info(InputIterator first,InputIterator last)
  {
    size_type n = this->number_of_vertices();
    std::vector<std::size_t> indices;
    std::vector<Point> points;
    std::vector<typename Tds::Vertex::Info> infos;
    std::size_t index=0;
    for (InputIterator it=first;it!=last;++it){
      Tuple_or_pair value=*it;
      points.push_back( top_get_first(value)  );
      infos.push_back ( top_get_second(value) );
      indices.push_back(index++);
    }

    typedef typename Pointer_property_map<Point>::type Pmap;
    typedef Spatial_sort_traits_adapter_2<Geom_traits,Pmap> Search_traits;

    spatial_sort(indices.begin(),
                 indices.end(),
                 Search_traits(make_property_map(points),geom_traits()));

    Vertex_handle v_hint;
    Face_handle hint;
    for (typename std::vector<std::size_t>::const_iterator
      it = indices.begin(), end = indices.end();
      it != end; ++it){
      v_hint = insert(points[*it], hint);
      if (v_hint!=Vertex_handle()){
        v_hint->info()=infos[*it];
        hint=v_hint->face();
      }
    }

    return this->number_of_vertices() - n;
  }

public:

  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first,
          InputIterator last,
          typename boost::enable_if<
            boost::is_convertible<
              typename internal::Get_iterator_value_type< InputIterator >::type,
              std::pair<Point,typename internal::Info_check<typename Tds::Vertex>::type>
            >
          >::type* =NULL
  )
  {
    return insert_with_info< std::pair<Point,typename internal::Info_check<typename Tds::Vertex>::type> >(first,last);
  }

  template <class  InputIterator_1,class InputIterator_2>
  std::ptrdiff_t
  insert( boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > first,
          boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > last,
          typename boost::enable_if<
            boost::mpl::and_<
              boost::is_convertible< typename std::iterator_traits<InputIterator_1>::value_type, Point >,
              boost::is_convertible< typename std::iterator_traits<InputIterator_2>::value_type, typename internal::Info_check<typename Tds::Vertex>::type >
            >
          >::type* =NULL
  )
  {
    return insert_with_info< boost::tuple<Point,typename internal::Info_check<typename Tds::Vertex>::type> >(first,last);
  }
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO


  template <class PointIterator, class IndicesIterator>
  std::size_t insert_constraints(PointIterator points_first,
                                 PointIterator points_beyond,
                                 IndicesIterator indices_first,
                                 IndicesIterator indices_beyond)
  {
    std::vector<Point> points(points_first, points_beyond);
    return internal::insert_constraints(*this,points, indices_first, indices_beyond);
  }


 template <class ConstraintIterator>
  std::size_t insert_constraints(ConstraintIterator first,
                                 ConstraintIterator beyond)
  {
    return internal::insert_constraints(*this,first,beyond);
  }


  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point  &p, 
			     OutputItFaces fit, 
			     OutputItBoundaryEdges eit,
			     Face_handle start = Face_handle()) const {
    CGAL_triangulation_precondition( dimension() == 2);
    int li;
    Locate_type lt;
    Face_handle fh = locate(p,lt,li, start);
    switch(lt) {
    case Ctr::OUTSIDE_AFFINE_HULL:
    case Ctr::VERTEX:
      return std::make_pair(fit,eit);
    case Ctr::FACE:
    case Ctr::EDGE:
    case Ctr::OUTSIDE_CONVEX_HULL:
      *fit++ = fh; //put fh in OutputItFaces
      std::pair<OutputItFaces,OutputItBoundaryEdges>
	pit = std::make_pair(fit,eit);
      pit = propagate_conflicts(p,fh,0,pit);
      pit = propagate_conflicts(p,fh,1,pit);
      pit = propagate_conflicts(p,fh,2,pit);
      return pit;
    }
    CGAL_triangulation_assertion(false);
    return std::make_pair(fit,eit);
  } 

  template <class OutputItFaces>
  OutputItFaces
  get_conflicts (const Point  &p, 
		 OutputItFaces fit, 
		 Face_handle start= Face_handle()) const {
    std::pair<OutputItFaces,Emptyset_iterator> pp = 
      get_conflicts_and_boundary(p,fit,Emptyset_iterator(),start);
    return pp.first;
  }

  template <class OutputItBoundaryEdges> 
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Point  &p, 
			    OutputItBoundaryEdges eit, 
			    Face_handle start= Face_handle()) const {
    std::pair<Emptyset_iterator, OutputItBoundaryEdges> pp = 
      get_conflicts_and_boundary(p,Emptyset_iterator(),eit,start);
    return pp.second;
  }


public:
// made  public for the need of Mesh_2
// but not documented
#ifdef CGAL_CDT2_USE_RECURSIVE_PROPAGATE_CONFLICTS
  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  propagate_conflicts (const Point  &p,
                      Face_handle fh, 
                      int i,
                      std::pair<OutputItFaces,OutputItBoundaryEdges>  pit)  const 
  {
     Face_handle fn = fh->neighbor(i);
     
     if ( fh->is_constrained(i) || ! test_conflict(p,fn)) {
       *(pit.second)++ = Edge(fn, fn->index(fh));
     } else {
       *(pit.first)++ = fn;
       int j = fn->index(fh);
       pit = propagate_conflicts(p,fn,ccw(j),pit);
       pit = propagate_conflicts(p,fn,cw(j), pit);
     }
     return pit;
  }
#else // NO CGAL_CDT2_USE_RECURSIVE_PROPAGATE_CONFLICTS
  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  non_recursive_propagate_conflicts (const Point  &p,
                                     Face_handle fh, 
                                     int i,
                                     std::pair<OutputItFaces,OutputItBoundaryEdges>  pit)  const 
  {
    std::stack<std::pair<Face_handle, int> > stack;
    stack.push(std::make_pair(fh, i));
    while(!stack.empty()) {
      const Face_handle fh = stack.top().first;
      const int i = stack.top().second;
      stack.pop();
      const Face_handle fn = fh->neighbor(i);
      if ( fh->is_constrained(i) || ! test_conflict(p,fn)) {
        *(pit.second)++ = Edge(fn, fn->index(fh));
      } else {
        *(pit.first)++ = fn;
        int j = fn->index(fh);
        stack.push(std::make_pair(fn, cw(j)));
        stack.push(std::make_pair(fn,ccw(j)));
      }
    }
    return pit;
  }

  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  propagate_conflicts (const Point  &p,
                      Face_handle fh, 
                      int i,
                      std::pair<OutputItFaces,OutputItBoundaryEdges>  pit,
                      int depth=0)  const 
  {
    if ( depth==100) return non_recursive_propagate_conflicts(p,fh,i,pit);
    Face_handle fn = fh->neighbor(i);
 
    if ( fh->is_constrained(i) || ! test_conflict(p,fn)) {
      *(pit.second)++ = Edge(fn, fn->index(fh));
    } else {
      *(pit.first)++ = fn;
      int j = fn->index(fh);
      pit = propagate_conflicts(p,fn,ccw(j),pit,depth+1);
      pit = propagate_conflicts(p,fn,cw(j), pit,depth+1);
    }
     return pit;
  }
#endif // NO CGAL_CDT2_USE_RECURSIVE_PROPAGATE_CONFLICTS




public:
 template <class OutputItFaces>
 OutputItFaces
 propagating_flip(List_edges & edges, 
		  OutputItFaces out = Emptyset_iterator()) {
  // makes the triangulation Delaunay by flipping 
  // List edges contains an initial list of edges to be flipped
  // Precondition : the output triangulation is Delaunay if the list 
  // edges contains all edges of the input triangulation that need to be
  // flipped (plus possibly others)
  int i, ii, indf, indn;
  Face_handle ni, f,ff;
  Edge ei,eni; 
  typename Ctr::Edge_set edge_set;
  typename Ctr::Less_edge less_edge;
  Edge e[4];
  typename List_edges::iterator itedge=edges.begin();

  // initialization of the set of edges to be flip
  while (itedge != edges.end()) {
    f=(*itedge).first;
    i=(*itedge).second;
    if (is_flipable(f,i)) {
      eni=Edge(f->neighbor(i),this->mirror_index(f,i));
      if (less_edge(*itedge,eni)) edge_set.insert(*itedge);
      else edge_set.insert(eni);
    }
    ++itedge;
  }

  // flip edges and updates the set of edges to be flipped
  while (!(edge_set.empty())) {
    f=(*(edge_set.begin())).first;
    indf=(*(edge_set.begin())).second;
 
    // erase from edge_set the 4 edges of the wing of the edge to be
    // flipped (edge_set.begin) , i.e. the edges of the faces f and
    // f->neighbor(indf) that are distinct from the edge to be flipped

    ni = f->neighbor(indf); 
    indn=this->mirror_index(f,indf);
    ei= Edge(f,indf);
    edge_set.erase(ei);
    e[0]= Edge(f,cw(indf));
    e[1]= Edge(f,ccw(indf));
    e[2]= Edge(ni,cw(indn));
    e[3]= Edge(ni,ccw(indn));

    for(i=0;i<4;i++) { 
      ff=e[i].first;
      ii=e[i].second;
      eni=Edge(ff->neighbor(ii),this->mirror_index(ff,ii));
      if (less_edge(e[i],eni)) {edge_set.erase(e[i]);}
      else { edge_set.erase(eni);} 
    } 

    // here is the flip 
    *out++ = f;
    *out++ = f->neighbor(indf);
    flip(f, indf); 
    

    //insert in edge_set the 4 edges of the wing of the edge that
    //have been flipped 
    e[0]= Edge(f,indf);
    e[1]= Edge(f,cw(indf));
    e[2]= Edge(ni,indn);
    e[3]= Edge(ni,cw(indn));

    for(i=0;i<4;i++) { 
      ff=e[i].first;
      ii=e[i].second;
      if (is_flipable(ff,ii)) {
	eni=Edge(ff->neighbor(ii),this->mirror_index(ff,ii));
	if (less_edge(e[i],eni)) { 
	  edge_set.insert(e[i]);}
	else {
	  edge_set.insert(eni);} 
      }
    } 
  }
  return out;
 }

 template <class OutputItFaces>
 OutputItFaces
 remove_constrained_edge(Face_handle f, int i, OutputItFaces out) {
  Ctr::remove_constrained_edge(f,i);
  if(dimension() == 2) {
    List_edges le;
    le.push_back(Edge(f,i));
    propagating_flip(le,out);
  }
  return out;  
 }
 
};


template < class Gt, class Tds, class Itag >
bool 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
is_flipable(Face_handle f, int i, bool perturb /* = true */) const
  // determines if edge (f,i) can be flipped 
{
  Face_handle ni = f->neighbor(i); 
  if (is_infinite(f) || is_infinite(ni)) return false; 
  if (f->is_constrained(i)) return false;
  return (side_of_oriented_circle(ni, f->vertex(i)->point(), perturb) 
                                        == ON_POSITIVE_SIDE);
}

template < class Gt, class Tds, class Itag >
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
flip (Face_handle& f, int i)
{
  Face_handle g = f->neighbor(i);
  int j = this->mirror_index(f,i);

  // save wings neighbors to be able to restore contraint status
  Face_handle f1 = f->neighbor(cw(i));
  int i1 = this->mirror_index(f,cw(i));
  Face_handle f2 = f->neighbor(ccw(i));
  int i2 = this->mirror_index(f,ccw(i));
  Face_handle f3 = g->neighbor(cw(j));
  int i3 = this->mirror_index(g,cw(j));
  Face_handle f4 = g->neighbor(ccw(j));
  int i4 = this->mirror_index(g,ccw(j));

  // The following precondition prevents the test suit 
  // of triangulation to work on constrained Delaunay triangulation
  //CGAL_triangulation_precondition(is_flipable(f,i));
  this->_tds.flip(f, i);
   
  // restore constraint status
  f->set_constraint(f->index(g), false);
  g->set_constraint(g->index(f), false);
  f1->neighbor(i1)->set_constraint(this->mirror_index(f1,i1),
				   f1->is_constrained(i1));
  f2->neighbor(i2)->set_constraint(this->mirror_index(f2,i2),
				   f2->is_constrained(i2));
  f3->neighbor(i3)->set_constraint(this->mirror_index(f3,i3),
				   f3->is_constrained(i3));
  f4->neighbor(i4)->set_constraint(this->mirror_index(f4,i4),
				   f4->is_constrained(i4));
  return;
}

template < class Gt, class Tds, class Itag >
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
flip_around(Vertex_handle va)
  // makes the triangles incident to vertex va Delaunay using flips
{
  if (dimension() <= 1) return;
  Face_handle f=va->face();
  Face_handle next;    
  Face_handle start(f);
  int i;
  do {
    i = f->index(va); // FRAGILE : DIM 1
    next = f->neighbor(ccw(i));  // turns ccw around a
    propagating_flip(f,i);
    f=next;
  } while(next != start);
}

template < class Gt, class Tds, class Itag >
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
flip_around(List_vertices& new_vertices)
{
  typename List_vertices::iterator itv=new_vertices.begin();
  for( ; itv != new_vertices.end(); itv++) {
    flip_around(*itv);
  }
  return;
}


#ifndef CGAL_CDT2_USE_RECURSIVE_PROPAGATING_FLIP
template <class Gt, class Tds, class Itag >
void
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
non_recursive_propagating_flip(Face_handle f , int i)
{
  std::stack<Edge> edges;
  const Vertex_handle& vp = f->vertex(i);
  edges.push(Edge(f,i));

  while(! edges.empty()){
    const Edge& e = edges.top();
    f = e.first;
    i = e.second;

    Face_handle ni = f->neighbor(i);
    flip(f,i);
    if ( !is_flipable(f,i) ) edges.pop();

    i = ni->index(vp);
    if ( is_flipable(ni,i) ) edges.push( Edge(ni,i) );
  }
}

template <class Gt, class Tds, class Itag >
void
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
propagating_flip(Face_handle f,int i, int depth)
{
  if (!is_flipable(f,i)) return;
#ifdef CGAL_CDT2_IMMEDIATELY_NON_RECURSIVE_PROPAGATING_FLIP
  non_recursive_propagating_flip(f,i);
#else
  int max_depth = 100;
  if(depth==max_depth){
    non_recursive_propagating_flip(f,i);
    return;
  }

  Face_handle ni = f->neighbor(i);
  flip(f, i); // flip for constrained triangulations
  propagating_flip(f,i, depth+1);
  i = ni->index(f->vertex(i));
  propagating_flip(ni,i, depth+1);
#endif
}
#else 
template < class Gt, class Tds, class Itag >
void
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
propagating_flip(Face_handle f,int i)
// similar to the corresponding function in Delaunay_triangulation_2.h
{
  if (!is_flipable(f,i)) return;
  Face_handle ni = f->neighbor(i);
  flip(f, i); // flip for constrained triangulations
  propagating_flip(f,i);
  i = ni->index(f->vertex(i));
  propagating_flip(ni,i);
}
#endif

 template < class Gt, class Tds, class Itag > 
 void
 Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>:: 
 propagating_flip(List_edges & edges) {
    propagating_flip(edges,Emptyset_iterator());
 }


template < class Gt, class Tds, class Itag >
inline bool
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
test_conflict(const Point& p, Face_handle fh) const
  // true if point P lies inside the circle circumscribing face fh
{
  // return true  if P is inside the circumcircle of fh
  // if fh is infinite, return true when p is in the positive
  // halfspace or on the boundary and in the  finite edge of fh
  Oriented_side os = side_of_oriented_circle(fh,p,true);
  if (os == ON_POSITIVE_SIDE) return true;
 
  if (os == ON_ORIENTED_BOUNDARY && is_infinite(fh)) {
    int i = fh->index(infinite_vertex());
    return collinear_between(fh->vertex(cw(i))->point(), p,
			     fh->vertex(ccw(i))->point() );
  }

  return false;
}

template < class Gt, class Tds, class Itag >
inline bool
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
test_conflict(Face_handle fh, const Point& p) const
  // true if point P lies inside the circle circumscribing face fh
{
  return test_conflict(p,fh);
}

template < class Gt, class Tds, class Itag >    
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
find_conflicts(const Point& p, std::list<Edge>& le, Face_handle hint) const
{
  // sets in le the counterclocwise list of the edges of the boundary of the 
  // union of the faces in conflict with p
  // an edge is represented by the incident face that is not in conflict with p
  get_boundary_of_conflicts(p, std::back_inserter(le), hint);
}
  
template < class Gt, class Tds, class Itag >  
inline 
typename Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
virtual_insert(const Point & a, Face_handle start)
  // virtual version of the insertion
{
  return insert(a,start);
}

template < class Gt, class Tds, class Itag >  
inline 
typename Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
virtual_insert(const Point& a,
	       Locate_type lt,
	       Face_handle loc, 
	       int li )
// virtual version of insert
{
  return insert(a,lt,loc,li);
}

template < class Gt, class Tds, class Itag >  
inline 
typename Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
insert(const Point & a, Face_handle start)
 // inserts a in the triangulation
// constrained edges are updated
// Delaunay property is restored
{
  Vertex_handle va= Ctr::insert(a, start);
  flip_around(va); 
  return va;
}

template < class Gt, class Tds, class Itag >
typename Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
insert(const Point& a, Locate_type lt, Face_handle loc, int li)
// insert a point p, whose localisation is known (lt, f, i)
// constrained edges are updated
// Delaunay property is restored
{
  Vertex_handle va= Ctr::insert(a,lt,loc,li);
  flip_around(va); 
  return va;
}

template < class Gt, class Tds, class Itag >
inline
typename Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
push_back(const Point &p)
{
  return insert(p);
}

template < class Gt, class Tds, class Itag >
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
triangulate_hole(List_faces& intersected_faces,
		 List_edges& conflict_boundary_ab,
		 List_edges& conflict_boundary_ba)
{
  List_edges new_edges;
  Ctr::triangulate_hole(intersected_faces,
		       conflict_boundary_ab,
		       conflict_boundary_ba,
		       new_edges);
  propagating_flip(new_edges);
  return;
}


template < class Gt, class Tds, class Itag >  
inline void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
remove(Vertex_handle v)
  // remove a vertex and updates the constrained edges of the new faces
  // precondition : there is no incident constraints
{
  CGAL_triangulation_precondition( v != Vertex_handle() );
  CGAL_triangulation_precondition( ! is_infinite(v));
  CGAL_triangulation_precondition( ! are_there_incident_constraints(v));
  if  (dimension() <= 1)    Ctr::remove(v);
  else  remove_2D(v);
  return;
}

// template < class Gt, class Tds, class Itag >  
// typename
// Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::Vertex_handle
// Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
// special_insert_in_edge(const Point & a, Face_handle f, int i)
//   // insert  point p in edge(f,i)
//   // bypass the precondition for point a to be in edge(f,i)
//   // update constrained status
//   // this member fonction is not robust with exact predicates 
//   // and approximate construction. Should be removed 
// {
//   Vertex_handle vh=Ctr::special_insert_in_edge(a,f,i);
//   flip_around(vh);
//   return vh;
// }

template < class Gt, class Tds, class Itag >  
void 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
remove_2D(Vertex_handle v)
{
 if (test_dim_down(v)) {  this->_tds.remove_dim_down(v);  }
  else {
     std::list<Edge> hole;
    make_hole(v, hole);
    std::list<Edge> shell=hole; //because hole will be emptied by fill_hole
    fill_hole_delaunay(hole);
    update_constraints(shell);
    delete_vertex(v);
  }
  return;
}

template < class Gt, class Tds, class Itag >
void
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
remove_incident_constraints(Vertex_handle v)
{
   List_edges iconstraints;
   if (are_there_incident_constraints(v,
				      std::back_inserter(iconstraints))) {
     Ctr::remove_incident_constraints(v);
     if (dimension()==2) propagating_flip(iconstraints);
   }
   return;
}

template < class Gt, class Tds, class Itag >
void
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
remove_constrained_edge(Face_handle f, int i)
{
  remove_constrained_edge(f,i,Emptyset_iterator());
  return;  
}


template < class Gt, class Tds, class Itag >  
bool 
Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>::
is_valid(bool verbose, int level) const
{
  bool result = Ctr::is_valid(verbose, level);
  CGAL_triangulation_assertion( result );

    Finite_faces_iterator fit= finite_faces_begin();
    for (; fit != finite_faces_end(); fit++) {
      for(int i=0;i<3;i++) {
	result = result && !is_flipable(fit,i, false);
	CGAL_triangulation_assertion( result );
      }
    }
    return result;
}



} //namespace CGAL
#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H
