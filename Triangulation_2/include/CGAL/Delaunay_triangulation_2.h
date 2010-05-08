// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Mariette Yvinec



#ifndef CGAL_DELAUNAY_TRIANGULATION_2_H
#define CGAL_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Triangulation_2.h>
#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                           Triangulation_vertex_base_2<Gt> > >
class Delaunay_triangulation_2 : public Triangulation_2<Gt,Tds>
{
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Point_2       Point;
  typedef typename Geom_traits::Segment_2     Segment;
  typedef typename Geom_traits::Triangle_2    Triangle;
  

  typedef typename Geom_traits::Orientation_2 Orientation_2;
  typedef typename Geom_traits::Compare_x_2   Compare_x;
  typedef typename Geom_traits::Compare_y_2   Compare_y;
  typedef typename Geom_traits::Side_of_oriented_circle_2 
                                              Side_of_oriented_circle;

  
  typedef Triangulation_2<Gt,Tds>                       Triangulation;
  typedef typename Triangulation::Locate_type           Locate_type;
  typedef typename Triangulation::Face_handle           Face_handle;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;
  typedef typename Triangulation::Edge                  Edge;
  typedef typename Triangulation::Edge_circulator       Edge_circulator;
  typedef typename Triangulation::Face_circulator       Face_circulator;
  typedef typename Triangulation::Vertex_circulator     Vertex_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Triangulation::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Triangulation::Finite_vertices_iterator 
                                                     Finite_vertices_iterator;

 
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::side_of_oriented_circle;
  using Triangulation::circumcenter;
  using Triangulation::collinear_between;
  using Triangulation::test_dim_down;
  using Triangulation::make_hole;
  using Triangulation::fill_hole_delaunay;
  using Triangulation::delete_vertex;
#endif


 Delaunay_triangulation_2(const Gt& gt = Gt())
  : Triangulation_2<Gt,Tds>(gt) {}
  
  Delaunay_triangulation_2(
	       const Delaunay_triangulation_2<Gt,Tds> &tr)
       : Triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition( is_valid() );  }
  
// CHECK -QUERY
  bool is_valid(bool verbose = false, int level = 0) const;

  Vertex_handle
  nearest_vertex(const Point& p, Face_handle f= Face_handle()) const;
  
  bool does_conflict(const Point  &p, Face_handle fh) const;// deprecated
  bool test_conflict(const Point  &p, Face_handle fh) const;
  bool find_conflicts(const Point  &p,                //deprecated
		      std::list<Face_handle>& conflicts,
		      Face_handle start= Face_handle() ) const;
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
   
 
  // DUAL
  Point dual (Face_handle f) const;
  Object dual(const Edge &e) const ;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;
  
  //INSERTION-REMOVAL
  Vertex_handle insert(const Point  &p, 
		       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  Vertex_handle push_back(const Point &p);

  void  remove(Vertex_handle v );
  void restore_Delaunay(Vertex_handle v);

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::geom_traits;
#endif

private:
  void propagating_flip(Face_handle& f,int i);
  void remove_2D(Vertex_handle v );

  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const;
  Vertex_handle nearest_vertex_1D(const Point& p) const;

  void  look_nearest_neighbor(const Point& p,
			      Face_handle f,
			      int i,
			      Vertex_handle& nn) const;

public:
  template < class Stream>
  Stream& draw_dual(Stream & ps)
    {
      Finite_edges_iterator eit= this->finite_edges_begin();
      for (; eit != this->finite_edges_end(); ++eit) {
	Object o = dual(eit);
	typename Geom_traits::Line_2  l;
	typename Geom_traits::Ray_2   r;
	Segment s;
	if (CGAL::assign(s,o)) ps << s;
	if (CGAL::assign(r,o)) ps << r;
	if (CGAL::assign(l,o)) ps << l;
      }
      return ps;
    }

  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
    {
      int n = this->number_of_vertices();

      std::vector<Point> points (first, last);
      spatial_sort (points.begin(), points.end(), geom_traits());
      Face_handle f;
      for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
              p != end; ++p)
          f = insert (*p, f)->face();

      return this->number_of_vertices() - n;
    }

  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point  &p, 
			     OutputItFaces fit, 
			     OutputItBoundaryEdges eit,
			     Face_handle start = Face_handle()) const {
    CGAL_triangulation_precondition( this->dimension() == 2);
    int li;
    Locate_type lt;
    Face_handle fh = this->locate(p,lt,li, start);
    switch(lt) {
    case Triangulation::OUTSIDE_AFFINE_HULL:
    case Triangulation::VERTEX:
      return std::make_pair(fit,eit);
    case Triangulation::FACE:
    case Triangulation::EDGE:
    case Triangulation::OUTSIDE_CONVEX_HULL:
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

private:
  template <class OutputItFaces, class OutputItBoundaryEdges> 
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  propagate_conflicts (const Point  &p,
		       Face_handle fh, 
		       int i,
		       std::pair<OutputItFaces,OutputItBoundaryEdges>
		       pit)  const {
    Face_handle fn = fh->neighbor(i);
    if (! test_conflict(p,fn)) {
      *(pit.second)++ = Edge(fn, fn->index(fh));
    } else {
      *(pit.first)++ = fn;
      int j = fn->index(fh);
      pit = propagate_conflicts(p,fn,ccw(j),pit);
      pit = propagate_conflicts(p,fn,cw(j), pit);
    }
    return pit;
  }

protected:

  void restore_edges(Vertex_handle v)
  {
    std::list<Edge> edges;
    Face_circulator fc = this->incident_faces(v), done(fc);
    int degree = 0;
    do {
      if((++degree) > 3) break;
    } while(++fc != done);
    fc = this->incident_faces(v);
    done = fc;
    if(degree == 3) {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
      } while(++fc != done);
    } else {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
        edges.push_back(Edge(fc, this->cw(i)));
      } while(++fc != done);
    }
    while(!edges.empty()) {
      const Edge &e = edges.front();
      Face_handle f = e.first;
      int i = e.second;
      edges.pop_front();
      if(this->is_infinite(f->vertex(i))) continue;
      Face_handle fi = f->neighbor(i);
      int mi = this->_tds.mirror_index(f, i);
      Vertex_handle vm = this->_tds.mirror_vertex(f, i);
      if(this->is_infinite(vm)) continue;
      if(this->side_of_oriented_circle(f, vm->point(),true) == ON_POSITIVE_SIDE) {
        this->_tds.flip(f, i);
        edges.push_back(Edge(f, i));
        edges.push_back(Edge(f, this->cw(i)));
        edges.push_back(Edge(fi, this->cw(mi)));
        edges.push_back(Edge(fi, mi));
      }
    }
  }


};

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
test_conflict(const Point  &p, Face_handle fh) const
{
  // return true  if P is inside the circumcircle of fh
  // if fh is infinite, return true when p is in the positive
  // halfspace or on the boundary and in the  finite edge of fh
  Oriented_side os = side_of_oriented_circle(fh,p,true);
  if (os == ON_POSITIVE_SIDE) return true;
 
  if (os == ON_ORIENTED_BOUNDARY && this->is_infinite(fh)) {
    int i = fh->index(this->infinite_vertex());
    return collinear_between(fh->vertex(cw(i))->point(), p,
			     fh->vertex(ccw(i))->point() );
  }

  return false;
}

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
does_conflict(const Point  &p, Face_handle fh) const
{
  return test_conflict(p,fh);
}

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
find_conflicts(const Point  &p, 
	       std::list<Face_handle>& conflicts,
	       Face_handle start ) const
{
  get_conflicts(p, std::back_inserter(conflicts), start);
  return (! conflicts.empty());
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_2<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  bool result = Triangulation_2<Gt,Tds>::is_valid(verbose, level);

  for( Finite_faces_iterator it = this->finite_faces_begin(); 
       it != this->finite_faces_end() ; it++) {
    for(int i=0; i<3; i++) {
      if ( ! this->is_infinite( this->mirror_vertex(it,i))) {
	result = result &&  ON_POSITIVE_SIDE != 
	  side_of_oriented_circle( it, this->mirror_vertex(it,i)->point(), false);
      }
      CGAL_triangulation_assertion( result );
    }
  }
  return result;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex(const Point  &p, Face_handle f) const
{
  switch (this->dimension()) {
  case 0:
    if (this->number_of_vertices() == 0) return Vertex_handle();
    if (this->number_of_vertices() == 1) return this->finite_vertex();
    //break;
  case 1:
    return nearest_vertex_1D(p);
    //break;      
  case 2:
    return nearest_vertex_2D(p,f);
    //break;
  }
  return Vertex_handle();
}
  
template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_2D(const Point& p, Face_handle f) const
{
  CGAL_triangulation_precondition(this->dimension() == 2);
  f = this->locate(p,f);

  typename Geom_traits::Compare_distance_2 
    compare_distance =  this->geom_traits().compare_distance_2_object();
  Vertex_handle nn =  !this->is_infinite(f->vertex(0)) ? f->vertex(0):f->vertex(1);
  if ( !this->is_infinite(f->vertex(1)) && compare_distance(p,
					    f->vertex(1)->point(),
					    nn->point()) == SMALLER) 
    nn=f->vertex(1);
  if ( !this->is_infinite(f->vertex(2)) && compare_distance(p,
					    f->vertex(2)->point(), 
					    nn->point()) == SMALLER) 
    nn=f->vertex(2);
       
  look_nearest_neighbor(p,f,0,nn);
  look_nearest_neighbor(p,f,1,nn);
  look_nearest_neighbor(p,f,2,nn);
  return nn;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>:: 
nearest_vertex_1D(const Point& p) const
{
  typename Geom_traits::Compare_distance_2 
    compare_distance =  this->geom_traits().compare_distance_2_object();
  Vertex_handle nn;
  
  Finite_vertices_iterator vit=this->finite_vertices_begin();
  nn = vit;
  for ( ; vit != this->finite_vertices_end(); ++vit){
    if (compare_distance(p, vit->point(), nn->point()) == SMALLER) 
      nn = vit;
  } 
  return nn;
}
  
template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
look_nearest_neighbor(const Point& p,
                      Face_handle f,
		      int i,
		      Vertex_handle& nn) const
{
  Face_handle  ni=f->neighbor(i);
  if ( ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p,true) ) return;

  typename Geom_traits::Compare_distance_2 
    compare_distance =  this->geom_traits().compare_distance_2_object();
  i = ni->index(f);
  if ( !this->is_infinite(ni->vertex(i)) &&
       compare_distance(p, 
	      ni->vertex(i)->point(),
	      nn->point())  == SMALLER)  nn=ni->vertex(i);
    
  // recursive exploration of triangles whose circumcircle contains p
  look_nearest_neighbor(p, ni, ccw(i), nn);
  look_nearest_neighbor(p, ni, cw(i), nn);
} 

//DUALITY
template<class Gt, class Tds>
inline
typename Delaunay_triangulation_2<Gt,Tds>::Point
Delaunay_triangulation_2<Gt,Tds>::
dual (Face_handle f) const
{
  CGAL_triangulation_precondition (this->dimension()==2);
  return circumcenter(f);
}

  
template < class Gt, class Tds >
Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Edge &e) const
{
  typedef typename Geom_traits::Line_2        Line;
  typedef typename Geom_traits::Ray_2         Ray;

  CGAL_triangulation_precondition (!this->is_infinite(e));
  if( this->dimension()== 1 ){
    const Point& p = (e.first)->vertex(cw(e.second))->point();
    const Point& q = (e.first)->vertex(ccw(e.second))->point();
    Line l  = this->geom_traits().construct_bisector_2_object()(p,q);
    return make_object(l);
  }
		    
  // dimension==2
  if( (!this->is_infinite(e.first)) &&
      (!this->is_infinite(e.first->neighbor(e.second))) ) {
    Segment s = this->geom_traits().construct_segment_2_object()
                          (dual(e.first),dual(e.first->neighbor(e.second)));
    return make_object(s);
  } 
  // one of the adjacent faces is infinite
  Face_handle f; int i;
  if (this->is_infinite(e.first)) {
    f=e.first->neighbor(e.second); i=f->index(e.first);
  }
  else {
    f=e.first; i=e.second;
  }
  const Point& p = f->vertex(cw(i))->point();
  const Point& q = f->vertex(ccw(i))->point();
  Line l = this->geom_traits().construct_bisector_2_object()(p,q);
  Ray r = this->geom_traits().construct_ray_2_object()(dual(f), l);
  return make_object(r);
}
  
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::  
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}
  
template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Finite_edges_iterator& ei) const
{
  return dual(*ei);
}
  
template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p,  Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = this->locate (p, lt, li, start);
  return insert(p, lt, loc, li);
}
  
template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  return insert(p);
}
  
template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v = Triangulation_2<Gt,Tds>::insert(p,lt,loc,li);
  restore_Delaunay(v);
  return(v);
}


template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
restore_Delaunay(Vertex_handle v)
{
  if(this->dimension() <= 1) return;

  Face_handle f=v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do {
    i = f->index(v);
    next = f->neighbor(ccw(i));  // turn ccw around v
    propagating_flip(f,i);
    f=next;
  } while(next != start);
  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !this->is_infinite(v));

  if ( this->dimension() <= 1) Triangulation::remove(v);
  else  remove_2D(v);

  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
propagating_flip(Face_handle& f,int i)
{
  Face_handle n = f->neighbor(i);
      
  if ( ON_POSITIVE_SIDE != 
       side_of_oriented_circle(n,  f->vertex(i)->point(), true) ) {          
    return;
  }
  this->flip(f, i);
  propagating_flip(f,i);
  i = n->index(f->vertex(i));
  propagating_flip(n,i);
}

template <class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) {  this->_tds.remove_dim_down(v);  }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole_delaunay(hole);
    delete_vertex(v);
  }
  return;       
}

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H
