// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec, Sebastien Loriot

#ifndef CGAL_INTERNAL_TRIANGULATION_EUCLIDEAN_TRAITS_PROJECTED_3_H
#define CGAL_INTERNAL_TRIANGULATION_EUCLIDEAN_TRAITS_PROJECTED_3_H

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/predicates/kernel_ftC2.h>

namespace CGAL { 

namespace internal {

//project Point_3 along coordinate dim
template <class R,int dim>
struct Projector;

//project onto yz
template <class R>
struct Projector<R,0>
{
  typedef typename R::Less_y_3                Less_x_2;
  typedef typename R::Less_z_3                Less_y_2;
  typedef typename R::Compare_y_3             Compare_x_2;
  typedef typename R::Compare_z_3             Compare_y_2;
  typedef typename R::Equal_y_3               Equal_x_2;
  typedef typename R::Equal_z_3               Equal_y_2;  
  
  static typename R::FT x(const typename R::Point_3& p) {return p.y();}
  static typename R::FT y(const typename R::Point_3& p) {return p.z();}
  static const int x_index=1;
  static const int y_index=2;
};
//project onto xz
template <class R>
struct Projector<R,1>
{
  typedef typename R::Less_x_3                Less_x_2;
  typedef typename R::Less_z_3                Less_y_2;
  typedef typename R::Compare_x_3             Compare_x_2;
  typedef typename R::Compare_z_3             Compare_y_2;  
  typedef typename R::Equal_x_3               Equal_x_2;
  typedef typename R::Equal_z_3               Equal_y_2;    
  static typename R::FT x(const typename R::Point_3& p) {return p.x();}
  static typename R::FT y(const typename R::Point_3& p) {return p.z();}
  static const int x_index=0;
  static const int y_index=2;  
};
//project onto xy
template <class R>
struct Projector<R,2>
{
  typedef typename R::Less_x_3                Less_x_2;
  typedef typename R::Less_y_3                Less_y_2;
  typedef typename R::Compare_x_3             Compare_x_2;
  typedef typename R::Compare_y_3             Compare_y_2;  
  typedef typename R::Equal_x_3               Equal_x_2;
  typedef typename R::Equal_y_3               Equal_y_2;    
  static typename R::FT x(const typename R::Point_3& p) {return p.x();}
  static typename R::FT y(const typename R::Point_3& p) {return p.y();}
  static const int x_index=0;
  static const int y_index=1;  
};
  
template <class R,int dim>
class Orientation_projected_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point &p) const { return Projector<R,dim>::y(p); }

  CGAL::Orientation operator()(const Point& p,
			       const Point& q,
			       const Point& r) const
    {
      return orientationC2(x(p), y(p), x(q), y(q), x(r), y(r));
    }
};

template <class R,int dim>
class Side_of_oriented_circle_projected_3 
{
public:
  typedef typename R::Point_3     Point; 
  typename R::FT x(const Point &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point &p) const { return Projector<R,dim>::y(p); }

  CGAL::Oriented_side operator() (const Point &p, 
				  const Point &q,
				  const Point &r, 
				  const Point &s) const
    {
      return side_of_oriented_circleC2(x(p), y(p),
				       x(q), y(q),
				       x(r), y(r),
				       x(s), y(s));
    }
};

template <class R,int dim>
class Compare_distance_projected_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;   
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  Comparison_result operator()(const Point_3& p,const Point_3& q,const Point_3& r) const
  {
    Point_2 p2 = project(p);
    Point_2 q2 = project(q);
    Point_2 r2 = project(r);
    return compare_distance_to_point(p2,q2,r2);
  }
};

template <class R,int dim>
class Squared_distance_projected_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Line_3    Line_3; 
  typedef typename R::Line_2    Line_2;
  typedef typename R::FT        RT;
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  RT operator()(const Line_3& l, const Point_3& p) const
  {
    Point_2 p2(project(p));
    Line_2 l2(project(l.point(0)), project(l.point(1)));
    return squared_distance(p2, l2);
  }
};

template <class R,int dim>
class  Intersect_projected_3
{
public:
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Segment_3 Segment_3;
  typedef typename R::Point_2   Point_2; 
  typedef typename R::Vector_2  Vector_2; 
  typedef typename R::Segment_2 Segment_2;
  typedef typename R::FT        FT;
  
  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }

  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  Object operator()(const Segment_3& s1, const Segment_3& s2) const
  {
    Point_2 s1_source = project(s1.source());
    Point_2 s1_target = project(s1.target());
    Point_2 s2_source = project(s2.source());
    Point_2 s2_target = project(s2.target());    
    Segment_2 s1_2(s1_source, s1_target);
    Segment_2 s2_2(s2_source, s2_target);
    CGAL_precondition(!s1_2.is_degenerate());
    CGAL_precondition(!s2_2.is_degenerate());
    
    //compute intersection points in projected plane
    //We know that none of the segment is degenerate
    Object o = intersection(s1_2,s2_2);
    const Point_2* pi=CGAL::object_cast<Point_2>(&o);
    if (pi==NULL) { //case of segment or empty
      const Segment_2* si=CGAL::object_cast<Segment_2>(&o);
      if (si==NULL) return Object();
      FT src[3],tgt[3];
      tgt[dim] = src[dim] = FT(0); //the third coordinate is arbitrarily set to 0 (not used in Constrained DT)
      src[Projector<R,dim>::x_index] = si->source().x();
      src[Projector<R,dim>::y_index] = si->source().y();
      tgt[Projector<R,dim>::x_index] = si->target().x();
      tgt[Projector<R,dim>::y_index] = si->target().y();      
      return make_object( Segment_3( Point_3(src[0],src[1],src[2]),Point_3(tgt[0],tgt[1],tgt[2]) ) );
    }
    FT coords[3];
    //compute the third coordinate of the projected intersection point onto 3D segments
    FT z1 = s1.source()[dim] + ( pi->x()-s1_source.x() ) / (s1_target.x() - s1_source.x()) * ( s1.target()[dim] - s1.source()[dim] );
    FT z2 = s2.source()[dim] + ( pi->x()-s2_source.x() ) / (s2_target.x() - s2_source.x()) * ( s2.target()[dim] - s2.source()[dim] );
    coords[dim] = (z1+z2) / FT(2);
    coords[Projector<R,dim>::x_index] = pi->x();
    coords[Projector<R,dim>::y_index] = pi->y();
    
    Point_3 res(coords[0],coords[1],coords[2]);
    CGAL_assertion(x(res)==pi->x() && y(res)==pi->y());
    return make_object(res);
  }
};

template <class R, int dim>
class Circumcenter_center_projected
{
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }
  
  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }
  
  Point_3 embed (const Point_2& p) const 
  {
    typename R::FT coords[3];
    coords[Projector<R,dim>::x_index]=p.x();
    coords[Projector<R,dim>::y_index]=p.y();
    coords[dim]=typename R::FT(0);
    return Point_3(coords[0],coords[1],coords[2]);
  }
  
public:
  Point_3 operator() (const Point_3& p1,const Point_3& p2) const
  {
    return embed( circumcenter(project(p1),project(p2)) );
  }

  Point_3 operator() (const Point_3& p1,const Point_3& p2,const Point_3& p3) const
  {
    return embed( circumcenter(project(p1),project(p2),project(p3)) );
  }
};

template <class R, int dim>
class Compute_area_projected
{
  typedef typename R::Point_3   Point_3; 
  typedef typename R::Point_2   Point_2;

  typename R::FT x(const Point_3 &p) const { return Projector<R,dim>::x(p); }
  typename R::FT y(const Point_3 &p) const { return Projector<R,dim>::y(p); }
  
  Point_2 project(const Point_3& p) const
  {
    return Point_2(x(p),y(p));
  }

  
public:
  typename R::FT operator() (const Point_3& p1,const Point_3& p2,const Point_3& p3) const
  {
    return R().compute_area_2_object() ( project(p1),project(p2),project(p3) );
  }
};

template < class R, int dim >
class Triangulation_euclidean_traits_projected_3 {
public:
  typedef Triangulation_euclidean_traits_projected_3<R,dim>   Traits;
  typedef R                                                   Rp;
  typedef typename R::FT                                      FT;
  typedef typename Rp::Point_3                                Point_2;
  typedef typename Rp::Segment_3                              Segment_2;
  typedef typename Rp::Triangle_3                             Triangle_2;
  typedef typename Rp::Line_3                                 Line_2;

  typedef typename Projector<R,dim>::Less_x_2                 Less_x_2;
  typedef typename Projector<R,dim>::Less_y_2                 Less_y_2;
  typedef typename Projector<R,dim>::Compare_x_2              Compare_x_2;
  typedef typename Projector<R,dim>::Compare_y_2              Compare_y_2;
  typedef Orientation_projected_3<Rp,dim>                     Orientation_2;
  typedef Side_of_oriented_circle_projected_3<Rp,dim>         Side_of_oriented_circle_2;
  typedef Compare_distance_projected_3<Rp,dim>                Compare_distance_2;
  typedef Squared_distance_projected_3<Rp,dim>                Compute_squared_distance_2;
  typedef Intersect_projected_3<Rp,dim>                       Intersect_2;
  typedef typename Rp::Construct_segment_3                    Construct_segment_2;
  typedef typename Rp::Construct_triangle_3                   Construct_triangle_2;
  typedef typename Rp::Construct_line_3                       Construct_line_2;

  //for natural_neighbor_coordinates_2
  typedef typename Projector<R,dim>::Equal_x_2                Equal_x_2;
  typedef Circumcenter_center_projected<Rp,dim>               Construct_circumcenter_2;
  typedef Compute_area_projected<Rp,dim>                      Compute_area_2;
  Construct_circumcenter_2 construct_circumcenter_2_object () const {return Construct_circumcenter_2();}
  Compute_area_2 compute_area_2_object () const {return Compute_area_2();}  
  

  // for compatibility with previous versions
  typedef Point_2      Point;
  typedef Segment_2    Segment;
  typedef Triangle_2   Triangle;

  Triangulation_euclidean_traits_projected_3(){}
  Triangulation_euclidean_traits_projected_3(
		   const Triangulation_euclidean_traits_projected_3&){}
  Triangulation_euclidean_traits_projected_3 &operator=(
	    const Triangulation_euclidean_traits_projected_3&){return *this;}

  typename Rp::FT x(const Point_2 &p) const { return Projector<R,dim>::x(p); }
  typename Rp::FT y(const Point_2 &p) const { return Projector<R,dim>::y(p); }
    
 
  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
    {return Side_of_oriented_circle_2();}

  Compare_distance_2
  compare_distance_2_object() const
  {
    return Compare_distance_2();
  }

  Compute_squared_distance_2
  compute_squared_distance_2_object () const
  {
    return Compute_squared_distance_2();
  }

  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2();
  }

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}
    
  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}
   
};
  

} } //namespace CGAL::internal

#endif // CGAL_INTERNAL_TRIANGULATION_EUCLIDEAN_TRAITS_PROJECTED_3_H
