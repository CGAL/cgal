// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Julia Floetotto
#ifndef CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H
#define CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H

#include <CGAL/tags.h>

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/predicates_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/constructions/constructions_for_voronoi_intersection_cartesian_2_3.h>

#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/predicates_for_voronoi_intersection_cartesian_2_3.h>
#include <CGAL/constructions/constructions_for_voronoi_intersection_cartesian_2_3.h>

#endif

CGAL_BEGIN_NAMESPACE 

template <class Traits>
class Orientation_with_normal_plane_2_3
{
public:
  typedef typename Traits::Point_3     Point;
  typedef typename Traits::Vector_3    Vector;


  Orientation_with_normal_plane_2_3(const Vector& _normal,
				    const Traits& _traits)
    :normal(_normal), traits(_traits){};
  
   Orientation operator() 
    (const Point &p,
     const Point  &q,
     const Point  &r) const
    {
      return traits.orientation_3_object()(p,q,q+normal,r);
    }
  
private:
  Vector  normal;
  const Traits&  traits;
};

template <class Point>
class Side_of_plane_centered_sphere_2_3
{
public:
  typedef Oriented_side  result_type;
  typedef typename Point::R    Rep;
  typedef typename Rep::Vector_3     Vector;
 
  typedef typename Rep::Plane_3   Plane;
  typedef typename Rep::Direction_3  Direction;

  Side_of_plane_centered_sphere_2_3(const Point& _p, 
				    const Vector& _normal)
    : a(_p), normal(_normal){};
  
  Oriented_side operator() 
    (const Point &p,
     const Point  &q,
     const Point  &r,
     const Point  &t) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q,r,t);
    }

private:
  Point   a;
  Vector  normal;
};




template <class Point>
class Side_of_plane_centered_sphere_degenerated_2_3
{
public:
  typedef typename Point::R::Vector_3     Vector;
   

  Side_of_plane_centered_sphere_degenerated_2_3(const Point& _p, 
						const Vector& _normal)
    : a(_p), normal(_normal){};
  
  Oriented_side operator() 
    (const Point &p,
     const Point &q,
     const Point &r) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q,r);
    }

  
  Oriented_side operator() 
    (const Point &p,
     const Point &q) const
    {
      return side_of_plane_centered_sphere(a,normal,p,q);
    }

private:
  Point   a;
  Vector  normal;
};

template <class Point>
class Construct_plane_centered_circumcenter_3
{
public:
  typedef typename Point::R Rep;
  typedef typename Rep::Vector_3     Vector;
 
  Construct_plane_centered_circumcenter_3(const Point& _p, 
					  const Vector& _normal)
    : a(_p), normal(_normal){};
  
  Point operator() 
    (const Point &p, const Point &q, const Point &r) const
    {
      return plane_centered_circumcenter_3(a,normal,p,q,r);
    }

private:
  Point   a;
  Vector  normal;
};


template <class Point>
class Construct_plane_intersected_bisector_3
{
public:
  typedef typename Point::R          Rep;
  typedef typename Rep::Vector_3     Vector;
  typedef typename Rep::Line_3       Line;
  
  Construct_plane_intersected_bisector_3(const Point& _a,
					 const Vector& _normal)
    : a(_a), normal(_normal){};
  
  Line operator() (const Point &p, const Point  &q) const
    {
      return plane_intersected_bisector_3(a,normal, p,q);
    }
  
private:
  Point a;
  Vector  normal;
};


template <class Point>
class Compare_first_projection_3
{
  //compares the projection of two points onto a second (non-trivial) 
  // vector in the projection plane:
public:
  typedef typename Point::R          Rep;
  typedef typename Rep::Vector_3     Vector;
  typedef typename Rep::FT           Coord_type;

  Compare_first_projection_3(const Vector& _normal)
    : normal(_normal){};
  
  Comparison_result operator() (const Point &p, const Point  &q) const
    {
      if(normal.x()!=Coord_type(0))
	return (Comparison_result) CGAL_NTS 
	  sign(Vector(normal.y(),-normal.x(),Coord_type(0))*(p-q));
      if(normal.y()!= Coord_type(0))
	return (Comparison_result) CGAL_NTS 
	  sign(Vector(-normal.y(),normal.x(),Coord_type(0))*(p-q));
      
      CGAL_assertion(normal.z()!= Coord_type(0));
      return (Comparison_result) CGAL_NTS 
	sign(Vector(-normal.z(),Coord_type(0),normal.x())*(p-q));
    }
  
private:
  Vector normal;
};
template <class Point>
class Compare_second_projection_3
{
  //compares the projection of two points onto a second (non-trivial) 
  // vector in the projection plane:
public:
  typedef typename Point::R          Rep;
  typedef typename Rep::FT           Coord_type;
  typedef typename Rep::Vector_3     Vector;
  
  Compare_second_projection_3(const Vector& _normal)
    : normal(_normal){};
  
  Comparison_result operator() (const Point &p, const Point  &q) const
    {
      if(normal.x()!=Coord_type(0))
	return (Comparison_result) CGAL_NTS 
	  sign(Vector(normal.z(),Coord_type(0),-normal.x())*(p-q));
      if(normal.y()!= Coord_type(0))
	return (Comparison_result) CGAL_NTS
	  sign(Vector(Coord_type(0),normal.z(),-normal.y())*(p-q));
      
      CGAL_assertion(normal.z()!= Coord_type(0));
      return (Comparison_result) CGAL_NTS 
	sign(Vector(Coord_type(0),-normal.z(),normal.y())*(p-q));
     }
  
private:
  Vector normal;
};

template <class R>
class Compute_area_3
{
  //squareroot of compute_squared_area_3:
  //-> if no sqrt is supported, cast to double
public:
  typedef typename  R::FT                     FT;
  typedef typename  R::Triangle_3             Triangle;
  typedef typename  R::Point_3                Point;
  typedef typename  R::Compute_squared_area_3 Compute_squared_area_3;

 //operators:
  FT operator() (const Triangle &t){
    return 
      cast_sqrt_to_double(t,
			  Number_type_traits<FT>::Has_sqrt);
  }
  FT operator() 
    (const Point& p,const Point& q,const Point& r){
    return 
      cast_sqrt_to_double(p,q,r,
			  Number_type_traits<FT>::Has_sqrt);
  }

private:
  //sqrt is supported by FT:
  inline FT 
  cast_sqrt_to_double(const Triangle &t,Tag_true)
    {
      return CGAL_NTS sqrt(Compute_squared_area_3()(t));
    }
  inline FT cast_sqrt_to_double
  (const Point& p,const Point& q,const Point& r, Tag_true)
    {
      return CGAL_NTS sqrt(Compute_squared_area_3()
			   (Triangle(p,q,r)));
    }
  //no sqrt supported:
  inline FT 
  cast_sqrt_to_double(const Triangle &t,Tag_false)
    {
      double approx = 
	to_double(Compute_squared_area_3()(t));
      return CGAL_NTS sqrt(approx);
    }
  inline FT cast_sqrt_to_double
  (const Point& p,const Point& q,const Point& r, Tag_false)
    {
      double approx = 
	to_double(Compute_squared_area_3()(Triangle(p,q,r)));
      return CGAL_NTS sqrt(approx);
    }

};

template < class _R>
class Voronoi_intersection_2_traits_3
{
public:
  typedef _R Rep;

  typedef typename Rep::RT                          Weight;
  typedef typename Rep::FT                          FT;

  //other types needed:
  typedef typename Rep::Point_3                     Point_2;
  typedef typename Rep::Segment_3                   Segment_2;
  typedef typename Rep::Triangle_3                  Triangle_2;
  typedef typename Rep::Line_3                      Line_2;
  typedef typename Rep::Ray_3                       Ray_2;
  typedef typename Rep::Direction_3                 Direction_2;
  typedef typename Rep::Vector_3                    Vector_2;
  typedef typename Rep::Construct_ray_3             Construct_ray_2;
  typedef typename Rep::Construct_direction_3
  Construct_direction_2;
  typedef typename Rep::Construct_triangle_3      Construct_triangle_2;

  typedef typename Rep::Compare_distance_3          Compare_distance_2;
  //if no sqrt is supported, it casts to double:
  typedef Compute_area_3<Rep>                       Compute_area_2;
  
  //the regular triangulation traits model:
  //Traits::Point_2 is a 3D point!!
  typedef  Point_2                        Weighted_point;
  typedef  Point_2                        Bare_point;
  
  
  //specific tests:
  typedef Orientation_with_normal_plane_2_3<Rep>          Orientation_2;
  typedef Side_of_plane_centered_sphere_2_3<Point_2>       Power_test_2;
  typedef Side_of_plane_centered_sphere_degenerated_2_3<Point_2>  
  Power_test_degenerated_2;
  
  typedef Construct_plane_centered_circumcenter_3<Point_2> 
  Construct_weighted_circumcenter_2;
  
  typedef Construct_plane_intersected_bisector_3<Point_2>    
  Construct_radical_axis_2;
  
  typedef Compare_first_projection_3<Point_2>              Compare_x_2;
  typedef Compare_second_projection_3<Point_2>             Compare_y_2;
  
  
  
  
  //instantiations and creation of functors:
  //for the triangulation:
  Orientation_2 
  orientation_2_object() const { 
    return Orientation_2(normal, Rep());}
  
  Power_test_2 
  power_test_2_object() const
    { return Power_test_2(a,normal);}
  
  Power_test_degenerated_2
  power_test_degenerated_2_object() const
    { return Power_test_degenerated_2(a,normal);}
  
  Compare_distance_2 compare_distance_2_object() const
    { return Compare_distance_2(); }
  
  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2(normal);}
  
  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2(normal);}

  //for the coordinate computation:
  Compute_area_2 compute_area_2_object() const
    { return Compute_area_2();}
  
  //for constructions of dual:
  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object() const
    {
      return Construct_weighted_circumcenter_2(a,normal);
    }
  
  Construct_radical_axis_2  construct_radical_axis_2_object() const
    {return Construct_radical_axis_2(a,normal);}
  
  Construct_ray_2  construct_ray_2_object() const
    {return Construct_ray_2();}
  
  Construct_direction_2  construct_direction_2_object() const
    {return Construct_direction_2();}

  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

  
public:  
    Voronoi_intersection_2_traits_3<_R>(const Point_2& _p = Point_2(), 
					const Vector_2& _normal = NULL_VECTOR) 
      : a(_p), normal(_normal){};
  
  Vector_2 get_normal() const {return normal;}
  Point_2  get_point() const {return a;}
  
private:
  //defining the intersection plane:
  Point_2 a;
  Vector_2 normal; 
};


 /********************************************/
//put the homogeneous or cartesian tag 
template <  class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point  &q,
           const Point  &r,
           const Point  &t)
{
  typedef typename Point::R::Rep_tag Tag;

  return( side_of_plane_centered_sphere(a,n,p,q,r,t, Tag()) );
};

template <   class Point, class Vector>
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point &q,
           const Point &r)
{
  typedef typename Point::R::Rep_tag Tag;
  return( side_of_plane_centered_sphere(a,n,p,q,r,Tag()) );
};

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point &q)
{
  typedef typename Point::R::RT  RT;
  Comparison_result r = Compare<RT>()(-CGAL_NTS square(Vector(a,p)*n), 
				      -CGAL_NTS
				      square(Vector(a,q)*n));
  
  if(r == LARGER)    return ON_NEGATIVE_SIDE;
  else if (r == SMALLER) return ON_POSITIVE_SIDE;
  return ON_ORIENTED_BOUNDARY;
};



template <  class Point, class Vector >
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n, 
			      /*defines the plane*/
			      const Point &p,
			      const Point  &q,
			      const Point  &r)
{ 
  typedef typename Point::R::Rep_tag Tag;

  
  return( plane_centered_circumcenter_3(a,n,p,q,r, Tag()));
};

template <  class Point, class Vector>
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n, 
			     /*defines the plane*/
			     const Point &p,
			     const Point  &q)
{ 
  typedef typename Point::R::Rep_tag Tag;
  
  return( plane_intersected_bisector_3(a,n,p,q, Tag()));
};

///-----------------------------------------------------------
//#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
template < class Point, class Vector>
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point  &q,
           const Point  &r,
           const Point  &t,
	   Cartesian_tag)
{
  return side_of_plane_centered_sphereC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z(), 
     t.x(), t.y(), t.z());
};

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point &q,
           const Point &r,
	   Cartesian_tag)
{
  return side_of_plane_centered_sphereC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z());
};


template < class Point, class Vector>
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n, 
			      /*defines the plane*/
			      const Point &p,
			      const Point  &q,
			      const Point  &r,
			      Cartesian_tag)
{
  typename Point::R::RT x,y,z;
  plane_centered_circumcenterC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z(),x,y,z);
  return (Point(x,y,z));
};

template < class Point,  class Vector>
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n, 
			     /*defines the plane*/
			     const Point &p,
			     const Point  &q,
			     Cartesian_tag)
{
  typename Point::R::RT x1,y1,z1, x2,y2,z2;
  typedef typename Point::R::Line_3 Line;
  bisector_plane_intersectionC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(), x1,y1,z1,x2,y2,z2);

  return Line(Point(x1,y1,z1), Point(x2,y2,z2));
};

//#endif // CGAL_CARTESIAN_H
//#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H

// The 3 following call the cartesian version over FT, because an 
// homogeneous special version has not yet been written.

template <class Point,  class Vector  >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point  &q,
           const Point  &r,
           const Point  &t,
	   Homogeneous_tag)
{
  
    return side_of_plane_centered_sphereC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z(), 
     t.x(), t.y(), t.z());
};

template < class Point, class Vector >
inline
Oriented_side
side_of_plane_centered_sphere(const Point &a, const Vector& n, 
/*defines the plane*/
	   const Point &p,
           const Point &q,
           const Point &r,
	   Homogeneous_tag)
{
  return side_of_plane_centered_sphereC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z());
};


template < class Point, class Vector>
inline
Point
plane_centered_circumcenter_3(const Point &a, const Vector& n, 
			      /*defines the plane*/
			      const Point &p,
			      const Point  &q,
			      const Point  &r,
			      Homogeneous_tag)
{
  typename Point::R::RT x,y,z;
  plane_centered_circumcenterC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),
     r.x(), r.y(), r.z(),x,y,z);
  return Point(x,y,z);
};


template < class Point,  class Vector>
inline
typename Point::R::Line_3
plane_intersected_bisector_3(const Point &a, const Vector& n, 
			     /*defines the plane*/
			     const Point &p,
			     const Point  &q,
			     Homogeneous_tag)
{

  typename Point::R::RT x1,y1,z1, x2,y2,z2;
  typedef typename Point::R::Line_3 Line;
  bisector_plane_intersectionC3
    (a.x(), a.y(), a.z(),
     n.x(), n.y(), n.z(),
     p.x(), p.y(), p.z(),
     q.x(), q.y(), q.z(),x1,y1,z1,x2,y2,z2);
  
  return Line(Point(x1,y1,z1), Point(x2,y2,z2));
};
//#endif // CGAL_HOMOGENEOUS_H



CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_INTERSECTION_2_TRAITS_3_H
/// /ATTENTION: cast the squareroot to double:
// template <class R>
// class Compute_area_3_double
// {
//   //squareroot of compute_squared_area_3 IN DOUBLE:
// public:
//   typedef typename  R::Triangle_3             Triangle;
//   typedef typename  R::Point_3                Point;
//   typedef typename  R::Compute_squared_area_3 Compute_squared_area_3;

//   typename R::FT operator() (const Triangle &t){
//     double approx = to_double(Compute_squared_area_3()(t));
 
//     return CGAL_NTS sqrt(approx);
//   }
//   typename R::FT operator() (const Point& p,const Point& q,
// 			     const Point& r){
   
//     double approx = 
//       to_double(Compute_squared_area_3()(Triangle(p,q,r)));
    
//     return CGAL_NTS sqrt(approx);
//   }
// };

// template <class R>
// class Compute_area_3
// {
//   //squareroot of compute_squared_area_3:
// public:
//   typedef typename  R::Triangle_3             Triangle;
//   typedef typename  R::Point_3                Point;
//   typedef typename  R::Compute_squared_area_3 Compute_squared_area_3;

//   typename R::FT operator() (const Triangle &t){
//     return CGAL_NTS sqrt(Compute_squared_area_3()(t));
//   }
  
//    typename R::FT operator() (const Point& p,const Point& q,const Point& r){
//      return CGAL_NTS sqrt(Compute_squared_area_3()(Triangle(p,q,r)));
//    }
// };
 //no sqrt supported:
  //typedef Compute_area_3_double<Rep>                Compute_area_2;
  //if sqrt is supported:
  //typedef Compute_area_3<Rep>                       Compute_area_2;
