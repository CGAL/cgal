// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, December 15
//
// file          : include/CGAL/Triangulation_euclidean_leda_traits_2.h
// package       : Point_set_2_tb (0.1)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.1
// revision_date : 15 December 2000 
// author(s)     : Matthias Baesken
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================

// LEDA traits classes for CGAL Delaunay triangulation

#ifndef CGAL_TRIANGULATION_EUCLIDEAN_LEDA_TRAITS_2_H
#define CGAL_TRIANGULATION_EUCLIDEAN_LEDA_TRAITS_2_H

// float kernel
#include <LEDA/point.h>
#include <LEDA/segment.h>
#include <LEDA/triangle.h>
#include <LEDA/circle.h>
#include <LEDA/line.h>
#include <LEDA/ray.h>
#include <LEDA/vector.h>

// rat kernel
#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_triangle.h>
#include <LEDA/rat_circle.h>
#include <LEDA/rat_line.h>
#include <LEDA/rat_ray.h>
#include <LEDA/rat_vector.h>


#include <CGAL/assertions.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>


//#define LEDA_TRIANG_DEBUG


template<class I>
class leda_distance {
public:  
  typedef leda_point Point;
  
  leda_distance(const I* = NULL) {}
  
  leda_distance(const Point& p0,const I* = NULL)
  { p[0]=p0; }

  leda_distance(const Point& p0, const Point& p1,const I* = NULL)
  { p[0]=p0; p[1]=p1; }


  leda_distance(const Point& p0, const Point& p1, const Point& p2,const I* = NULL)
  { p[0]=p0; p[1]=p1; p[2]=p2; }

  void set_point(int i, const Point& q)
  {
    CGAL_precondition( ((unsigned int) i) < 3 );
    p[i] = q;
  }

  Point get_point(int i) const
  {
    CGAL_precondition( ((unsigned int) i) < 3 );
    return p[i];
  }

  CGAL::Comparison_result compare() const
  {
    return (CGAL::Comparison_result) (::cmp_distances(p[0], p[1], p[0], p[2]));
  }

private:
  Point p[3];
};

template<class I>
class leda_rat_distance {
public:
  typedef leda_rat_point Point;

  leda_rat_distance(const I* = NULL) {}

  leda_rat_distance(const Point& p0,const I* = NULL)
  { p[0]=p0; }

  leda_rat_distance(const Point& p0, const Point& p1,const I* = NULL)
  { p[0]=p0; p[1]=p1; }

  leda_rat_distance(const Point& p0, const Point& p1, const Point& p2,const I* = NULL)
  { p[0]=p0; p[1]=p1; p[2]=p2; }

  void set_point(int i, const Point& q)
  {
    CGAL_precondition( ((unsigned int) i) < 3 );
    p[i] = q;
  }

  Point get_point(int i) const
  {
    CGAL_precondition( ((unsigned int) i) < 3 );
    return p[i];
  }

  CGAL::Comparison_result compare() const
  {
    return (CGAL::Comparison_result) (::cmp_distances(p[0], p[1], p[0], p[2]));
  }
private:
  Point p[3];
};


CGAL_BEGIN_NAMESPACE 

// rational kernel ...

class Leda_rat_compare_x_2 {
public:
    Comparison_result operator()(const leda_rat_point &p, const leda_rat_point &q) const
    {
        return (CGAL::Comparison_result)(leda_rat_point::cmp_x(p, q));
    }
};

class Leda_rat_compare_y_2 {
public:
    Comparison_result operator()(const leda_rat_point &p, const leda_rat_point &q) const
    {
        return (CGAL::Comparison_result)(leda_rat_point::cmp_y(p, q));        
    }
};

class Leda_rat_orientation_2 {
public:
    Orientation operator()(const leda_rat_point &p,const leda_rat_point &q, const leda_rat_point &r) const
    {
        return (CGAL::Orientation)(::orientation(p, q, r));
    }
};

class Leda_rat_side_of_oriented_circle_2 {
public:
    Oriented_side operator()(const leda_rat_point &p,const leda_rat_point &q,
					const leda_rat_point &r,
					const leda_rat_point &s) const
    {
      return (CGAL::Oriented_side)(::side_of_circle(p, q, r, s));
    }
};

class Leda_rat_construct_circumcenter_2 {
public:
    leda_rat_point operator()(const leda_rat_point &p, const leda_rat_point &q, const leda_rat_point &r) const
    {
        leda_rat_circle C(p,q,r);
        return C.center();
    }
};

class Leda_rat_construct_bisector_2 {
public:
    leda_rat_line operator()(const leda_rat_point &p1, leda_rat_point& p2) const
    {
      return leda_rat_line(::p_bisector(p1,p2));
    }
};

class Leda_rat_construct_midpoint {
public:
    leda_rat_point operator()(const leda_rat_point &p, const leda_rat_point &q) const
    {
        return ::midpoint(p,q);        
    }
};

class Leda_rat_construct_segment_2 {
public:
    leda_rat_segment operator()(const leda_rat_point &p, const leda_rat_point &q) const
    {
        return leda_rat_segment(p,q);        
    }
};

class Leda_rat_construct_triangle_2 {
public:
    leda_rat_triangle operator()(const leda_rat_point &p, const leda_rat_point &q, const leda_rat_point &r) const
    {
        return leda_rat_triangle(p,q,r);        
    }
};

class Leda_rat_less_distance_to_point_2 {
    leda_rat_point _p;
public:
    Leda_rat_less_distance_to_point_2( const leda_rat_point& p) : _p(p) {}

    bool operator()( const leda_rat_point& p1, const leda_rat_point& p2) const
    {
      int c = ::cmp_distances(_p,p1,_p,p2);
      if (c==-1) return true; else return false;
    }
};

class Leda_rat_construct_direction_of_line_2  {
public:
    leda_rat_vector  operator()( const leda_rat_line& l)
    { return (l.point2() - l.point1()); }
};

class Leda_rat_construct_ray_2 {
public:
    leda_rat_ray operator()(const leda_rat_point& p, const leda_rat_vector& d)
    { return leda_rat_ray(p, p+d); }
}; 


// float kernel ...

class Leda_float_compare_x_2 {
public:
    Comparison_result operator()(const leda_point &p, const leda_point &q) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "compare_x_2!\n"; std::cout.flush();
#endif
        return (CGAL::Comparison_result)(leda_point::cmp_x(p, q));
    }
};

class Leda_float_compare_y_2 {
public:
    Comparison_result operator()(const leda_point &p, const leda_point &q) const
    {
#if defined LEDA_TRIANG_DEBUG 
        std::cout << "compare_y_2!\n"; std::cout.flush();   
#endif    
        return (CGAL::Comparison_result)(leda_point::cmp_y(p, q));        
    }
};

class Leda_float_orientation_2 {
public:
    Orientation operator()(const leda_point &p,const leda_point &q, const leda_point &r) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "orientation_2!\n"; std::cout.flush();
#endif    
        return (CGAL::Orientation)(::orientation(p, q, r));
    }
};

class Leda_float_side_of_oriented_circle_2 {
public:
    Oriented_side operator()(const leda_point &p,const leda_point &q,
					const leda_point &r,
					const leda_point &s) const
    {
#if defined LEDA_TRIANG_DEBUG 
        std::cout << "side_of_oriented_circle_2!\n"; std::cout.flush();  
	if (::collinear(p, q, r))  { std::cout << "  circle is degenerate!\n"; std::cout.flush(); }
#endif    
      return (CGAL::Oriented_side)(::side_of_circle(p, q, r, s));
    }
};

class Leda_float_construct_circumcenter_2 {
public:
    leda_point operator()(const leda_point &p, const leda_point &q, const leda_point &r) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_circumcenter_2!\n"; std::cout.flush();
#endif    
        if (::collinear(p,q,r)){
	  std::cout << "collinear!\n"; std::cout.flush();
	  return leda_point();
        }
	else {
         leda_circle C(p,q,r);
         return C.center();
	} 
    }
};

class Leda_float_construct_bisector_2 {
public:
    leda_line operator()(const leda_point& p1, const leda_point& p2 ) const
    {
#if defined LEDA_TRIANG_DEBUG   
        std::cout << "construct_bisector_2!\n"; std::cout.flush(); 
#endif    
      return leda_line(::p_bisector(p1,p2));
    }
};

class Leda_float_construct_midpoint {
public:
    leda_point operator()(const leda_point &p, const leda_point &q) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_midpoint!\n"; std::cout.flush();
#endif    
        return ::midpoint(p,q);        
    }
};

class Leda_float_construct_segment_2 {
public:
    leda_segment operator()(const leda_point &p, const leda_point &q) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_segment_2!\n"; std::cout.flush();
#endif    
        return leda_segment(p,q);        
    }
};

class Leda_float_construct_triangle_2 {
public:
    leda_triangle operator()(const leda_point &p, const leda_point &q, const leda_point &r) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_triangle_2!\n"; std::cout.flush();
#endif    
        return leda_triangle(p,q,r);        
    }
};

class Leda_float_less_distance_to_point_2 {
    leda_point _p;
public:
    Leda_float_less_distance_to_point_2( const leda_point& p) : _p(p) {}

    bool operator()( const leda_point& p1, const leda_point& p2) const
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "less_distance_to_point_2!\n"; std::cout.flush();
#endif    
      int c = ::cmp_distances(_p,p1,_p,p2);
      if (c==-1) return true; else return false;
    }
};

class Leda_float_construct_direction_of_line_2  {
public:
    leda_vector  operator()( const leda_line& l)
    { 
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_direction_of_line_2!\n"; std::cout.flush();
#endif    
      return (l.point2() - l.point1()); 
    }
};

class Leda_float_construct_ray_2 {
public:
    leda_ray operator()(const leda_point& p, const leda_vector& d)
    {
#if defined LEDA_TRIANG_DEBUG    
        std::cout << "construct_ray_2!\n"; std::cout.flush();
#endif     
      return leda_ray(p,d); 
    }
}; 


class Triangulation_euclidean_leda_rat_traits_2 {
public:
  typedef leda_rational Rep;
  typedef leda_rat_point        Point_2;
  typedef leda_rat_segment      Segment_2;
  typedef leda_rat_triangle     Triangle_2;
  typedef leda_rat_circle       Circle_2;
  typedef leda_rat_line         Line_2;
  typedef leda_rat_ray          Ray_2;
  typedef leda_rat_vector       Direction_2;
  typedef leda_rat_distance<Triangulation_euclidean_leda_rat_traits_2> Distance;
  
  // new
  typedef Leda_rat_compare_x_2                Compare_x_2;
  typedef Leda_rat_compare_y_2                Compare_y_2;
  typedef Leda_rat_orientation_2              Orientation_2;
  typedef Leda_rat_side_of_oriented_circle_2  Side_of_oriented_circle_2;
  typedef Leda_rat_construct_circumcenter_2   Construct_circumcenter_2;
  typedef Leda_rat_construct_bisector_2       Construct_bisector_2;
  typedef Leda_rat_construct_segment_2        Construct_segment_2;
  typedef Leda_rat_construct_triangle_2       Construct_triangle_2;
  typedef Leda_rat_construct_midpoint         Construct_midpoint;
  typedef Leda_rat_less_distance_to_point_2   Less_distance_to_point_2;   
  
  typedef Leda_rat_construct_direction_of_line_2  Construct_direction_of_line_2;
  typedef Leda_rat_construct_ray_2                Construct_ray_2;
  

  Triangulation_euclidean_leda_rat_traits_2() {}
  Triangulation_euclidean_leda_rat_traits_2(const Triangulation_euclidean_leda_rat_traits_2 &) {}
  Triangulation_euclidean_leda_rat_traits_2 &operator=
      (const Triangulation_euclidean_leda_rat_traits_2 &)
  {return *this;}
 
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
 
  Construct_circumcenter_2
  construct_circumcenter_2_object() const
    { return Construct_circumcenter_2();}

  Construct_bisector_2
  construct_bisector_2_object() const
    {return Construct_bisector_2();}
  
  Construct_midpoint
  construct_midpoint_object() const
    {return Construct_midpoint();}

  Construct_segment_2
  construct_segment_2_object() const
  {return Construct_segment_2(); }

  Construct_triangle_2
  construct_triangle_2_object() const
  {return Construct_triangle_2(); }

  Less_distance_to_point_2
  less_distance_to_point_2_object(const Point_2& p) const
    {return Less_distance_to_point_2(p);}
    
  Construct_direction_of_line_2
  construct_direction_of_line_2_object() const
  { return  Construct_direction_of_line_2(); }
  
  Construct_ray_2
  construct_ray_2_object() const
  { return  Construct_ray_2(); }
};


class Triangulation_euclidean_leda_float_traits_2 {
public:
  typedef double Rep;
  typedef leda_point      Point_2;
  typedef leda_segment    Segment_2;
  typedef leda_triangle   Triangle_2;
  typedef leda_circle     Circle_2;
  typedef leda_line       Line_2;
  typedef leda_ray        Ray_2;
  typedef leda_vector     Direction_2;
  typedef leda_distance<Triangulation_euclidean_leda_float_traits_2> Distance;

  // new
  typedef Leda_float_compare_x_2                Compare_x_2;
  typedef Leda_float_compare_y_2                Compare_y_2;
  typedef Leda_float_orientation_2              Orientation_2;
  typedef Leda_float_side_of_oriented_circle_2  Side_of_oriented_circle_2;
  typedef Leda_float_construct_circumcenter_2   Construct_circumcenter_2;
  typedef Leda_float_construct_bisector_2       Construct_bisector_2;
  typedef Leda_float_construct_segment_2        Construct_segment_2;
  typedef Leda_float_construct_triangle_2       Construct_triangle_2;  
  typedef Leda_float_construct_midpoint         Construct_midpoint;
  typedef Leda_float_less_distance_to_point_2   Less_distance_to_point_2; 
  
  typedef Leda_float_construct_direction_of_line_2  Construct_direction_of_line_2;
  typedef Leda_float_construct_ray_2                Construct_ray_2;    
  

  Triangulation_euclidean_leda_float_traits_2() {}
  Triangulation_euclidean_leda_float_traits_2(const Triangulation_euclidean_leda_float_traits_2 &) {}
  Triangulation_euclidean_leda_float_traits_2 &operator=
      (const Triangulation_euclidean_leda_float_traits_2 &)
  {return *this;}
 
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
 
  Construct_circumcenter_2
  construct_circumcenter_2_object() const
    { return Construct_circumcenter_2();}

  Construct_bisector_2
  construct_bisector_2_object() const
    {return Construct_bisector_2();}
  
  Construct_midpoint
  construct_midpoint_object() const
    {return Construct_midpoint();}

  Construct_segment_2
  construct_segment_2_object() const
  {return Construct_segment_2(); }

  Construct_triangle_2
  construct_triangle_2_object() const
  {return Construct_triangle_2(); }

  Less_distance_to_point_2
  less_distance_to_point_2_object(const Point_2& p) const
    {return Less_distance_to_point_2(p);}
    
  Construct_direction_of_line_2
  construct_direction_of_line_2_object() const
  { return  Construct_direction_of_line_2(); }
  
  Construct_ray_2
  construct_ray_2_object() const
  { return  Construct_ray_2(); }    
};


CGAL_END_NAMESPACE 

#endif // CGAL_TRIANGULATION_EUCLIDEAN_TRAITS_2_H
