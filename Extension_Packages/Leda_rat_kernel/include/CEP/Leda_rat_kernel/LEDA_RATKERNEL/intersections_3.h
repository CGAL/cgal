#ifndef CEP_LEDA_RAT_INTERSECTIONS_3_H
#define CEP_LEDA_RAT_INTERSECTIONS_3_H

// LEDA rational kernel intersection objects ...
// 3d intersections ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

#include <LEDA/d3_rat_plane.h>
#include <LEDA/d3_rat_line.h>
#include <LEDA/d3_rat_ray.h>
#include <LEDA/d3_rat_segment.h>

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE

class Assign_leda_rat_3 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

  template<class T>
  bool operator()(T& t, CGAL::Object& obj) const
  { 
    return CGAL::assign(t, obj);
  }
};

// computations of 3d intersections ...

template<class HELP_KERNEL>
class CGAL_intersect_leda_rat_3 {
public:
  typedef Arity_tag< 2 >     Arity;
  typedef CGAL::Object       result_type;
  
  typedef typename HELP_KERNEL::Intersect_3  Intersect_3;
  
  template<class T1, class T2>
  CGAL::Object operator()(const T1& obj1, const T2& obj2) const
  {
     leda_to_cgal_3 conv;
     Intersect_3 inter;
     
     // result "stores" a CGAL object ...
     CGAL::Object result = inter(conv(obj1), conv(obj2));
     
     // convert it back ...
     cgal_to_leda_3 conv_back;
     
     return conv_back(result);
  }  
};


/*
// computations of 2d intersections ...
class Intersect_leda_rat_3 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef CGAL::Object   result_type;

  CGAL::Object operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_plane& pl2) const
  {
    CGAL::Object obj;
    
    leda_d3_rat_point i1,i2;
    int res = pl.intersection(pl2, i1, i2);
    
    if (res == 0) return obj;  // no intersection ...
    if (res == 1) { // result is a line ...
      leda_d3_rat_line l(i1,i2); 
      return CGAL::make_object(l);
    }
    
    // result is whole plane ...
    return CGAL::make_object(pl);
  }

  CGAL::Object operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_line& l) const
  {
    CGAL::Object obj;
  
    leda_d3_rat_point inter;
    leda_d3_rat_point p1 = l.point1();
    leda_d3_rat_point p2 = l.point2();    
    
    int res = pl.intersection(p1,p2, inter);
    if (res == 0) return obj;  // no intersection ...    
    if (res == 1) { // point
       return CGAL::make_object(inter);
    }
    
    // line
    return CGAL::make_object(l);
  }  

  // -----------------------------------------------------------------------------------------------
  // first do intersection test;
  // if intersection is sure, compute the result using line intersection ...
  // -----------------------------------------------------------------------------------------------
  
  CGAL::Object operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_ray& r) const
  {
    CGAL::Object obj;
    leda_d3_rat_point p1 = r.point1(); // source
    leda_d3_rat_point p2 = r.point2(); // target
  
    int s1 = pl.side_of(p1);
    int s2 = pl.side_of(p2);
    
    if (s1==0 && s2==0) return CGAL::make_object(r); // whole ray is in the plane ...
    if (s1 != s2) return this->operator()(pl,leda_d3_rat_line(p1,p2));  
    
    // source and other point are on same side ...
    leda_rational r1 = pl.sqr_dist(p1);
    leda_rational r2 = pl.sqr_dist(p2); 
    
    if (r1 <= r2) return obj; // the ray is going away from the plane
    return this->operator()(pl,leda_d3_rat_line(p1,p2));    
  }  

  CGAL::Object operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_segment& s) const
  {
    CGAL::Object obj;
    int s1 = pl.side_of(s.source());
    int s2 = pl.side_of(s.target());
    
    if (s1==0 && s2==0) return CGAL::make_object(s);
    if (s1 != s2) return this->operator()(pl,leda_d3_rat_line(s));
    return obj; // empty intersection   
  }  

  // -----------------------------------------------------------------------------------------------
  
};
*/

CGAL_END_NAMESPACE


#endif




