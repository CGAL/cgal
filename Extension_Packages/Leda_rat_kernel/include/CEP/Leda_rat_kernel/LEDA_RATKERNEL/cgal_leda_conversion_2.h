// conversion from CGAL to LEDA kernel (2d) and other direction ...

#ifndef CEP_LEDA_RAT_CONVERSIONS_2_H
#define CEP_LEDA_RAT_CONVERSIONS_2_H

#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>

// LEDA types ...
#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_line.h>
#include <LEDA/rat_circle.h>
#include <LEDA/rat_ray.h>
#include <LEDA/rat_triangle.h>
#include <LEDA/rat_rectangle.h>
#include <LEDA/rat_vector.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/rat_direction.h>

#if defined(CGAL_COMPATIBLE_CIRCLES)
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_rat_circle.h>
#endif


CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------------------
//  ... to CGAL:
//-----------------------------------------------------------------------------------


// attention : add special segment here ...

struct leda_to_cgal_2 {

  typedef CGAL::Homogeneous<leda_integer> HELP_KERNEL;
  
  HELP_KERNEL::Point_2 operator()(const leda_rat_point& p) const
  { return HELP_KERNEL::Point_2(p.X(), p.Y(), p.W()); }

  HELP_KERNEL::Segment_2 operator()(const leda_rat_segment& s) const
  {    
    HELP_KERNEL::Point_2 p1 = this->operator()(s.source());
    HELP_KERNEL::Point_2 p2 = this->operator()(s.target());
    return HELP_KERNEL::Segment_2(p1,p2); 
  }

  HELP_KERNEL::Line_2 operator()(const leda_rat_line& l) const
  {     
    HELP_KERNEL::Point_2 p1 = this->operator()(l.point1());
    HELP_KERNEL::Point_2 p2 = this->operator()(l.point2());
    return HELP_KERNEL::Line_2(p1,p2); 
  }

  HELP_KERNEL::Ray_2 operator()(const leda_rat_ray& r) const
  {    
    HELP_KERNEL::Point_2 p1 = this->operator()(r.point1());
    HELP_KERNEL::Point_2 p2 = this->operator()(r.point2());
    return HELP_KERNEL::Ray_2(p1,p2); 
  }

  HELP_KERNEL::Triangle_2 operator()(const leda_rat_triangle& t) const
  {   
    HELP_KERNEL::Point_2 p1 = this->operator()(t.point1());
    HELP_KERNEL::Point_2 p2 = this->operator()(t.point2());
    HELP_KERNEL::Point_2 p3 = this->operator()(t.point3());
    
    return HELP_KERNEL::Triangle_2(p1,p2,p3); 
  }

  HELP_KERNEL::Vector_2 operator()(const leda_rat_vector& v) const
  {     
    return HELP_KERNEL::Vector_2(v.X(),v.Y(),v.W()); 
  }

  HELP_KERNEL::Direction_2 operator()(const LEDA_NAMESPACE_NAME::rat_direction& d) const
  {       
    HELP_KERNEL::Vector_2 v = this->operator()(d.get_vector());
  
    return HELP_KERNEL::Direction_2(v); 
  }

  HELP_KERNEL::Iso_rectangle_2 operator()(const leda_rat_rectangle& r) const
  {    
    HELP_KERNEL::Point_2 p1 = this->operator()(r.lower_left());
    HELP_KERNEL::Point_2 p2 = this->operator()(r.upper_right());
    return HELP_KERNEL::Iso_rectangle_2(p1,p2); 
  }

  HELP_KERNEL::Circle_2 operator()(const leda_rat_circle& c) const
  {    
    HELP_KERNEL::Point_2 p1 = this->operator()(c.point1());
    HELP_KERNEL::Point_2 p2 = this->operator()(c.point2());
    HELP_KERNEL::Point_2 p3 = this->operator()(c.point3());
    
    return HELP_KERNEL::Circle_2(p1,p2,p3); 
  }

#if defined(CGAL_COMPATIBLE_CIRCLES)
  HELP_KERNEL::Circle_2 operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c) const
  {    
    HELP_KERNEL::Point_2 center = this->operator()(c.center());
    leda_rational        rad = c.sqr_radius();
    
    return HELP_KERNEL::Circle_2(center, CGAL::Quotient<leda_integer>(rad.numerator(),rad.denominator())); 
  }
#endif

  // object conversion; the Object stores a leda 2d object ...
  CGAL::Object  operator()(const CGAL::Object& leda_obj) const
  {
    leda_rat_point p;
    
    if (CGAL::assign(p, leda_obj)){
      return CGAL::make_object(this->operator()(p));
    }
    
    leda_rat_segment s;
    
    if (CGAL::assign(s, leda_obj)){
      return CGAL::make_object(this->operator()(s));
    }    
    
    leda_rat_line l;
    
    if (CGAL::assign(l, leda_obj)){
      return CGAL::make_object(this->operator()(l));
    }    
    
    leda_rat_ray r;
    
    if (CGAL::assign(r, leda_obj)){
      return CGAL::make_object(this->operator()(r));
    }    
    
    leda_rat_triangle t;
    
    if (CGAL::assign(t, leda_obj)){
      return CGAL::make_object(this->operator()(t));
    }    
    
    leda_rat_vector v;
    
    if (CGAL::assign(v, leda_obj)){
      return CGAL::make_object(this->operator()(v));
    }    
    
    LEDA_NAMESPACE_NAME::rat_direction d;
    
    if (CGAL::assign(d, leda_obj)){
      return CGAL::make_object(this->operator()(d));
    }    
    
    leda_rat_rectangle rect;
    
    if (CGAL::assign(rect, leda_obj)){
      return CGAL::make_object(this->operator()(rect));
    }    
    
    leda_rat_circle c;
    
    if (CGAL::assign(c, leda_obj)){
      return CGAL::make_object(this->operator()(c));
    }    
    
#if defined(CGAL_COMPATIBLE_CIRCLES)
    LEDA_NAMESPACE_NAME::cgal_rat_circle circ;
    
    if (CGAL::assign(circ, leda_obj)){
      return CGAL::make_object(this->operator()(circ));
    }    
#endif    

    // this might happen for empty objects ...
    
    return leda_obj;
  }
};

//-----------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------
//  ... to LEDA:
//-----------------------------------------------------------------------------------

struct cgal_to_leda_2 {

  typedef CGAL::Homogeneous<leda_integer> HELP_KERNEL;

  leda_rat_point operator()(const HELP_KERNEL::Point_2& p) const
  { return leda_rat_point(p.hx(), p.hy(), p.hw()); }

  leda_rat_segment operator()(const HELP_KERNEL::Segment_2& s) const
  {    
    leda_rat_point p1 = this->operator()(s.source());
    leda_rat_point p2 = this->operator()(s.target());
    return leda_rat_segment(p1,p2); 
  }

  leda_rat_line operator()(const HELP_KERNEL::Line_2& l) const
  {    
    leda_rat_point p1 = this->operator()(l.point(0));
    leda_rat_point p2 = this->operator()(l.point(1));
    return leda_rat_line(p1,p2); 
  }

  leda_rat_ray operator()(const HELP_KERNEL::Ray_2& r) const
  {    
    leda_rat_point p1 = this->operator()(r.point(0));
    leda_rat_point p2 = this->operator()(r.point(1));
    return leda_rat_ray(p1,p2); 
  }

  leda_rat_triangle operator()(const HELP_KERNEL::Triangle_2& t) const
  {   
    leda_rat_point p1 = this->operator()(t.vertex(0));
    leda_rat_point p2 = this->operator()(t.vertex(1));
    leda_rat_point p3 = this->operator()(t.vertex(2));
    
    return leda_rat_triangle(p1,p2,p3); 
  }

  leda_rat_vector operator()(const HELP_KERNEL::Vector_2& v) const
  {     
    return leda_rat_vector(v.hx(),v.hy(),v.hw()); 
  }

  LEDA_NAMESPACE_NAME::rat_direction operator()(const HELP_KERNEL::Direction_2& d) const
  {       
    leda_rat_vector v = this->operator()(d.vector());
  
    return LEDA_NAMESPACE_NAME::rat_direction(v); 
  }

  leda_rat_rectangle operator()(const HELP_KERNEL::Iso_rectangle_2& r) const
  {    
    leda_rat_point p1 = this->operator()(r.min());
    leda_rat_point p2 = this->operator()(r.max());
    return leda_rat_rectangle(p1,p2); 
  }

  // conversion to leda circles does not work ...

#if defined(CGAL_COMPATIBLE_CIRCLES)
  LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const HELP_KERNEL::Circle_2& c) const
  {    
    leda_rat_point center = this->operator()(c.center());
    CGAL::Quotient<leda_integer>  q = c.squared_radius();
    leda_rational  rad(q.numerator(), q.denominator() );
    
    return LEDA_NAMESPACE_NAME::cgal_rat_circle(center, rad); 
  }
#endif

  // object conversion; the Object stores a cgal 2d object ...
  CGAL::Object  operator()(const CGAL::Object& cgal_obj) const
  {
    HELP_KERNEL::Point_2 p;
    
    if (CGAL::assign(p, cgal_obj)){
      return CGAL::make_object(this->operator()(p));
    }
    
    HELP_KERNEL::Segment_2 s;
    
    if (CGAL::assign(s, cgal_obj)){
      return CGAL::make_object(this->operator()(s));
    }    
    
    HELP_KERNEL::Line_2 l;
    
    if (CGAL::assign(l, cgal_obj)){
      return CGAL::make_object(this->operator()(l));
    }    
    
    HELP_KERNEL::Ray_2 r;
    
    if (CGAL::assign(r, cgal_obj)){
      return CGAL::make_object(this->operator()(r));
    }    
    
    HELP_KERNEL::Triangle_2 t;
    
    if (CGAL::assign(t, cgal_obj)){
      return CGAL::make_object(this->operator()(t));
    }    
    
    HELP_KERNEL::Vector_2 v;
    
    if (CGAL::assign(v, cgal_obj)){
      return CGAL::make_object(this->operator()(v));
    }    
    
    HELP_KERNEL::Direction_2 d;
    
    if (CGAL::assign(d, cgal_obj)){
      return CGAL::make_object(this->operator()(d));
    }    
    
    HELP_KERNEL::Iso_rectangle_2 rect;
    
    if (CGAL::assign(rect, cgal_obj)){
      return CGAL::make_object(this->operator()(rect));
    }    
    
    // conversion to LEDA circles will not work ...
    
#if defined(CGAL_COMPATIBLE_CIRCLES)
    HELP_KERNEL::Circle_2 circ;
    
    if (CGAL::assign(circ, cgal_obj)){
      return CGAL::make_object(this->operator()(circ));
    }    
#endif    

    // this might happen for empty objects ...
    
    return cgal_obj;  
  }
  
};


//-----------------------------------------------------------------------------------

CGAL_END_NAMESPACE

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif


#endif
