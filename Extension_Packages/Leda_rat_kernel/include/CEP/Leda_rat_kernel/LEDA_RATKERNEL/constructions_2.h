#ifndef CEP_LEDA_RAT_CONSTRUCTIONS_2_H
#define CEP_LEDA_RAT_CONSTRUCTIONS_2_H

// LEDA rational kernel construction objects ...

// 2d constructions ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/support_functions.h>

#include <LEDA/rat_point.h>
#include <LEDA/rat_vector.h>
#include <LEDA/rat_line.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_triangle.h>
#include <LEDA/rat_rectangle.h>
#include <LEDA/rat_ray.h>
#include <LEDA/rat_circle.h>


#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

CGAL_BEGIN_NAMESPACE


class Construct_leda_rat_point {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_point           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_origin;
#endif   
  
  //----------------------------------------------------------------------------------
  //undocumented
  leda_rat_point operator()() const
  { 
   leda_rat_point p;
   return p;
  } 
  
  leda_rat_point operator()(const leda_rational& x, const leda_rational& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rational&,const leda_rational&>(Construct_leda_rat_point::ev_leda_rational, x, y);
#endif     
    leda_rat_point p(x,y);
    return p;
  }  
  
  leda_rat_point operator()(const leda_integer& x, const leda_integer& y, const leda_integer& w) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_integer&,const leda_integer&, const leda_integer&> \
       (Construct_leda_rat_point::ev_leda_integer, x, y, w);
#endif   
    leda_rat_point p(x,y,w);
    return p;
  } 
  
  //---------------------------------------------------------------------------------- 
  leda_rat_point operator()(const CGAL::Origin& orig) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const CGAL::Origin&> \
       (Construct_leda_rat_point::ev_origin, orig);
#endif    
    return leda_rat_point(0,0,1); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_point::ev_leda_rational;
CGAL::event Construct_leda_rat_point::ev_leda_integer;
CGAL::event Construct_leda_rat_point::ev_origin;
#endif   

class Construct_leda_rat_vector {
public:
  typedef leda_rat_vector           result_type;
  typedef Arity_tag< 2 > Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_leda_rat_point;
  static CGAL::event ev_null_vector;
#endif   

  //----------------------------------------------------------------------------------   
  //undocumented
  leda_rat_vector operator()() const
  { 
   leda_rat_vector v(2);
   return v;
  }  
  
  leda_rat_vector operator()(const leda_rational& x, const leda_rational& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rational&,const leda_rational&>(Construct_leda_rat_vector::ev_leda_rational, x, y);
#endif   
    leda_rat_vector v(x,y);
    return v;
  }  
  
  leda_rat_vector operator()(const leda_integer& x, const leda_integer& y, const leda_integer& w) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_integer&,const leda_integer&> \
      (Construct_leda_rat_vector::ev_leda_integer, x, y, w);
#endif    
    leda_rat_vector v(x,y,w);
    return v;
  } 
  
  //----------------------------------------------------------------------------------   

  leda_rat_vector operator()(const leda_rat_point& a, const leda_rat_point& b) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&,const leda_rat_point&> \
      (Construct_leda_rat_vector::ev_leda_rat_point, a, b);
#endif   
    return b-a; 
  }
   
  leda_rat_vector operator()(const CGAL::Null_vector& nv) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const CGAL::Null_vector&> \
      (Construct_leda_rat_vector::ev_null_vector, nv);
#endif     
    return leda_rat_vector(2); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_vector::ev_leda_rational;
CGAL::event Construct_leda_rat_vector::ev_leda_integer;
CGAL::event Construct_leda_rat_vector::ev_leda_rat_point;
CGAL::event Construct_leda_rat_vector::ev_null_vector;
#endif   

class Construct_leda_rat_direction {
public:
  typedef Arity_tag< 1 > Arity;
  typedef LEDA_NAMESPACE_NAME::rat_direction           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_rat_vector;
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;    
  static CGAL::event ev_leda_rat_segment;
#endif   
  
  //undocumented
  LEDA_NAMESPACE_NAME::rat_direction operator()() const
  { 
   LEDA_NAMESPACE_NAME::rat_direction d(2);
   return d;
  } 
  
  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rational& x, const leda_rational& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rational&,const leda_rational&> \
      (Construct_leda_rat_direction::ev_leda_rational, x, y);
#endif   
    LEDA_NAMESPACE_NAME::rat_direction d(x,y);
    return d;
  }        

  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rat_vector& v) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_vector&> \
      (Construct_leda_rat_direction::ev_leda_rat_vector, v);
#endif   
    return LEDA_NAMESPACE_NAME::rat_direction(v); 
  }
    
  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rat_line& l) const
  {   
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_line&> \
      (Construct_leda_rat_direction::ev_leda_rat_line, l);
#endif   
      leda_rat_vector v = l.point2()-l.point1();
      return LEDA_NAMESPACE_NAME::rat_direction(v); 
  }
    
  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rat_ray& r) const
  {  
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_ray&> \
      (Construct_leda_rat_direction::ev_leda_rat_ray, r);
#endif   
     leda_rat_vector v = r.point2()-r.point1(); 
     return LEDA_NAMESPACE_NAME::rat_direction(v);
  }
    
  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rat_segment& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_segment&>(Construct_leda_rat_direction::ev_leda_rat_segment, s);
#endif  
    leda_rat_vector v = s.end()-s.start();
    return LEDA_NAMESPACE_NAME::rat_direction(v); 
  }            
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event  Construct_leda_rat_direction::ev_leda_rational;
CGAL::event  Construct_leda_rat_direction::ev_leda_rat_vector;
CGAL::event  Construct_leda_rat_direction::ev_leda_rat_line;
CGAL::event  Construct_leda_rat_direction::ev_leda_rat_ray;    
CGAL::event  Construct_leda_rat_direction::ev_leda_rat_segment;
#endif


// attention - we need a special functor for special segments ...
// member template is here not possible ...

class Construct_leda_rat_segment {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_segment           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif    
  
  //undocumented
  leda_rat_segment operator()() const
  { 
   leda_rat_segment s;
   return s;
  }   

  leda_rat_segment operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_segment::ev_leda_rat_point, p1, p2);
#endif    
    return leda_rat_segment(p1,p2); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_segment::ev_leda_rat_point;
#endif   

class Construct_leda_rat_line {
public:
  typedef leda_rat_line           result_type;
  typedef Arity_tag< 2 >          Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_leda_rat_point; 
  static CGAL::event ev_leda_rat_point_direction; 
  static CGAL::event ev_leda_rat_segment;
  static CGAL::event ev_leda_rat_ray;      
#endif   
  
  //undocumented
  leda_rat_line operator()() const
  { 
   leda_rat_line l;
   return l;
  }    

  leda_rat_line operator()(const leda_integer& a, const leda_integer& b, const leda_integer& c) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_integer&, const leda_integer&, const leda_integer&> \
      (Construct_leda_rat_line::ev_leda_integer, a, b, c);
#endif   
  
      if (b == 0){ // par. to y - axis
       leda_rat_point p1(-c,1,a);
       leda_rat_point p2(-c,0,a);
       return leda_rat_line(p1,p2);      
      }
      if (a == 0){ // par. to x - axis
       leda_rat_point p1(0,-c,b);
       leda_rat_point p2(1,-c,b); 
       return leda_rat_line(p1,p2);        
      }
      // a == 0 and c == 0 not allowed
      
      leda_rat_point p1(0,-c,b);
      leda_rat_point p2(-c,0,a);
      
      return leda_rat_line(p1,p2);
  }

  leda_rat_line operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_line::ev_leda_rat_point, p1, p2);
#endif  
    return leda_rat_line(p1,p2); 
  }

  leda_rat_line operator()(const leda_rat_point& p1, const LEDA_NAMESPACE_NAME::rat_direction& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const LEDA_NAMESPACE_NAME::rat_direction&> \
      (Construct_leda_rat_line::ev_leda_rat_point_direction, p1, dir);
#endif  
    return leda_rat_line(p1, dir.get_vector()); 
  }

  leda_rat_line operator()(const leda_rat_segment& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_segment&> \
      (Construct_leda_rat_line::ev_leda_rat_segment, s);
#endif    
    return leda_rat_line(s.start(), s.end()); 
  }

  leda_rat_line operator()(const leda_rat_ray& r) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_ray&> \
      (Construct_leda_rat_line::ev_leda_rat_ray, r);
#endif     
    return leda_rat_line(r); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_line::ev_leda_integer;
CGAL::event Construct_leda_rat_line::ev_leda_rat_point; 
CGAL::event Construct_leda_rat_line::ev_leda_rat_point_direction; 
CGAL::event Construct_leda_rat_line::ev_leda_rat_segment;
CGAL::event Construct_leda_rat_line::ev_leda_rat_ray;      
#endif 


class Construct_leda_rat_ray {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_ray           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
  static CGAL::event ev_leda_rat_point_direction;      
#endif    
  
  //undocumented
  leda_rat_ray operator()() const
  { 
   leda_rat_ray r;
   return r;
  }    

  leda_rat_ray operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_ray::ev_leda_rat_point, p1, p2);
#endif        
    return leda_rat_ray(p1,p2); 
  }
    
  leda_rat_ray operator()(const leda_rat_point& p1, const LEDA_NAMESPACE_NAME::rat_direction& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const LEDA_NAMESPACE_NAME::rat_direction&> \
      (Construct_leda_rat_ray::ev_leda_rat_point_direction, p1, dir);
#endif    
    return leda_rat_ray(p1, dir.get_vector()); 
  }    
};


#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_ray::ev_leda_rat_point; 
CGAL::event Construct_leda_rat_ray::ev_leda_rat_point_direction;      
#endif 


#if defined(CGAL_COMPATIBLE_CIRCLES)
class Construct_leda_rat_circle {
public:
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle           result_type;
  typedef Arity_tag< 3 >          Arity;  

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_rational_orientation; 
  static CGAL::event ev_leda_rat_point_point_point; 
  static CGAL::event ev_leda_rat_point_point_orientation;
  static CGAL::event ev_leda_rat_point_orientation;       
#endif  
    
  //undocumented
  LEDA_NAMESPACE_NAME::cgal_rat_circle operator()() const
  { 
   LEDA_NAMESPACE_NAME::cgal_rat_circle  c;
   return c;
  }      
    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const leda_rat_point& center, 
                                                    const leda_rational&  sqrad, 
						    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rational&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_rational_orientation, center, sqrad, ori);
#endif    
      return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,sqrad,ori);
    }    

    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const leda_rat_point& p1, 
                                                    const leda_rat_point& p2, 
						    const leda_rat_point& p3) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_point, p1, p2, p3);
#endif    
     int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
     leda_rat_circle C(p1,p2,p3);
     leda_rat_point  center = C.center();
     
     CGAL::Orientation cg_ori;
     switch(ori){
        case -1: { cg_ori = CGAL::RIGHTTURN; break; }
	case  0: { cg_ori = CGAL::COLLINEAR; break; }
	case  1: { cg_ori = CGAL::LEFTTURN; break; }
     }
    
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,p1,cg_ori);
    }
    
    // diameter version ...    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const leda_rat_point& p1, const leda_rat_point& p2, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_orientation, p1, p2, ori);
#endif        
     leda_rat_point m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(m,p1,ori);
    }
    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const leda_rat_point& p1, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_orientation, p1, ori);
#endif     
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(p1,ori);
    }        
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_rational_orientation; 
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_point_point; 
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_point_orientation;
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_orientation;       
#endif 


#else
// use LEDA circles ...
// in this case we cannot provide some constructions ...

class Construct_leda_rat_circle {
public:
  typedef leda_rat_circle           result_type;
  typedef Arity_tag< 3 >          Arity;  

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_point_point; 
  static CGAL::event ev_leda_rat_point_point_orientation;
  static CGAL::event ev_leda_rat_point_orientation;       
#endif  
  
  //undocumented
  leda_rat_circle operator()() const
  { 
   leda_rat_circle  c;
   return c;
  }   

    leda_rat_circle operator()(const leda_rat_point& p1, const leda_rat_point& p2, const leda_rat_point& p3) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_point, p1, p2, p3);
#endif     
     return leda_rat_circle(p1,p2,p3);
    }
    
    leda_rat_circle operator()(const leda_rat_point& p1, const leda_rat_point& p2, CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_orientation, p1, p2, ori);
#endif    
     // p1p2 is the diameter ...
     leda_rat_point m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     leda_rat_circle C(m,p1);
     
     // check circle orientation ...
     int cori = (ori ==  CGAL::COUNTERCLOCKWISE) ?  -1 : 1;
     
     if (C.orientation() == cori) return C;
     else C = C.reverse();
     
     return C;
    }
    
    leda_rat_circle operator()(const leda_rat_point& p1, CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_orientation, p1, ori);
#endif     
     return leda_rat_circle(p1,p1);
    }        
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_point_point; 
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_point_orientation;
CGAL::event Construct_leda_rat_circle::ev_leda_rat_point_orientation;       
#endif

#endif

class Construct_leda_rat_triangle {
public:
  typedef Arity_tag< 3 > Arity;
  typedef leda_rat_triangle           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
#endif    

  //undocumented
  leda_rat_triangle operator()() const
  { 
   leda_rat_triangle  t;
   return t;
  }  

  leda_rat_triangle operator()(const leda_rat_point& p1, const leda_rat_point& p2, const leda_rat_point& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_triangle::ev_leda_rat_point, p1, p2, p3);
#endif   
   return leda_rat_triangle(p1,p2,p3);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_triangle::ev_leda_rat_point; 
#endif 

class Construct_leda_rat_rectangle {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_rectangle           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
#endif    

  //undocumented
  leda_rat_rectangle operator()() const
  { 
   leda_rat_rectangle r;
   return r;
  }  

  leda_rat_rectangle operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
      (Construct_leda_rat_rectangle::ev_leda_rat_point, p1, p2);
#endif   
   return leda_rat_rectangle(p1,p2);
  }     
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_rectangle::ev_leda_rat_point; 
#endif  

class Construct_leda_rat_object {
public:
  typedef Arity_tag< 1 > Arity;
  typedef CGAL::Object           result_type;


  template<class T>
  CGAL::Object operator()(const T& obj) const
  {
     return CGAL::make_object(obj);
  }     
};

class Construct_leda_rat_scaled_vector_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_vector           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector_integer;
  static CGAL::event ev_leda_rat_vector_quotient;   
#endif     

  leda_rat_vector operator()(const leda_rat_vector& v, const leda_integer& scale) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_vector&, const leda_integer&> \
      (Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_integer, v, scale);
#endif   
  
   return scale * v;
  }      

  leda_rat_vector operator()(const leda_rat_vector& v, const CGAL::Quotient<leda_integer>& scale) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_vector&, const CGAL::Quotient<leda_integer>& > \
      (Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_quotient, v, scale);
#endif  
   leda_rational fkt(scale.numerator(), scale.denominator());
   return fkt * v;
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_integer;
CGAL::event Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_quotient;   
#endif 

class Construct_leda_rat_translated_point_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_vector;
#endif 

  leda_rat_point operator()(const leda_rat_point& p, const leda_rat_vector& v) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_point&, const leda_rat_vector&> \
      (Construct_leda_rat_translated_point_2::ev_leda_rat_point_vector, p, v);
#endif     
    return p.translate(v);
  }    
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_translated_point_2::ev_leda_rat_point_vector;
#endif 


class Construct_leda_rat_point_on_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_int;
  static CGAL::event ev_leda_rat_ray_int;
  static CGAL::event ev_leda_rat_segment_int;    
#endif   

  leda_rat_point operator()(const leda_rat_line& l, int i =0) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_line&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_line_int, l, i);
#endif   
    leda_rat_segment s = l.seg();
    leda_rat_vector  v = s.to_vector();
      
    return s.start() + (leda_integer(i) * v);   
  }

  leda_rat_point operator()(const leda_rat_ray& r, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_ray&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_ray_int, r, i);
#endif    
    if (i==0) return r.source();
    return r.point2(); // return a point different from the source ...
  }
    
  leda_rat_point operator()(const leda_rat_segment& s, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_segment&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_segment_int, s, i);
#endif    
    i = i % 2;
    if (i==0) return s.start();
      
    return s.end();
  }        
    
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_point_on_2::ev_leda_rat_line_int;
CGAL::event Construct_leda_rat_point_on_2::ev_leda_rat_ray_int;
CGAL::event Construct_leda_rat_point_on_2::ev_leda_rat_segment_int;    
#endif   

class Construct_leda_rat_projected_point_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_point; 
#endif  

  // orthogonal projection of p onto l
  leda_rat_point operator()(const leda_rat_line& l, const leda_rat_point& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_line&, const leda_rat_point&> \
      (Construct_leda_rat_projected_point_2::ev_leda_rat_line_point, l, p);
#endif  
    leda_rat_segment s = l.perpendicular(p);
    return s.end();
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_projected_point_2::ev_leda_rat_line_point; 
#endif

class Construct_leda_rat_vertex_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_triangle_int;
  static CGAL::event ev_leda_rat_rectangle_int;
  static CGAL::event ev_leda_rat_segment_int;    
#endif    

  leda_rat_point operator()(const leda_rat_segment& s, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_segment&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_segment_int, s, i);
#endif   
    i = i % 2;
      
    if (i==0) return s.start();
    return s.end();
  }
    
  leda_rat_point operator()(const leda_rat_rectangle& r, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_rectangle&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_rectangle_int, r, i);
#endif   
    i = i % 4;
    return r[i+1];
  }
    
  leda_rat_point operator()(const leda_rat_triangle& t, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_triangle&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_triangle_int, t, i);
#endif    
    i = i % 3;
    return t[i+1];
  }    
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_vertex_2::ev_leda_rat_triangle_int;
CGAL::event Construct_leda_rat_vertex_2::ev_leda_rat_rectangle_int;
CGAL::event Construct_leda_rat_vertex_2::ev_leda_rat_segment_int;    
#endif 

class Construct_leda_rat_supporting_line_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_line  result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_segment;    
#endif    

  leda_rat_line operator()(const leda_rat_ray& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_ray&> \
      (Construct_leda_rat_supporting_line_2::ev_leda_rat_ray, r);
#endif    
     return leda_rat_line(r.point1(), r.point2());
  }
    
  leda_rat_line operator()(const leda_rat_segment& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_segment&> \
      (Construct_leda_rat_supporting_line_2::ev_leda_rat_segment, s);
#endif   
     return leda_rat_line(s.start(), s.end());
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_supporting_line_2::ev_leda_rat_ray;
CGAL::event Construct_leda_rat_supporting_line_2::ev_leda_rat_segment;    
#endif 

class Construct_leda_rat_perpendicular_direction_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef LEDA_NAMESPACE_NAME::rat_direction   result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction_orientation;   
#endif    

  LEDA_NAMESPACE_NAME::rat_direction operator()(const LEDA_NAMESPACE_NAME::rat_direction& d,
                                                CGAL::Orientation ori) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const LEDA_NAMESPACE_NAME::rat_direction&, CGAL::Orientation> \
        (Construct_leda_rat_perpendicular_direction_2::ev_leda_rat_direction_orientation, d, ori);
#endif  
  
      leda_rat_vector v = d.get_vector();
  
      if (ori == CGAL::COUNTERCLOCKWISE)
      {
        return v.rotate90(1);
      }
      else { // clockwise ...
        return v.rotate90(-1);
      }
      
      // collinear is not allowed
      return LEDA_NAMESPACE_NAME::rat_direction(v);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_perpendicular_direction_2::ev_leda_rat_direction_orientation;   
#endif

class Construct_leda_rat_perpendicular_vector_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_vector           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector_orientation;   
#endif  

  leda_rat_vector operator()(const leda_rat_vector& v, CGAL::Orientation ori) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const leda_rat_vector&, CGAL::Orientation> \
        (Construct_leda_rat_perpendicular_vector_2::ev_leda_rat_vector_orientation, v, ori);
#endif  
  
      if (ori == CGAL::COUNTERCLOCKWISE)
      {
        return v.rotate90(1);
      }
      else { // clockwise ...
        return v.rotate90(-1);
      }
      
      // collinear is not allowed
      return v;
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_perpendicular_vector_2::ev_leda_rat_vector_orientation;   
#endif

class Construct_leda_rat_perpendicular_line_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_line           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_point;   
#endif 

  leda_rat_line operator()(const leda_rat_line& l, const leda_rat_point& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const leda_rat_line&, const leda_rat_point&> \
        (Construct_leda_rat_perpendicular_line_2::ev_leda_rat_line_point, l, p);
#endif  
  
    // construct perp. line through p; rotation ccw by 90 degrees
    int ori = l.side_of(p);
      
    if (ori == 0) { // special case: collinear
      return l.rotate90(p,1);
    }
      
    leda_rat_segment s = l.perpendicular(p);
      
    if (ori == +1){
      return leda_rat_line(s.reverse());
    }
      
    return leda_rat_line(s);    
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_perpendicular_line_2::ev_leda_rat_line_point;   
#endif

class Construct_leda_rat_midpoint_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif  

  leda_rat_point operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
       (Construct_leda_rat_midpoint_2::ev_leda_rat_point, p1, p2);
#endif  
     return LEDA_NAMESPACE_NAME::midpoint(p1,p2);
  }
};    

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_midpoint_2::ev_leda_rat_point;   
#endif  

class Construct_leda_rat_center_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_point           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle;   
#endif  

#if defined(CGAL_COMPATIBLE_CIRCLES)
  leda_rat_point operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& C) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&> \
       (Construct_leda_rat_center_2::ev_leda_rat_circle, C);
#endif    
     return C.center();
  }
#else
  leda_rat_point operator()(const leda_rat_circle& C) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&> \
       (Construct_leda_rat_center_2::ev_leda_rat_circle, C);
#endif   
     return C.center();
  }
#endif    
};  

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_center_2::ev_leda_rat_circle;   
#endif 

class Construct_leda_rat_centroid_2 {
public:
  typedef leda_rat_point           result_type;
  typedef Arity_tag< 3 >          Arity;  
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_point_point;   
  static CGAL::event ev_leda_rat_point_point_point_point;     
#endif  

  leda_rat_point operator()(const leda_rat_point& p1, const leda_rat_point& p2, const leda_rat_point& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
       (Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point, p1, p2, p3);
#endif   
  
     // sum up coordinates, divide by 3
     leda_rational x = (p1.xcoord() + p2.xcoord() + p3.xcoord() )/ leda_rational(3.0);
     leda_rational y = (p1.ycoord() + p2.ycoord() + p3.ycoord() )/ leda_rational(3.0);
       
     return leda_rat_point(x,y);
  }

  leda_rat_point operator()(const leda_rat_point& p1, const leda_rat_point& p2, 
                            const leda_rat_point& p3, const leda_rat_point& p4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
       (Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point_point, p1, p2, p3, p4);
#endif     
  
     // sum up coordinates, divide by 4
     leda_rational x = (p1.xcoord() + p2.xcoord() + p3.xcoord() + p4.xcoord() )/ leda_rational(4.0);
     leda_rational y = (p1.ycoord() + p2.ycoord() + p3.ycoord() + p4.ycoord() )/ leda_rational(4.0);
       
     return leda_rat_point(x,y);       
  }
}; 

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point;   
CGAL::event Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point_point;     
#endif  

class Construct_leda_rat_circumcenter_2 {
public:
  typedef Arity_tag< 3 > Arity;
  typedef leda_rat_point           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif

  leda_rat_point operator()(const leda_rat_point& p1, const leda_rat_point& p2, const leda_rat_point& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&, const leda_rat_point&> \
       (Construct_leda_rat_circumcenter_2::ev_leda_rat_point, p1, p2, p3);
#endif  
  
    leda_rat_circle C(p1,p2,p3);
    return C.center();
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_circumcenter_2::ev_leda_rat_point;   
#endif

class Construct_leda_rat_bisector_2 {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_line           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif

  leda_rat_line operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_point&, const leda_rat_point&> \
       (Construct_leda_rat_bisector_2::ev_leda_rat_point, p1, p2);
#endif   
     return LEDA_NAMESPACE_NAME::p_bisector(p1,p2);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_bisector_2::ev_leda_rat_point;   
#endif

class Construct_leda_rat_opposite_direction_2 {
public:

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction;
#endif  

  typedef Arity_tag< 1 > Arity;
  typedef LEDA_NAMESPACE_NAME::rat_direction           result_type;

  LEDA_NAMESPACE_NAME::rat_direction operator()(const LEDA_NAMESPACE_NAME::rat_direction& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const LEDA_NAMESPACE_NAME::rat_direction&>(Construct_leda_rat_opposite_direction_2::ev_leda_rat_direction, dir);
#endif    
      leda_rat_vector v = dir.get_vector(); 
      return LEDA_NAMESPACE_NAME::rat_direction(-v);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event  Construct_leda_rat_opposite_direction_2::ev_leda_rat_direction;
#endif

// attention - we need a special functor for a special segment ...
// (for result type)

class Construct_leda_rat_opposite_segment_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_segment           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_segment;   
#endif
  
  leda_rat_segment operator()(const leda_rat_segment& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_segment&> \
       (Construct_leda_rat_opposite_segment_2::ev_leda_rat_segment, s);
#endif 
    return s.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_segment_2::ev_leda_rat_segment;   
#endif

class Construct_leda_rat_opposite_ray_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_ray           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_ray;   
#endif

  leda_rat_ray operator()(const leda_rat_ray& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_ray&> \
       (Construct_leda_rat_opposite_ray_2::ev_leda_rat_ray, r);
#endif   
    return r.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_ray_2::ev_leda_rat_ray;   
#endif

class Construct_leda_rat_opposite_line_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_line           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line;   
#endif

  leda_rat_line operator()(const leda_rat_line& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_line&> \
       (Construct_leda_rat_opposite_line_2::ev_leda_rat_line, l);
#endif    
    return l.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_line_2::ev_leda_rat_line;   
#endif

class Construct_leda_rat_opposite_triangle_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_triangle       result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_triangle;   
#endif

  leda_rat_triangle operator()(const leda_rat_triangle& t) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_triangle&> \
       (Construct_leda_rat_opposite_triangle_2::ev_leda_rat_triangle, t);
#endif   
    return t.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_triangle_2::ev_leda_rat_triangle;   
#endif

class Construct_leda_rat_opposite_circle_2 {
public:
  typedef Arity_tag< 1 > Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle;   
#endif  
  
#if defined(CGAL_COMPATIBLE_CIRCLES)
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle         result_type;

  LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&> \
       (Construct_leda_rat_opposite_circle_2::ev_leda_rat_circle, c);
#endif  
    leda_rat_point p  = c.center();
    leda_rational  sq = c.sqr_radius();
    CGAL::Orientation ori = c.orientation();
    
    ori = reverse_orientation(ori);
  
    return LEDA_NAMESPACE_NAME::cgal_rat_circle(p,sq,ori); 
  }
#else  
  typedef leda_rat_circle         result_type;

  leda_rat_circle operator()(const leda_rat_circle& c) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_circle&> \
       (Construct_leda_rat_opposite_circle_2::ev_leda_rat_circle, c);
#endif  
    return c.reverse(); 
  }
#endif  
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_circle_2::ev_leda_rat_circle;   
#endif 

class Construct_leda_rat_opposite_vector_2 {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_vector         result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector;   
#endif

  leda_rat_vector operator()(const leda_rat_vector& v) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_vector&> \
       (Construct_leda_rat_opposite_vector_2::ev_leda_rat_vector, v);
#endif  
    return -v; 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
CGAL::event Construct_leda_rat_opposite_vector_2::ev_leda_rat_vector;   
#endif

CGAL_END_NAMESPACE

#endif




