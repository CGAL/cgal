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

template<class K>
class Construct_leda_rat_point {

  typedef typename K::Point_2    Point_2;
  typedef typename K::FT         FT;
  typedef typename K::RT         RT;    

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_2           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_origin;
#endif   
  
  //----------------------------------------------------------------------------------
  //undocumented
  Point_2 operator()() const
  { 
   Point_2 p;
   return p;
  } 
  
  Point_2 operator()(const FT& x, const FT& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const FT&,const FT&>(Construct_leda_rat_point::ev_leda_rational, x, y);
#endif     
    Point_2 p(x,y);
    return p;
  }  
  
  Point_2 operator()(const RT& x, const RT& y, const RT& w) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const RT&,const RT&, const RT&> \
       (Construct_leda_rat_point::ev_leda_integer, x, y, w);
#endif   
    Point_2 p(x,y,w);
    return p;
  } 
  
  //---------------------------------------------------------------------------------- 
  Point_2 operator()(const CGAL::Origin& orig) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const CGAL::Origin&> \
       (Construct_leda_rat_point::ev_origin, orig);
#endif    
    return Point_2(0,0,1); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_point<K>::ev_leda_rational;
template<class K> CGAL::event Construct_leda_rat_point<K>::ev_leda_integer;
template<class K> CGAL::event Construct_leda_rat_point<K>::ev_origin;
#endif   

template<class K>
class Construct_leda_rat_vector {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Vector_2   Vector_2; 
  typedef typename K::FT         FT;
  typedef typename K::RT         RT;       

public:
  typedef Vector_2           result_type;
  typedef Arity_tag< 2 > Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_leda_rat_point;
  static CGAL::event ev_null_vector;
#endif   

  //----------------------------------------------------------------------------------   
  //undocumented
  Vector_2 operator()() const
  { 
   Vector_2 v(2);
   return v;
  }  
  
  Vector_2 operator()(const FT& x, const FT& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const FT&,const FT&>(Construct_leda_rat_vector::ev_leda_rational, x, y);
#endif   
    Vector_2 v(x,y);
    return v;
  }  
  
  Vector_2 operator()(const RT& x, const RT& y, const RT& w) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const RT&,const RT&,const RT&> \
      (Construct_leda_rat_vector::ev_leda_integer, x, y, w);
#endif    
    Vector_2 v(x,y,w);
    return v;
  } 
  
  //----------------------------------------------------------------------------------   

  Vector_2 operator()(const Point_2& a, const Point_2& b) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&,const Point_2&> \
      (Construct_leda_rat_vector::ev_leda_rat_point, a, b);
#endif   
    return b-a; 
  }
   
  Vector_2 operator()(const CGAL::Null_vector& nv) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const CGAL::Null_vector&> \
      (Construct_leda_rat_vector::ev_null_vector, nv);
#endif     
    return Vector_2(2); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_vector<K>::ev_leda_rational;
template<class K> CGAL::event Construct_leda_rat_vector<K>::ev_leda_integer;
template<class K> CGAL::event Construct_leda_rat_vector<K>::ev_leda_rat_point;
template<class K> CGAL::event Construct_leda_rat_vector<K>::ev_null_vector;
#endif   

template<class K>
class Construct_leda_rat_direction {

  typedef typename K::Direction_2  Direction_2;
  typedef typename K::Vector_2   Vector_2; 
  typedef typename K::FT         FT;
  typedef typename K::RT         RT; 
  typedef typename K::Line_2     Line_2;
  typedef typename K::Ray_2      Ray_2;
  typedef typename K::Segment_2  Segment_2;          

public:
  typedef Arity_tag< 1 > Arity;
  typedef Direction_2           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rational;
  static CGAL::event ev_leda_rat_vector;
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;    
  static CGAL::event ev_leda_rat_segment;
#endif   
  
  //undocumented
  Direction_2 operator()() const
  { 
   Direction_2 d(2);
   return d;
  } 
  
  Direction_2 operator()(const FT& x, const FT& y) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const FT&,const FT&> \
      (Construct_leda_rat_direction::ev_leda_rational, x, y);
#endif   
    Direction_2 d(x,y);
    return d;
  }        

  Direction_2 operator()(const Vector_2& v) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Vector_2&> \
      (Construct_leda_rat_direction::ev_leda_rat_vector, v);
#endif   
    return Direction_2(v); 
  }
    
  Direction_2 operator()(const Line_2& l) const
  {   
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Line_2&> \
      (Construct_leda_rat_direction::ev_leda_rat_line, l);
#endif   
      Vector_2 v = l.point2()-l.point1();
      return Direction_2(v); 
  }
    
  Direction_2 operator()(const Ray_2& r) const
  {  
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&> \
      (Construct_leda_rat_direction::ev_leda_rat_ray, r);
#endif   
     Vector_2 v = r.point2()-r.point1(); 
     return Direction_2(v);
  }
    
  Direction_2 operator()(const Segment_2& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&>(Construct_leda_rat_direction::ev_leda_rat_segment, s);
#endif  
    Vector_2 v = s.end()-s.start();
    return Direction_2(v); 
  }            
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Construct_leda_rat_direction<K>::ev_leda_rational;
template<class K> CGAL::event  Construct_leda_rat_direction<K>::ev_leda_rat_vector;
template<class K> CGAL::event  Construct_leda_rat_direction<K>::ev_leda_rat_line;
template<class K> CGAL::event  Construct_leda_rat_direction<K>::ev_leda_rat_ray;    
template<class K> CGAL::event  Construct_leda_rat_direction<K>::ev_leda_rat_segment;
#endif


template<class K>
class Construct_leda_rat_segment {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Segment_2 Segment_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Segment_2      result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif    
  
  //undocumented
  Segment_2 operator()() const
  { 
   Segment_2 s;
   return s;
  }   

  Segment_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&> \
      (Construct_leda_rat_segment::ev_leda_rat_point, p1, p2);
#endif    
    return Segment_2(p1,p2); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_segment<K>::ev_leda_rat_point;
#endif   

template<class K>
class Construct_leda_rat_line {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Line_2    Line_2; 
  typedef typename K::RT        RT;
  typedef typename K::Direction_2  Direction_2;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Ray_2     Ray_2;         

public:
  typedef Line_2           result_type;
  typedef Arity_tag< 2 >   Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_integer;
  static CGAL::event ev_leda_rat_point; 
  static CGAL::event ev_leda_rat_point_direction; 
  static CGAL::event ev_leda_rat_segment;
  static CGAL::event ev_leda_rat_ray;      
#endif   
  
  //undocumented
  Line_2 operator()() const
  { 
   Line_2 l;
   return l;
  }    

  Line_2 operator()(const RT& a, const RT& b, const RT& c) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const RT&, const RT&, const RT&> \
      (Construct_leda_rat_line::ev_leda_integer, a, b, c);
#endif   
  
      if (b == 0){ // par. to y - axis
       Point_2 p1(-c,1,a);
       Point_2 p2(-c,0,a);
       return Line_2(p1,p2);      
      }
      if (a == 0){ // par. to x - axis
       Point_2 p1(0,-c,b);
       Point_2 p2(1,-c,b); 
       return Line_2(p1,p2);        
      }
      // a == 0 and c == 0 not allowed
      
      Point_2 p1(0,-c,b);
      Point_2 p2(-c,0,a);
      
      return Line_2(p1,p2);
  }

  Line_2 operator()(const Point_2& p1, const Point_2& p2) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&> \
      (Construct_leda_rat_line::ev_leda_rat_point, p1, p2);
#endif  
    return Line_2(p1,p2); 
  }

  Line_2 operator()(const Point_2& p1, const Direction_2& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Direction_2&> \
      (Construct_leda_rat_line::ev_leda_rat_point_direction, p1, dir);
#endif  
    return Line_2(p1, dir.get_vector()); 
  }

  Line_2 operator()(const Segment_2& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&> \
      (Construct_leda_rat_line::ev_leda_rat_segment, s);
#endif    
    return Line_2(s.start(), s.end()); 
  }

  Line_2 operator()(const Ray_2& r) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&> \
      (Construct_leda_rat_line::ev_leda_rat_ray, r);
#endif     
    return Line_2(r); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_line<K>::ev_leda_integer;
template<class K> CGAL::event Construct_leda_rat_line<K>::ev_leda_rat_point; 
template<class K> CGAL::event Construct_leda_rat_line<K>::ev_leda_rat_point_direction; 
template<class K> CGAL::event Construct_leda_rat_line<K>::ev_leda_rat_segment;
template<class K> CGAL::event Construct_leda_rat_line<K>::ev_leda_rat_ray;      
#endif 

template<class K>
class Construct_leda_rat_ray {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Direction_2  Direction_2;
  typedef typename K::Ray_2     Ray_2; 

public:
  typedef Arity_tag< 2 >  Arity;
  typedef Ray_2           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
  static CGAL::event ev_leda_rat_point_direction;      
#endif    
  
  //undocumented
  Ray_2 operator()() const
  { 
   Ray_2 r;
   return r;
  }    

  Ray_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&> \
      (Construct_leda_rat_ray::ev_leda_rat_point, p1, p2);
#endif        
    return Ray_2(p1,p2); 
  }
    
  Ray_2 operator()(const Point_2& p1, const Direction_2& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Direction_2&> \
      (Construct_leda_rat_ray::ev_leda_rat_point_direction, p1, dir);
#endif    
    return Ray_2(p1, dir.get_vector()); 
  }    
};


#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_ray<K>::ev_leda_rat_point; 
template<class K> CGAL::event Construct_leda_rat_ray<K>::ev_leda_rat_point_direction;      
#endif 


#if defined(CGAL_COMPATIBLE_CIRCLES)
template<class K>
class Construct_leda_rat_circle {

  typedef typename K::Point_2   Point_2;
  typedef typename K::FT        FT;

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
    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& center, 
                                                    const FT&  sqrad, 
						    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const FT&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_rational_orientation, center, sqrad, ori);
#endif    
      return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,sqrad,ori);
    }    

    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, 
                                                    const Point_2& p2, 
						    const Point_2& p3) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&, const Point_2&> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_point, p1, p2, p3);
#endif    
     int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
     leda_rat_circle C(p1,p2,p3);
     Point_2  center = C.center();
     
     CGAL::Orientation cg_ori;
     switch(ori){
        case -1: { cg_ori = CGAL::RIGHTTURN; break; }
	case  0: { cg_ori = CGAL::COLLINEAR; break; }
	case  1: { cg_ori = CGAL::LEFTTURN; break; }
     }
    
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,p1,cg_ori);
    }
    
    // diameter version ...    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, const Point_2& p2, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_orientation, p1, p2, ori);
#endif        
     Point_2 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(m,p1,ori);
    }
    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_orientation, p1, ori);
#endif     
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(p1,ori);
    }        
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_rational_orientation; 
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_point_point; 
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_point_orientation;
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_orientation;       
#endif 


#else
// use LEDA circles ...
// in this case we cannot provide some constructions ...

template<class K>
class Construct_leda_rat_circle {

  typedef typename K::Point_2   Point_2;
  typedef typename K::FT        FT;

public:
  typedef leda_rat_circle         result_type;
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

    leda_rat_circle operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, const Point_2&> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_point, p1, p2, p3);
#endif     
     return leda_rat_circle(p1,p2,p3);
    }
    
    leda_rat_circle operator()(const Point_2& p1, const Point_2& p2, CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_point_orientation, p1, p2, ori);
#endif    
     // p1p2 is the diameter ...
     Point_2 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     leda_rat_circle C(m,p1);
     
     // check circle orientation ...
     int cori = (ori ==  CGAL::COUNTERCLOCKWISE) ?  -1 : 1;
     
     if (C.orientation() == cori) return C;
     else C = C.reverse();
     
     return C;
    }
    
    leda_rat_circle operator()(const Point_2& p1, CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, CGAL::Orientation> \
      (Construct_leda_rat_circle::ev_leda_rat_point_orientation, p1, ori);
#endif     
     return leda_rat_circle(p1,p1);
    }        
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_point_point; 
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_point_orientation;
template<class K> CGAL::event Construct_leda_rat_circle<K>::ev_leda_rat_point_orientation;       
#endif

#endif

template<class K>
class Construct_leda_rat_triangle {

  typedef typename K::Point_2     Point_2;
  typedef typename K::Triangle_2  Triangle_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef Triangle_2     result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
#endif    

  //undocumented
  Triangle_2 operator()() const
  { 
   Triangle_2  t;
   return t;
  }  

  Triangle_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&, const Point_2&> \
      (Construct_leda_rat_triangle::ev_leda_rat_point, p1, p2, p3);
#endif   
   return Triangle_2(p1,p2,p3);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_triangle<K>::ev_leda_rat_point; 
#endif 

template<class K>
class Construct_leda_rat_rectangle {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;

public:
  typedef Arity_tag< 2 >    Arity;
  typedef Iso_rectangle_2   result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point; 
#endif    

  //undocumented
  Iso_rectangle_2 operator()() const
  { 
   Iso_rectangle_2 r;
   return r;
  }  

  Iso_rectangle_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Point_2&> \
      (Construct_leda_rat_rectangle::ev_leda_rat_point, p1, p2);
#endif   
   return Iso_rectangle_2(p1,p2);
  }     
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_rectangle<K>::ev_leda_rat_point; 
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

template<class K>
class Construct_leda_rat_scaled_vector_2 {

  typedef typename K::Vector_2  Vector_2;
  typedef typename K::RT        RT;
  typedef typename K::FT        FT;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_2       result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector_integer;
  static CGAL::event ev_leda_rat_vector_quotient;   
#endif     

  Vector_2 operator()(const Vector_2& v, const RT& scale) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Vector_2&, const RT&> \
      (Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_integer, v, scale);
#endif   
  
   return scale * v;
  }      

  Vector_2 operator()(const Vector_2& v, const CGAL::Quotient<RT>& scale) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Vector_2&, const CGAL::Quotient<RT>& > \
      (Construct_leda_rat_scaled_vector_2::ev_leda_rat_vector_quotient, v, scale);
#endif  
   FT fkt(scale.numerator(), scale.denominator());
   return fkt * v;
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_scaled_vector_2<K>::ev_leda_rat_vector_integer;
template<class K> CGAL::event Construct_leda_rat_scaled_vector_2<K>::ev_leda_rat_vector_quotient;   
#endif 


template<class K>
class Construct_leda_rat_translated_point_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Vector_2      Vector_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_vector;
#endif 

  Point_2 operator()(const Point_2& p, const Vector_2& v) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&, const Vector_2&> \
      (Construct_leda_rat_translated_point_2::ev_leda_rat_point_vector, p, v);
#endif     
    return p.translate(v);
  }    
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_translated_point_2<K>::ev_leda_rat_point_vector;
#endif 

template<class K>
class Construct_leda_rat_point_on_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Vector_2      Vector_2;  
  typedef typename K::Line_2        Line_2;
  typedef typename K::Ray_2         Ray_2;  
  typedef typename K::Segment_2     Segment_2; 
  typedef typename K::RT            RT;     

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_int;
  static CGAL::event ev_leda_rat_ray_int;
  static CGAL::event ev_leda_rat_segment_int;    
#endif   

  Point_2 operator()(const Line_2& l, int i = 0) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Line_2&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_line_int, l, i);
#endif   
    Segment_2 s = l.seg();
    Vector_2  v = s.to_vector();
      
    return s.start() + (RT(i) * v);   
  }

  Point_2 operator()(const Ray_2& r, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_ray_int, r, i);
#endif    
    if (i==0) return r.source();
    return r.point2(); // return a point different from the source ...
  }
    
  Point_2 operator()(const Segment_2& s, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&, int> \
      (Construct_leda_rat_point_on_2::ev_leda_rat_segment_int, s, i);
#endif    
    i = i % 2;
    if (i==0) return s.start();
      
    return s.end();
  }        
    
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_point_on_2<K>::ev_leda_rat_line_int;
template<class K> CGAL::event Construct_leda_rat_point_on_2<K>::ev_leda_rat_ray_int;
template<class K> CGAL::event Construct_leda_rat_point_on_2<K>::ev_leda_rat_segment_int;    
#endif   

template<class K> 
class Construct_leda_rat_projected_point_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Line_2        Line_2;
  typedef typename K::Segment_2     Segment_2; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_point; 
#endif  

  // orthogonal projection of p onto l
  Point_2 operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Line_2&, const Point_2&> \
      (Construct_leda_rat_projected_point_2::ev_leda_rat_line_point, l, p);
#endif  
    Segment_2 s = l.perpendicular(p);
    return s.end();
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_projected_point_2<K>::ev_leda_rat_line_point; 
#endif

template<class K>
class Construct_leda_rat_vertex_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Triangle_2    Triangle_2;
  typedef typename K::Segment_2     Segment_2; 
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_triangle_int;
  static CGAL::event ev_leda_rat_rectangle_int;
  static CGAL::event ev_leda_rat_segment_int;    
#endif    

  Point_2 operator()(const Segment_2& s, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_segment_int, s, i);
#endif   
    i = i % 2;
      
    if (i==0) return s.start();
    return s.end();
  }
    
  Point_2 operator()(const Iso_rectangle_2& r, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Iso_rectangle_2&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_rectangle_int, r, i);
#endif   
    i = i % 4;
    return r[i+1];
  }
    
  Point_2 operator()(const Triangle_2& t, int i) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Triangle_2&, int> \
      (Construct_leda_rat_vertex_2::ev_leda_rat_triangle_int, t, i);
#endif    
    i = i % 3;
    return t[i+1];
  }    
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_vertex_2<K>::ev_leda_rat_triangle_int;
template<class K> CGAL::event Construct_leda_rat_vertex_2<K>::ev_leda_rat_rectangle_int;
template<class K> CGAL::event Construct_leda_rat_vertex_2<K>::ev_leda_rat_segment_int;    
#endif 

template<class K>
class Construct_leda_rat_supporting_line_2 {

  typedef typename K::Line_2        Line_2;
  typedef typename K::Ray_2         Ray_2;
  typedef typename K::Segment_2     Segment_2; 

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_2  result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_segment;    
#endif    

  Line_2 operator()(const Ray_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&> \
      (Construct_leda_rat_supporting_line_2::ev_leda_rat_ray, r);
#endif    
     return Line_2(r.point1(), r.point2());
  }
    
  Line_2 operator()(const Segment_2& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&> \
      (Construct_leda_rat_supporting_line_2::ev_leda_rat_segment, s);
#endif   
     return Line_2(s.start(), s.end());
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_supporting_line_2<K>::ev_leda_rat_ray;
template<class K> CGAL::event Construct_leda_rat_supporting_line_2<K>::ev_leda_rat_segment;    
#endif 

template<class K>
class Construct_leda_rat_perpendicular_direction_2 {

  typedef typename K::Direction_2        Direction_2;
  typedef typename K::Vector_2           Vector_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Direction_2    result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction_orientation;   
#endif    

  Direction_2 operator()(const Direction_2& d, CGAL::Orientation ori) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Direction_2&, CGAL::Orientation> \
        (Construct_leda_rat_perpendicular_direction_2::ev_leda_rat_direction_orientation, d, ori);
#endif  
  
      Vector_2 v = d.get_vector();
  
      if (ori == CGAL::COUNTERCLOCKWISE)
      {
        return v.rotate90(1);
      }
      else { // clockwise ...
        return v.rotate90(-1);
      }
      
      // collinear is not allowed
      return Direction_2(v);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_perpendicular_direction_2<K>::ev_leda_rat_direction_orientation;   
#endif

template<class K>
class Construct_leda_rat_perpendicular_vector_2 {

  typedef typename K::Vector_2           Vector_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_2       result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector_orientation;   
#endif  

  Vector_2 operator()(const Vector_2& v, CGAL::Orientation ori) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Vector_2&, CGAL::Orientation> \
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
template<class K> CGAL::event Construct_leda_rat_perpendicular_vector_2<K>::ev_leda_rat_vector_orientation;   
#endif

template<class K>
class Construct_leda_rat_perpendicular_line_2 {

  typedef typename K::Line_2           Line_2;
  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Line_2         result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_point;   
#endif 

  Line_2 operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Line_2&, const Point_2&> \
        (Construct_leda_rat_perpendicular_line_2::ev_leda_rat_line_point, l, p);
#endif  
  
    // construct perp. line through p; rotation ccw by 90 degrees
    int ori = l.side_of(p);
      
    if (ori == 0) { // special case: collinear
      return l.rotate90(p,1);
    }
      
    Segment_2 s = l.perpendicular(p);
      
    if (ori == +1){
      return Line_2(s.reverse());
    }
      
    return Line_2(s);    
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_perpendicular_line_2<K>::ev_leda_rat_line_point;   
#endif

template<class K>
class Construct_leda_rat_midpoint_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif  

  Point_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&> \
       (Construct_leda_rat_midpoint_2::ev_leda_rat_point, p1, p2);
#endif  
     return LEDA_NAMESPACE_NAME::midpoint(p1,p2);
  }
};    

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_midpoint_2<K>::ev_leda_rat_point;   
#endif  

template<class K>
class Construct_leda_rat_center_2 {

  typedef typename K::Point_2          Point_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_2        result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle;   
#endif  

#if defined(CGAL_COMPATIBLE_CIRCLES)
  Point_2 operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& C) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&> \
       (Construct_leda_rat_center_2::ev_leda_rat_circle, C);
#endif    
     return C.center();
  }
#else
  Point_2 operator()(const leda_rat_circle& C) const
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
template<class K> CGAL::event Construct_leda_rat_center_2<K>::ev_leda_rat_circle;   
#endif 

template<class K>
class Construct_leda_rat_centroid_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::FT               FT;

public:
  typedef Point_2           result_type;
  typedef Arity_tag< 3 >    Arity;  
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_point_point;   
  static CGAL::event ev_leda_rat_point_point_point_point;     
#endif  

  Point_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, const Point_2&> \
       (Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point, p1, p2, p3);
#endif   
  
     // sum up coordinates, divide by 3
     FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() )/ FT(3.0);
     FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() )/ FT(3.0);
       
     return Point_2(x,y);
  }

  Point_2 operator()(const Point_2& p1, const Point_2& p2, 
                     const Point_2& p3, const Point_2& p4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, const Point_2&, const Point_2&> \
       (Construct_leda_rat_centroid_2::ev_leda_rat_point_point_point_point, p1, p2, p3, p4);
#endif     
  
     // sum up coordinates, divide by 4
     FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() + p4.xcoord() )/ FT(4.0);
     FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() + p4.ycoord() )/ FT(4.0);
       
     return Point_2(x,y);       
  }
}; 

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_centroid_2<K>::ev_leda_rat_point_point_point;   
template<class K> CGAL::event Construct_leda_rat_centroid_2<K>::ev_leda_rat_point_point_point_point;     
#endif  

template<class K>
class Construct_leda_rat_circumcenter_2 {

  typedef typename K::Point_2          Point_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef Point_2        result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif

  Point_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&, const Point_2&> \
       (Construct_leda_rat_circumcenter_2::ev_leda_rat_point, p1, p2, p3);
#endif  
  
    leda_rat_circle C(p1,p2,p3);
    return C.center();
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_circumcenter_2<K>::ev_leda_rat_point;   
#endif

template<class K>
class Construct_leda_rat_bisector_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Line_2           Line_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Line_2         result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;   
#endif

  Line_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&, const Point_2&> \
       (Construct_leda_rat_bisector_2::ev_leda_rat_point, p1, p2);
#endif   
     return LEDA_NAMESPACE_NAME::p_bisector(p1,p2);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_bisector_2<K>::ev_leda_rat_point;   
#endif

template<class K>
class Construct_leda_rat_opposite_direction_2 {

  typedef typename K::Vector_2         Vector_2;
  typedef typename K::Direction_2      Direction_2;

public:

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction;
#endif  

  typedef Arity_tag< 1 > Arity;
  typedef Direction_2    result_type;

  Direction_2 operator()(const Direction_2& dir) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Direction_2&>(Construct_leda_rat_opposite_direction_2::ev_leda_rat_direction, dir);
#endif    
      Vector_2 v = dir.get_vector(); 
      return Direction_2(-v);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Construct_leda_rat_opposite_direction_2<K>::ev_leda_rat_direction;
#endif


template<class K>
class Construct_leda_rat_opposite_segment_2 {

  typedef typename K::Segment_2  Segment_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Segment_2      result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_segment;   
#endif
  
  Segment_2 operator()(const Segment_2& s) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&> \
       (Construct_leda_rat_opposite_segment_2::ev_leda_rat_segment, s);
#endif 
    return s.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_opposite_segment_2<K>::ev_leda_rat_segment;   
#endif

template<class K>
class Construct_leda_rat_opposite_ray_2 {

  typedef typename K::Ray_2  Ray_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Ray_2          result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_ray;   
#endif

  Ray_2 operator()(const Ray_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&> \
       (Construct_leda_rat_opposite_ray_2::ev_leda_rat_ray, r);
#endif   
    return r.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_opposite_ray_2<K>::ev_leda_rat_ray;   
#endif

template<class K>
class Construct_leda_rat_opposite_line_2 {

  typedef typename K::Line_2  Line_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_2         result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line;   
#endif

  Line_2 operator()(const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Line_2&> \
       (Construct_leda_rat_opposite_line_2::ev_leda_rat_line, l);
#endif    
    return l.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_opposite_line_2<K>::ev_leda_rat_line;   
#endif

template<class K>
class Construct_leda_rat_opposite_triangle_2 {

  typedef typename K::Triangle_2   Triangle_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Triangle_2     result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_triangle;   
#endif

  Triangle_2 operator()(const Triangle_2& t) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Triangle_2&> \
       (Construct_leda_rat_opposite_triangle_2::ev_leda_rat_triangle, t);
#endif   
    return t.reverse(); 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_opposite_triangle_2<K>::ev_leda_rat_triangle;   
#endif

template<class K>
class Construct_leda_rat_opposite_circle_2 {

  typedef typename K::Point_2  Point_2;
  typedef typename K::FT       FT;

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
    Point_2 p  = c.center();
    FT  sq = c.sqr_radius();
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
template<class K> CGAL::event Construct_leda_rat_opposite_circle_2<K>::ev_leda_rat_circle;   
#endif 

template<class K>
class Construct_leda_rat_opposite_vector_2 {

  typedef typename K::Vector_2 Vector_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Vector_2       result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_vector;   
#endif

  Vector_2 operator()(const Vector_2& v) const
  { 
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Vector_2&> \
       (Construct_leda_rat_opposite_vector_2::ev_leda_rat_vector, v);
#endif  
    return -v; 
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Construct_leda_rat_opposite_vector_2<K>::ev_leda_rat_vector;   
#endif

CGAL_END_NAMESPACE

#endif




