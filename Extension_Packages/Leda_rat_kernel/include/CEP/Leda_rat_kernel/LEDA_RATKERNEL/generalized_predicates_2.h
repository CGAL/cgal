#ifndef CEP_LEDA_RAT_GENERALIZED_PREDICATES_2_H
#define CEP_LEDA_RAT_GENERALIZED_PREDICATES_2_H

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>
#include <CGAL/functional_base.h>

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

// ------------------------------------------------------------------------------------------------------------
//  all predicates are templated with the kernel
// ------------------------------------------------------------------------------------------------------------

template<class K>
class Predicate_leda_rat_angle_2 {

  typedef typename K::Point_2  Point_2;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::FT       FT;

public:
  typedef Arity_tag< 3 > Arity;
  typedef CGAL::Angle    result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif      

// return sign of scalar (dot) product of the two defined vectors
  CGAL::Angle operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_angle_2::ev_leda_rat_point, p1, p2, p3);
#endif    
      Vector_2 v1 = p1-p2;
      Vector_2 v2 = p3-p2;
  
      FT s_prod = v1*v2;
   
      if (s_prod == 0) return CGAL::RIGHT;
      if (s_prod >  0) return CGAL::ACUTE;
      return CGAL::OBTUSE;  
  }

};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_angle_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K>
class Predicate_leda_rat_equal_2 {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Vector_2   Vector_2;
  typedef typename K::Direction_2 Direction_2;
  typedef typename K::Line_2      Line_2;
  typedef typename K::Ray_2       Ray_2;
  typedef typename K::Segment_2   Segment_2;
  typedef typename K::Triangle_2  Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT          FT;          

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
  static CGAL::event ev_leda_rat_segment;
  static CGAL::event ev_leda_rat_vector;
  static CGAL::event ev_leda_rat_direction;
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_circle;
  static CGAL::event ev_leda_rat_triangle;
  static CGAL::event ev_leda_rat_rectangle;
#endif    

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_point, p1, p2);
#endif    
    return ( p1 == p2 );
  }
    
  bool operator()(const Vector_2& v1, const Vector_2& v2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Vector_2&,const Vector_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_vector, v1, v2);
#endif   
    return ( v1 == v2 );
  } 
  
  bool operator()(const Direction_2& d1, const Direction_2& d2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Direction_2&,const Direction_2&> \
        (Predicate_leda_rat_equal_2::ev_leda_rat_direction, d1, d2);
#endif   
    Vector_2 v1 = d1.get_vector();
    Vector_2 v2 = d2.get_vector();
    
    CGAL_precondition( (v1.dim()==2) && (v2.dim()==2));
    
    FT xc1 = v1[0];
    FT xc2 = v2[0];
    FT yc1 = v1[1];
    FT yc2 = v2[1];
    
    // equal directions ???  
    FT q1,q2;
      
    if (xc2==0) {
      if (xc1!=0) return false;
      else q1 = 0;
    }
    else q1 = xc1/xc2;
    
    if (yc2==0) {
      if (yc1!=0) return false;
      else q2 = 0;
    }
    else q2 = yc1/yc2;
    
    if (q1==q2) return true;
    return false;
  }
  
  bool operator()(const Line_2& l1, const Line_2& l2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Line_2&,const Line_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_line, l1, l2);
#endif   
    return (l1  ==  l2 );
  }  
  
  bool operator()(const Ray_2& r1, const Ray_2& r2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Ray_2&,const Ray_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_ray, r1, r2);
#endif   
    return (r1  ==  r2);
  }   
  
  bool operator()(const Segment_2& s1, const Segment_2& s2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Segment_2&,const Segment_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_segment, s1, s2);
#endif      
    return (s1 == s2);
  } 

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c1, 
                  const LEDA_NAMESPACE_NAME::cgal_rat_circle& c2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const LEDA_NAMESPACE_NAME::cgal_rat_circle&> \
    (Predicate_leda_rat_equal_2::ev_leda_rat_circle, c1, c2);
#endif  
    return (c1 == c2);
  } 
#else 
  bool operator()(const leda_rat_circle& c1, const leda_rat_circle& c2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const leda_rat_circle&,const leda_rat_circle&>(Predicate_leda_rat_equal_2::ev_leda_rat_circle, c1, c2);
#endif   
    return (c1 == c2);
  }  
#endif       
  
  bool operator()(const Triangle_2& t1, const Triangle_2& t2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Triangle_2&,const Triangle_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_triangle, t1, t2);
#endif   
#if (__LEDA__ <= 420)  
    return (  ((Triangle_2&) t1) == ((Triangle_2&) t2));
#else
    return (t1 == t2);
#endif    
  } 
  
  bool operator()(const Iso_rectangle_2& r1, const Iso_rectangle_2& r2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Iso_rectangle_2&,const Iso_rectangle_2&>(Predicate_leda_rat_equal_2::ev_leda_rat_rectangle, r1, r2);
#endif   
    return (r1 == r2);
  }         
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_point;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_segment;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_vector;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_direction;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_line;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_ray;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_circle;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_triangle;
template<class K> CGAL::event  Predicate_leda_rat_equal_2<K>::ev_leda_rat_rectangle;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_xy_2 {
  typedef typename K::Point_2   Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;  
#endif

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_less_xy_2::ev_leda_rat_point, p1, p2);
#endif   
    return ( __My_Point_2::cmp_xy(p1,p2)  <  0 );
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_xy_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_yx_2 {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;  
#endif

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_less_yx_2::ev_leda_rat_point, p1, p2);
#endif   
    return ( __My_Point_2::cmp_yx(p1,p2)  <  0 );
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_yx_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_equal_x_2 {
  typedef typename K::Point_2    Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;  
#endif  

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_equal_x_2::ev_leda_rat_point, p1, p2);
#endif   
    return (__My_Point_2::cmp_x(p1,p2) == 0);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_equal_x_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_equal_y_2 {
  typedef typename K::Point_2    Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;  
#endif 

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_equal_y_2::ev_leda_rat_point, p1, p2);
#endif     
    return (__My_Point_2::cmp_y(p1,p2) == 0);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_equal_y_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_equal_xy_2 {
  typedef typename K::Point_2   Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;  
#endif   

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
    CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_equal_xy_2::ev_leda_rat_point, p1, p2);
#endif   
    return (__My_Point_2::cmp_xy(p1,p2) == 0);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_equal_xy_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_x_2 {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif   

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_less_x_2::ev_leda_rat_point, p1, p2);
#endif    
    if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return false;
    return (__My_Point_2::cmp_x(p1,p2) < 0);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_x_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_y_2 {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif    

  bool operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_less_y_2::ev_leda_rat_point, p1, p2);
#endif  
    if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return false;
    return (__My_Point_2::cmp_y(p1,p2) < 0);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_y_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_x_2 {

  typedef typename K::Point_2 Point_2;
  typedef typename K::Line_2  Line_2;
  
  typedef typename K::Point_2   __My_Point_2;  

public:

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
  static CGAL::event ev_leda_rat_point_line_line;
  static CGAL::event ev_leda_rat_line_line_line;
  static CGAL::event ev_leda_rat_line_line_line_line;
#endif 

  typedef Comparison_result       result_type;
  typedef Arity_tag< 2 >          Arity;  

  Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_compare_x_2::ev_leda_rat_point, p1, p2);
#endif   
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return EQUAL;  
     return ( (Comparison_result) __My_Point_2::cmp_x(p1,p2));     
  }
  
  Comparison_result operator()(const Point_2& p, const Line_2& l1, const Line_2& l2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_x_2::ev_leda_rat_point_line_line, p, l1, l2);
#endif   
    Point_2 inter;
    
    l1.intersection(l2,inter);
    
    return ( (Comparison_result) __My_Point_2::cmp_x(p,inter));
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, const Line_2& l3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_x_2::ev_leda_rat_line_line_line, l1, l2, l3);
#endif    
    // compute intersections
    Point_2 p1,p2;
    
    l1.intersection(l2, p1);
    l1.intersection(l3, p2);
        
    return ( (Comparison_result) __My_Point_2::cmp_x(p1,p2));
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, 
                               const Line_2& l3, const Line_2& l4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_x_2::ev_leda_rat_line_line_line_line, l1, l2, l3, l4);
#endif   
    // compute intersections
    Point_2 p1,p2;
    
    l1.intersection(l2, p1);
    l3.intersection(l4, p2);
        
    return ( (Comparison_result) __My_Point_2::cmp_x(p1,p2));  
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_x_2<K>::ev_leda_rat_point;
template<class K> CGAL::event  Predicate_leda_rat_compare_x_2<K>::ev_leda_rat_point_line_line;
template<class K> CGAL::event  Predicate_leda_rat_compare_x_2<K>::ev_leda_rat_line_line_line;
template<class K> CGAL::event  Predicate_leda_rat_compare_x_2<K>::ev_leda_rat_line_line_line_line;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_x_at_y_2 {

  typedef typename K::Point_2         Point_2;
  typedef typename K::Line_2          Line_2;
  typedef typename K::Segment_2       Segment_2;    
  
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Comparison_result           result_type;
  typedef Arity_tag< 2 >              Arity;  

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_line;  
  static CGAL::event ev_leda_rat_point_line_line;
  static CGAL::event ev_leda_rat_line_line_line;
  static CGAL::event ev_leda_rat_line_line_line_line;  
#endif  

  // prec.: l is not horizontal (?)

  Comparison_result operator()(const Point_2& p, const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Line_2&>(Predicate_leda_rat_compare_x_at_y_2::ev_leda_rat_point_line, p, l);
#endif   
#if (__LEDA__ <= 420)
    int ori = - ::orientation(l,p);
#else
    int ori = - l.side_of(p);
#endif    
    
    // has point1 larger y - coord than point2 ????
    if (__My_Point_2::cmp_y(l.point1(),l.point2()) == 1) ori = -ori;
    return ((Comparison_result) ori);
  }
  
  // compare x - coords of horizontal projections of p onto h1 and h2 
  Comparison_result operator()(const Point_2& p, const Line_2& h1, const Line_2& h2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&,const Line_2&, const Line_2&> \
      (Predicate_leda_rat_compare_x_at_y_2::ev_leda_rat_point_line_line, p, h1, h2);
#endif   
       Point_2  pnew(p.Y(), p.X(), p.W());
       
       Point_2 p1 = h1.point1(), p2 = h1.point2();
       Point_2 p3 = h2.point1(), p4 = h2.point2();       
       
       Segment_2 s1(Point_2(p1.Y(), p1.X(), p1.W()), Point_2(p2.Y(), p2.X(), p2.W()) );
       Segment_2 s2(Point_2(p3.Y(), p3.X(), p3.W()), Point_2(p4.Y(), p4.X(), p4.W()) ); 
       
       int res = LEDA_NAMESPACE_NAME::cmp_segments_at_xcoord(s1,s2,pnew);
       
       return ((Comparison_result) res);  
  } 
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, const Line_2& h) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Line_2&,const Line_2&, const Line_2&> \
      (Predicate_leda_rat_compare_x_at_y_2::ev_leda_rat_line_line_line, l1, l2, h);
#endif  
       Point_2 inter;
       
       // compute intersection of l1/l2
       l1.intersection(l2,inter);
       
       // compare result ...
#if (__LEDA__ <= 420)
       int ori = - ::orientation(h, inter);
#else       
       int ori = - h.side_of(inter);
#endif       
    
       // has point1 larger y - coord than point2 ????
       if (__My_Point_2::cmp_y(h.point1(),h.point2()) == 1) ori = -ori;
       return ((Comparison_result) ori);       
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2,
                               const Line_2& h1, const Line_2& h2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Line_2&,const Line_2&, const Line_2&, const Line_2&> \
      (Predicate_leda_rat_compare_x_at_y_2::ev_leda_rat_line_line_line_line, l1, l2, h1, h2);
#endif    
       Point_2 p;
       
       // compute intersection of l1/l2
       l1.intersection(l2,p);  
  
       Point_2  pnew(p.Y(), p.X(), p.W());
       
       Point_2 p1 = h1.point1(), p2 = h1.point2();
       Point_2 p3 = h2.point1(), p4 = h2.point2();       
       
       Segment_2 s1(Point_2(p1.Y(), p1.X(), p1.W()), Point_2(p2.Y(), p2.X(), p2.W()) );
       Segment_2 s2(Point_2(p3.Y(), p3.X(), p3.W()), Point_2(p4.Y(), p4.X(), p4.W()) ); 
       
       int res = LEDA_NAMESPACE_NAME::cmp_segments_at_xcoord(s1,s2,pnew);
       
       return ((Comparison_result) res);  
  }      
};


#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Predicate_leda_rat_compare_x_at_y_2<K>::ev_leda_rat_point_line;  
template<class K> CGAL::event Predicate_leda_rat_compare_x_at_y_2<K>::ev_leda_rat_point_line_line;
template<class K> CGAL::event Predicate_leda_rat_compare_x_at_y_2<K>::ev_leda_rat_line_line_line;
template<class K> CGAL::event Predicate_leda_rat_compare_x_at_y_2<K>::ev_leda_rat_line_line_line_line;  
#endif  

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_y_2 {

  typedef typename K::Point_2         Point_2;
  typedef typename K::Line_2          Line_2;  
  
  typedef typename K::Point_2         __My_Point_2;  

public:
  typedef Comparison_result           result_type;
  typedef Arity_tag< 2 > Arity;  
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
  static CGAL::event ev_leda_rat_point_line_line;
  static CGAL::event ev_leda_rat_line_line_line;
  static CGAL::event ev_leda_rat_line_line_line_line; 
#endif   

  Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_compare_y_2::ev_leda_rat_point, p1, p2);
#endif   
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return EQUAL;
     return ( (Comparison_result) __My_Point_2::cmp_y(p1,p2));     
  }
  
  Comparison_result operator()(const Point_2& p, const Line_2& l1, const Line_2& l2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_y_2::ev_leda_rat_point_line_line, p, l1, l2);
#endif     
    Point_2 inter;
    
    l1.intersection(l2,inter);
    
    return ( (Comparison_result) __My_Point_2::cmp_y(p,inter));
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, const Line_2& l3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_y_2::ev_leda_rat_line_line_line, l1, l2, l3);
#endif   
  
    // compute intersections
    Point_2 p1,p2;
    
    l1.intersection(l2, p1);
    l1.intersection(l3, p2);
        
    return ( (Comparison_result) __My_Point_2::cmp_y(p1,p2));
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, 
                               const Line_2& l3, const Line_2& l4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&,const Line_2&> \
      (Predicate_leda_rat_compare_y_2::ev_leda_rat_line_line_line_line, l1, l2, l3, l4);
#endif   
  
    // compute intersections
    Point_2 p1,p2;
    
    l1.intersection(l2, p1);
    l3.intersection(l4, p2);
        
    return ( (Comparison_result) __My_Point_2::cmp_y(p1,p2));  
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_y_2<K>::ev_leda_rat_point;
template<class K> CGAL::event  Predicate_leda_rat_compare_y_2<K>::ev_leda_rat_point_line_line;
template<class K> CGAL::event  Predicate_leda_rat_compare_y_2<K>::ev_leda_rat_line_line_line;
template<class K> CGAL::event  Predicate_leda_rat_compare_y_2<K>::ev_leda_rat_line_line_line_line;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_xy_2 {

  typedef typename K::Point_2  Point_2;
  
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif

  Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&>(Predicate_leda_rat_compare_xy_2::ev_leda_rat_point, p1, p2);
#endif   
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return EQUAL;
     return ( (Comparison_result) __My_Point_2::cmp_xy(p1,p2));
  }
    
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_xy_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_y_at_x_2 {

  typedef typename K::Point_2         Point_2;
  typedef typename K::Segment_2       Segment_2;
  typedef typename K::Line_2          Line_2;    
  
  typedef typename K::Point_2         __My_Point_2;  

public:
  typedef Comparison_result           result_type;
  typedef Arity_tag< 2 > Arity;
 
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_segment_segment;
  static CGAL::event ev_leda_rat_point_segment; 
  static CGAL::event ev_leda_rat_point_line;
  static CGAL::event ev_leda_rat_point_line_line;
  static CGAL::event ev_leda_rat_line_line_line;
  static CGAL::event ev_leda_rat_line_line_line_line;
#endif  
  
  // ------------------------------------------------------------------------------------------
  // this was copied from planar map LEDA traits ...

  bool curve_is_in_x_range(const Segment_2 & cv, const Point_2 & p) const
  {
      return
        !(((__My_Point_2::cmp_x(p, cv.source()) < 0) && (__My_Point_2::cmp_x(p, cv.target()) < 0)) ||
         ((__My_Point_2::cmp_x(p, cv.source()) > 0) && (__My_Point_2::cmp_x(p, cv.target()) > 0)));
  }
	
  bool curve_is_in_y_range(const Segment_2 & cv, const Point_2 & p) const
  { 
      return
        !(((__My_Point_2::cmp_y(p, cv.source()) < 0) && (__My_Point_2::cmp_y(p, cv.target()) < 0)) ||
         ((__My_Point_2::cmp_y(p, cv.source()) > 0) && (__My_Point_2::cmp_y(p, cv.target()) > 0)));
  }
  
    Comparison_result operator()(const Point_2 & q,
                                 const Segment_2 & cv1, const Segment_2 & cv2) const
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Segment_2&,const Segment_2&> \
         (Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_point_segment_segment,q, cv1, cv2);
#endif    
    
      // precondition ???
      if ((!curve_is_in_x_range(cv1, q)) || (!curve_is_in_x_range(cv2, q))) return EQUAL;
		
      Segment_2 cv1_ = cv1;
      Segment_2 cv2_ = cv2;
      if (__My_Point_2::cmp_xy(cv1.source(), cv1.target()) > 0)
        cv1_ = cv1.reversal();
      if (__My_Point_2::cmp_xy(cv2.source(), cv2.target()) > 0)
        cv2_ = cv2.reversal();
  		
      //  vertical curves ...
      
      if (cv1_.is_vertical()) {
        if (cv2_.is_vertical()) {
          // both cv1 and cv2 are vertical
          int res = __My_Point_2::cmp_y(cv1_.target(), cv2_.source());
          return ((res < 0) ? CGAL::SMALLER : ((res > 0) ? CGAL::LARGER : CGAL::EQUAL));
        }

        // only cv1 is vertical.
        if (LEDA_NAMESPACE_NAME::orientation(cv2_.source(), cv2_.target(), cv1_.source()) > 0)
          return CGAL::LARGER;
                      
        if (LEDA_NAMESPACE_NAME::orientation(cv2_.source(), cv2_.target(), cv1_.target()) < 0)
          return CGAL::SMALLER;
  
        return CGAL::EQUAL;
      }
                  
      if (cv2_.is_vertical()) {
        if (LEDA_NAMESPACE_NAME::orientation(cv1_.source(), cv1_.target(), cv2_.source()) > 0 )
          return CGAL::SMALLER;
                      
        if (LEDA_NAMESPACE_NAME::orientation(cv1_.source(), cv1_.target(), cv2_.target()) < 0)
          return CGAL::LARGER;
  
        return CGAL::EQUAL;  
      }
                    
      int res = LEDA_NAMESPACE_NAME::cmp_segments_at_xcoord(cv1_, cv2_, q);
      return ((res < 0) ? CGAL::SMALLER : ((res > 0) ? CGAL::LARGER : CGAL::EQUAL));
    }
  // ------------------------------------------------------------------------------------------    
    
    // precondition: p is in the x range of s ...

    Comparison_result operator()(const Point_2 & p, const Segment_2 &  s)
    {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Segment_2&>(Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_point_segment, p, s);
#endif          
      // what happens when prec. is not fulfilled ?
    
      // this handles trivial segments as well ...
      if (s.is_vertical()) {
        if ((__My_Point_2::cmp_y(p, s.source()) < 0) && (__My_Point_2::cmp_y(p, s.target()) < 0)) return CGAL::SMALLER;	  
        if ((__My_Point_2::cmp_y(p, s.source()) > 0) && (__My_Point_2::cmp_y(p, s.target()) > 0)) return CGAL::LARGER;
        return CGAL::EQUAL;
      }
      
      int ori = 0;
      
      if (__My_Point_2::cmp_x(s.source(), s.target()) < 0){
         ori = LEDA_NAMESPACE_NAME::orientation(s.source(), s.target(), p);
      }
      else ori = LEDA_NAMESPACE_NAME::orientation(s.target(), s.source(), p);
			
      if (ori < 0) return CGAL::SMALLER;
      if (ori > 0) return CGAL::LARGER;
      return  CGAL::EQUAL;
    }  
  
  // ------------------------------------------------------------------------------------------  

  // prec.: l is not vertical

  Comparison_result operator()(const Point_2& p, const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Line_2&> \
         (Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_point_line, p, l);
#endif  
    // add - here ??
#if (__LEDA__ <= 420)
    int ori = ::orientation(l,p);
#else  
    int ori = l.side_of(p);
#endif
    
    // has point1 larger x - coord than point2 ????
    if (__My_Point_2::cmp_x(l.point1(),l.point2()) == 1) ori = -ori;
    return ((Comparison_result) ori);
  }
  
  // compare y-coord of vertical projection of p onto l1 and l2 ...
  Comparison_result operator()(const Point_2& p, const Line_2& l1, const Line_2& l2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Line_2&,const Line_2&> \
         (Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_point_line_line, p, l1, l2);
#endif    
    Segment_2 s1 = l1.seg();
    Segment_2 s2 = l2.seg();
    
    return ( (Comparison_result) LEDA_NAMESPACE_NAME::cmp_segments_at_xcoord(s1,s2,p));
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, const Line_2& l3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&> \
         (Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_line_line_line, l1, l2, l3);
#endif   
  
    // compute intersection of l1/l2 ...    
    Point_2 p;
    
    l1.intersection(l2, p);
#if (__LEDA__ <= 420)
    int ori = ::orientation(l3,p);
#else    
    int ori = l3.side_of(p);
#endif    
    
    // has point1 larger x - coord than point2 ????
    if (__My_Point_2::cmp_x(l3.point1(),l3.point2()) == 1) ori = -ori;
    return ((Comparison_result) ori);    
  }
  
  Comparison_result operator()(const Line_2& l1, const Line_2& l2, 
                               const Line_2& l3, const Line_2& l4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&,const Line_2&,const Line_2&> \
         (Predicate_leda_rat_compare_y_at_x_2::ev_leda_rat_line_line_line_line, l1, l2, l3, l4);
#endif   
    // compute intersection of l1/l2 ...
    Point_2 p;
    
    l1.intersection(l2, p);
    
    Segment_2 s3 = l3.seg();
    Segment_2 s4 = l4.seg();
    
    return ( (Comparison_result) LEDA_NAMESPACE_NAME::cmp_segments_at_xcoord(s3,s4,p));      
  }      
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_point_segment_segment;
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_point_segment; 
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_point_line;
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_point_line_line;
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_line_line_line;
template<class K> CGAL::event Predicate_leda_rat_compare_y_at_x_2<K>::ev_leda_rat_line_line_line_line;
#endif  

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_distance_2 {

  typedef typename K::Point_2   Point_2;

public:
   typedef Arity_tag< 3 > Arity;
   typedef Comparison_result           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif 

   Comparison_result operator()(const Point_2& p, const Point_2& q,
                                const Point_2& r) const
   {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
        (Predicate_leda_rat_compare_distance_2::ev_leda_rat_point, p, q, r);
#endif     
     return  ((Comparison_result) p.cmp_dist(q,r));
   }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_distance_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_angle_with_x_axis_2 {

  typedef typename K::Direction_2   Direction_2;
  typedef typename K::Vector_2      Vector_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction;
#endif    

   Comparison_result operator()(const Direction_2& d1, 
                                const Direction_2& d2) const
   {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Direction_2&,const Direction_2&> \
       (Predicate_leda_rat_compare_angle_with_x_axis_2::ev_leda_rat_direction, d1, d2);
#endif    
     Vector_2 v1 = d1.get_vector();
     Vector_2 v2 = d2.get_vector();
     int res2 = LEDA_NAMESPACE_NAME::compare_by_angle(v1,v2);
     
     if (res2>0) res2=1;
     if (res2<0) res2=-1;  
     
     return  ((Comparison_result) res2);
   }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_angle_with_x_axis_2<K>::ev_leda_rat_direction;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_compare_slope_2 {

  typedef typename K::Segment_2      Segment_2;
  typedef typename K::Line_2         Line_2;  

public:
  typedef Arity_tag<2> Arity;
  typedef CGAL::Comparison_result    result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_segment;
  static CGAL::event ev_leda_rat_line;  
#endif  

  CGAL::Comparison_result operator()(const Segment_2& s1, const Segment_2& s2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&,const Segment_2&> \
        (Predicate_leda_rat_compare_slope_2::ev_leda_rat_segment, s1, s2);
#endif      
     return (CGAL::Comparison_result)(LEDA_NAMESPACE_NAME::cmp_slopes(s1,s2));
  }
  
  CGAL::Comparison_result operator()(const Line_2& l1, const Line_2& l2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Line_2&> \
        (Predicate_leda_rat_compare_slope_2::ev_leda_rat_line, l1, l2);
#endif   
     return (CGAL::Comparison_result)(LEDA_NAMESPACE_NAME::cmp_slopes(l1,l2));
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_compare_slope_2<K>::ev_leda_rat_segment;
template<class K> CGAL::event  Predicate_leda_rat_compare_slope_2<K>::ev_leda_rat_line;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_distance_to_point_2 {

  typedef typename K::Point_2   Point_2;
  
  typedef typename K::Point_2   __My_Point_2;  

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif   

   bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
   {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
      (Predicate_leda_rat_less_distance_to_point_2::ev_leda_rat_point, p, q, r);
#endif    
     return  (p.cmp_dist(q,r) == -1);
   }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_distance_to_point_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_signed_distance_to_line_2 {

  typedef typename K::Point_2 Point_2;

public:
  typedef Arity_tag< 4 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif  

  bool operator()(const Point_2& p1, const Point_2& p2, 
                  const Point_2& p3, const Point_2& p4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&,const Point_2&> \
      (Predicate_leda_rat_less_signed_distance_to_line_2::ev_leda_rat_point, p1, p2, p3, p4);
#endif   
     // Achtung: CGAL - Kernel doc wird wohl hier geaendert;
     // Fall von gleicher Distanz dann anders zu handlen !!!
  
     return ( ((Comparison_result) LEDA_NAMESPACE_NAME::cmp_signed_dist(p1,p2,p3,p4)) == SMALLER);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_signed_distance_to_line_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_less_rotate_ccw_2 {

  typedef typename K::Point_2  Point_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif   

  bool operator()(const Point_2& p1, const Point_2& p2, 
                  const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_less_rotate_ccw_2::ev_leda_rat_point, p1, p2, p3);
#endif     
     int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
     if (ori == -1) return false;
     
     if (ori == 0) {  // distance comparison
        if (p1.cmp_dist(p2,p3) > 0) return true;
	else return false;
     }
     
     return true;
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_less_rotate_ccw_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_leftturn_2 {

  typedef typename K::Point_2   Point_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif   

  bool operator()(const Point_2& p1, const Point_2& p2, 
                  const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_leftturn_2::ev_leda_rat_point, p1, p2, p3);
#endif   
     return LEDA_NAMESPACE_NAME::left_turn(p1,p2,p3);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_leftturn_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_collinear_2 {

  typedef typename K::Point_2  Point_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif  
  
  bool operator()(const Point_2& p1, const Point_2& p2, 
                  const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_collinear_2::ev_leda_rat_point, p1, p2, p3);
#endif   
     return LEDA_NAMESPACE_NAME::collinear(p1,p2,p3);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_collinear_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_orientation_2 {

  typedef typename K::Point_2   Point_2;

public:
  typedef Arity_tag< 3 >  Arity;
  typedef Orientation     result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif 

  Orientation operator()(const Point_2& p1, const Point_2& p2, 
                         const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_orientation_2::ev_leda_rat_point, p1, p2, p3);
#endif   
     return (Orientation) LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_orientation_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_side_of_oriented_circle_2 {

  typedef typename K::Point_2  Point_2;

public:
  typedef Arity_tag< 4 > Arity;
  typedef Oriented_side  result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif    

  Oriented_side operator()(const Point_2& p1, const Point_2& p2, 
                           const Point_2& p3, const Point_2& p4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_side_of_oriented_circle_2::ev_leda_rat_point, p1, p2, p3, p4);
#endif   
     return (Oriented_side) LEDA_NAMESPACE_NAME::side_of_circle(p1,p2,p3,p4);
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_side_of_oriented_circle_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_side_of_bounded_circle_2 {

  typedef typename K::Point_2  Point_2;

public:
  typedef Bounded_side           result_type;
  typedef Arity_tag< 3 > Arity;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point_point_point_point;
  static CGAL::event ev_leda_rat_point_point_point;  
#endif    

  // precondition: p1,p2,p3 must not be collinear ...
  Bounded_side operator()(const Point_2& p1, const Point_2& p2, 
                          const Point_2& p3, const Point_2& p4) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_side_of_bounded_circle_2::ev_leda_rat_point_point_point_point, p1, p2, p3, p4);
#endif   
     int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
     return (Bounded_side) (LEDA_NAMESPACE_NAME::side_of_circle(p1,p2,p3,p4) * ori);
  }
  
  // p1 and p2 is diameter ...
  Bounded_side operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_side_of_bounded_circle_2::ev_leda_rat_point_point_point, p1, p2, p3);
#endif   
     // compute center ...
     Point_2 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
     // compare distances ...
     return (Bounded_side) m.cmp_dist(p1,p3);
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_side_of_bounded_circle_2<K>::ev_leda_rat_point_point_point_point;
template<class K> CGAL::event  Predicate_leda_rat_side_of_bounded_circle_2<K>::ev_leda_rat_point_point_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_is_horizontal_2 {

  typedef typename K::Line_2      Line_2;
  typedef typename K::Ray_2       Ray_2;
  typedef typename K::Segment_2   Segment_2;    

public:
  typedef Arity_tag< 1 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_segment;    
#endif 

  bool operator()(const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&> \
       (Predicate_leda_rat_is_horizontal_2::ev_leda_rat_line, l);
#endif   
      return l.is_horizontal();
  }
  
  bool operator()(const Ray_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Ray_2&> \
       (Predicate_leda_rat_is_horizontal_2::ev_leda_rat_ray, r);
#endif     
      return r.is_horizontal();  
  }
  
  bool operator()(const Segment_2& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&> \
       (Predicate_leda_rat_is_horizontal_2::ev_leda_rat_segment, s);
#endif   
      return s.is_horizontal();    
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_is_horizontal_2<K>::ev_leda_rat_line;
template<class K> CGAL::event  Predicate_leda_rat_is_horizontal_2<K>::ev_leda_rat_ray;
template<class K> CGAL::event  Predicate_leda_rat_is_horizontal_2<K>::ev_leda_rat_segment;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_is_in_x_range_2 {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Segment_2  Segment_2; 
  
  typedef typename K::Point_2    __My_Point_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event     ev_leda_rat_point_segment;
#endif  
  
  bool operator()(const Point_2& p, const Segment_2& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Point_2&, const Segment_2&> \
        (Predicate_leda_rat_is_in_x_range_2::ev_leda_rat_point_segment, p, s);
#endif   
      // two compare operations ...
      int res1 = __My_Point_2::cmp_x(s.start(), p);
      if (res1 == 0) return true;
      int res2 = __My_Point_2::cmp_x(s.end(), p);
      if (res2 == 0 || res1 != res2) return true;
      
      return false;
  }    

};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event Predicate_leda_rat_is_in_x_range_2<K>::ev_leda_rat_point_segment;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_is_vertical_2 {

  typedef typename K::Line_2    Line_2;
  typedef typename K::Ray_2     Ray_2;
  typedef typename K::Segment_2 Segment_2;    

public:
  typedef Arity_tag< 1 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_segment;
#endif  

  bool operator()(const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Line_2&>(Predicate_leda_rat_is_vertical_2::ev_leda_rat_line, l);
#endif  
      return l.is_vertical();
  }
  
  bool operator()(const Ray_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Ray_2&>(Predicate_leda_rat_is_vertical_2::ev_leda_rat_ray, r);
#endif   
      return r.is_vertical();  
  }
  
  bool operator()(const Segment_2& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
      CGAL::occur<const Segment_2&>(Predicate_leda_rat_is_vertical_2::ev_leda_rat_segment, s);
#endif   
      return s.is_vertical();    
      
      // (test: use float approximation)
      // return (s.dxD() == 0.0);
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_is_vertical_2<K>::ev_leda_rat_line;
template<class K> CGAL::event  Predicate_leda_rat_is_vertical_2<K>::ev_leda_rat_ray;
template<class K> CGAL::event  Predicate_leda_rat_is_vertical_2<K>::ev_leda_rat_segment;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_is_degenerate_2 {

  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::Line_2           Line_2;
  typedef typename K::Ray_2            Ray_2;
  typedef typename K::Segment_2        Segment_2;
  typedef typename K::Triangle_2       Triangle_2; 

public:
  typedef Arity_tag< 1 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle;
  static CGAL::event ev_leda_rat_rectangle;
  static CGAL::event ev_leda_rat_line;
  static CGAL::event ev_leda_rat_ray;
  static CGAL::event ev_leda_rat_segment;
  static CGAL::event ev_leda_rat_triangle;        
#endif   

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_circle, c);
#endif   
      return c.is_degenerate();  
  }
#else  
  bool operator()(const leda_rat_circle& c) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_circle, c);
#endif   
      return c.is_degenerate();  
  }
#endif  

  bool operator()(const Iso_rectangle_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Iso_rectangle_2&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_rectangle, r);
#endif   
      return r.is_degenerate();
  }
  
  bool operator()(const Line_2& l) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_line, l);
#endif   
      return l.seg().is_trivial();
  }  
  
  bool operator()(const Ray_2& r) const 
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Ray_2&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_ray, r);
#endif   
      return (r.point1() == r.point2());  
  }
  
  bool operator()(const Segment_2& s) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_segment, s);
#endif   
      return s.is_trivial();    
  }  
  
  bool operator()(const Triangle_2& t) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&> (Predicate_leda_rat_is_degenerate_2::ev_leda_rat_triangle, t);
#endif   
      return t.is_degenerate();
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_line;
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_ray;
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_segment;
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_circle;
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_rectangle;
template<class K> CGAL::event  Predicate_leda_rat_is_degenerate_2<K>::ev_leda_rat_triangle;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_2 {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Line_2    Line_2;
  typedef typename K::Ray_2     Ray_2;
  typedef typename K::Segment_2 Segment_2;      

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_line_point;
  static CGAL::event ev_leda_rat_ray_point;
  static CGAL::event ev_leda_rat_segment_point;    
#endif   

  bool operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_2::ev_leda_rat_line_point, l, p);
#endif   
      return l.contains(p);
  }  
  
  bool operator()(const Ray_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Ray_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_2::ev_leda_rat_ray_point, r, p);
#endif  
      return r.contains(p);  
  }
  
  bool operator()(const Segment_2& s, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_2::ev_leda_rat_segment_point, s, p);
#endif  
      return s.contains(p);    
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_2<K>::ev_leda_rat_line_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_2<K>::ev_leda_rat_ray_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_2<K>::ev_leda_rat_segment_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_collinear_has_on_2 {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Ray_2      Ray_2;
  typedef typename K::Segment_2  Segment_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_ray_point;
  static CGAL::event ev_leda_rat_segment_point;    
#endif 
  
  bool operator()(const Ray_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Ray_2&,const Point_2&> \
       (Predicate_leda_rat_collinear_has_on_2::ev_leda_rat_ray_point, r, p);
#endif  
      return r.contains(p);  
  }
  
  bool operator()(const Segment_2& s, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Segment_2&,const Point_2&> \
       (Predicate_leda_rat_collinear_has_on_2::ev_leda_rat_segment_point, s, p);
#endif   
      return s.contains(p);    
  }  
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_collinear_has_on_2<K>::ev_leda_rat_ray_point;
template<class K> CGAL::event  Predicate_leda_rat_collinear_has_on_2<K>::ev_leda_rat_segment_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_bounded_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT               FT;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_rectangle_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif   

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_bounded_side_2::ev_leda_rat_circle_point, c, p);
#endif  
     Point_2 center = c.center();
     FT sq      = c.squared_radius();
     FT d       = center.sqr_dist(p);
     
     if (d < sq) return true;
     return false;
      
  }
#else
  bool operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_bounded_side_2::ev_leda_rat_circle_point, c, p);
#endif    
      return c.inside(p);  
  }
#endif  

  bool operator()(const Iso_rectangle_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Iso_rectangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_bounded_side_2::ev_leda_rat_rectangle_point, r, p);
#endif   
      return r.inside(p);
  }

  bool operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_bounded_side_2::ev_leda_rat_triangle_point, t, p);
#endif   
      return t.inside(p);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_bounded_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_bounded_side_2<K>::ev_leda_rat_rectangle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_bounded_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_unbounded_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT               FT;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_rectangle_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif     

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_unbounded_side_2::ev_leda_rat_circle_point, c, p);
#endif  
     Point_2 center = c.center();
     FT      sq      = c.squared_radius();
     FT      d       = center.sqr_dist(p);
     
     if (d > sq) return true;
     return false;
      
  }
#else
  bool operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_unbounded_side_2::ev_leda_rat_circle_point, c, p);
#endif     
      return c.outside(p);  
  }
#endif

  bool operator()(const Iso_rectangle_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Iso_rectangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_unbounded_side_2::ev_leda_rat_rectangle_point, r, p);
#endif  
      return r.outside(p);
  }

  bool operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_unbounded_side_2::ev_leda_rat_triangle_point, t, p);
#endif   
      return t.outside(p);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_unbounded_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_unbounded_side_2<K>::ev_leda_rat_rectangle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_unbounded_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_boundary_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT               FT;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_rectangle_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif  

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_boundary_2::ev_leda_rat_circle_point, c, p);
#endif  
     Point_2 center = c.center();
     FT      sq      = c.squared_radius();
     FT      d       = center.sqr_dist(p);
     
     if (d == sq) return true;
     return false;
  }
#else
  bool operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_boundary_2::ev_leda_rat_circle_point, c, p);
#endif   
      return c.contains(p);  
  }
#endif  

  bool operator()(const Iso_rectangle_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Iso_rectangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_boundary_2::ev_leda_rat_rectangle_point, r, p);
#endif   
      return r.contains(p);
  }

  bool operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_boundary_2::ev_leda_rat_triangle_point, t, p);
#endif  
      return t.on_boundary(p);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_boundary_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_boundary_2<K>::ev_leda_rat_rectangle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_boundary_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_positive_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Line_2           Line_2;
  typedef typename K::FT               FT;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_line_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif  

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_positive_side_2::ev_leda_rat_circle_point, c, p);
#endif   
     Point_2           center = c.center();
     CGAL::Orientation ori    = c.orientation();  
     FT                sq     = c.squared_radius();
     FT                d      = center.sqr_dist(p);
     
     bool res;
     if (d < sq) res = true;
     else res = false;
     
     // test ...
     if (ori == CLOCKWISE) res = !res;
     
     return res; 
  }
#else
  bool operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_positive_side_2::ev_leda_rat_circle_point, c, p);
#endif   
      return (c.side_of(p) == 1);  
  }
#endif

  bool operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_positive_side_2::ev_leda_rat_line_point, l, p);
#endif  

#if (__LEDA__ <= 420)
      return (::orientation(l,p) == 1); 
#else  
      return (l.side_of(p) == 1);   
#endif          
  }

  bool operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_positive_side_2::ev_leda_rat_triangle_point, t, p);
#endif    
      return (t.side_of(p) == 1);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_positive_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_positive_side_2<K>::ev_leda_rat_line_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_positive_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_has_on_negative_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Line_2           Line_2;
  typedef typename K::FT               FT; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_line_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif   

#if defined(CGAL_COMPATIBLE_CIRCLES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_negative_side_2::ev_leda_rat_circle_point, c, p);
#endif   
     Point_2           center = c.center();
     CGAL::Orientation ori    = c.orientation();  
     FT                sq     = c.squared_radius();
     FT                d      = center.sqr_dist(p);
     
     bool res;
     if (d > sq) res = true;
     else res = false;
     
     // testen ...
     if (ori == CLOCKWISE) res = !res;
     
     return res; 
  }
#else
  bool operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_has_on_negative_side_2::ev_leda_rat_circle_point, c, p);
#endif  
      return (c.side_of(p) == -1);  
  }
#endif

  bool operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_negative_side_2::ev_leda_rat_line_point, l, p);
#endif  
#if (__LEDA__ <= 420)  
      return (::orientation(l,p) == -1); 
#else
      return (l.side_of(p) == -1); 
#endif      
  }

  bool operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_has_on_negative_side_2::ev_leda_rat_triangle_point, t, p);
#endif   
      return (t.side_of(p) == -1);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_has_on_negative_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_negative_side_2<K>::ev_leda_rat_line_point;
template<class K> CGAL::event  Predicate_leda_rat_has_on_negative_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_oriented_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Line_2           Line_2;
  typedef typename K::FT               FT; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef Oriented_side           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_line_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif    

#if defined(CGAL_COMPATIBLE_CIRCLES)
  CGAL::Oriented_side operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_oriented_side_2::ev_leda_rat_circle_point, c, p);
#endif  
  
     Point_2           center = c.center();
     CGAL::Orientation ori    = c.orientation();  
     FT                sq     = c.squared_radius();
     FT                d      = center.sqr_dist(p);
     
     if (d == sq) return CGAL::ON_ORIENTED_BOUNDARY;
     
     bool res;
     if (d < sq) res = true;
     else res = false;
     
     // testen ...
     if (ori == CLOCKWISE) res = !res;
     
     if (res) return CGAL::ON_POSITIVE_SIDE;
     return CGAL::ON_NEGATIVE_SIDE; 
  }
#else
  CGAL::Oriented_side operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_oriented_side_2::ev_leda_rat_circle_point, c, p);
#endif  
      return  (CGAL::Oriented_side) c.side_of(p);  
  }
#endif

  CGAL::Oriented_side operator()(const Line_2& l, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Line_2&,const Point_2&> \
       (Predicate_leda_rat_oriented_side_2::ev_leda_rat_line_point, l, p);
#endif
#if (__LEDA__ <= 420)  
      return  (CGAL::Oriented_side) ::orientation(l,p);
#else    
      return  (CGAL::Oriented_side) l.side_of(p);
#endif       
  }

  CGAL::Oriented_side operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_oriented_side_2::ev_leda_rat_triangle_point, t, p);
#endif    
      return (CGAL::Oriented_side) t.side_of(p);
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_oriented_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_oriented_side_2<K>::ev_leda_rat_line_point;
template<class K> CGAL::event  Predicate_leda_rat_oriented_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_bounded_side_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::FT               FT; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef Bounded_side           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_circle_point;
  static CGAL::event ev_leda_rat_rectangle_point;
  static CGAL::event ev_leda_rat_triangle_point;    
#endif    

#if defined(CGAL_COMPATIBLE_CIRCLES)
  CGAL::Bounded_side operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const LEDA_NAMESPACE_NAME::cgal_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_bounded_side_2::ev_leda_rat_circle_point, c, p);
#endif   
     Point_2 center = c.center();
     FT      sq     = c.squared_radius();
     FT      d      = center.sqr_dist(p);
     
     if (d==sq) return CGAL::ON_BOUNDARY;
     if (d <sq) return CGAL::ON_BOUNDED_SIDE;
     return CGAL::ON_UNBOUNDED_SIDE;
  }
#else
  CGAL::Bounded_side operator()(const leda_rat_circle& c, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const leda_rat_circle&,const Point_2&> \
       (Predicate_leda_rat_bounded_side_2::ev_leda_rat_circle_point, c, p);
#endif  
      return (CGAL::Bounded_side) (c.side_of(p) * c.orientation());  
  }
#endif

  CGAL::Bounded_side operator()(const Iso_rectangle_2& r, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Iso_rectangle_2&,const Point_2&> \
       (Predicate_leda_rat_bounded_side_2::ev_leda_rat_rectangle_point, r, p);
#endif   
      if (r.inside(p)) return CGAL::ON_BOUNDED_SIDE;
      if (r.contains(p)) return CGAL::ON_BOUNDARY;
      return CGAL::ON_UNBOUNDED_SIDE;
  }

  CGAL::Bounded_side operator()(const Triangle_2& t, const Point_2& p) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Triangle_2&,const Point_2&> \
       (Predicate_leda_rat_bounded_side_2::ev_leda_rat_triangle_point, t, p);
#endif  
      LEDA_NAMESPACE_NAME::region_kind rk = t.region_of(p);

      if (rk == LEDA_NAMESPACE_NAME::BOUNDED_REGION) return CGAL::ON_BOUNDED_SIDE;
      if (rk == LEDA_NAMESPACE_NAME::ON_REGION)      return CGAL::ON_BOUNDARY;
      return CGAL::ON_UNBOUNDED_SIDE;      
  } 
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_bounded_side_2<K>::ev_leda_rat_circle_point;
template<class K> CGAL::event  Predicate_leda_rat_bounded_side_2<K>::ev_leda_rat_rectangle_point;
template<class K> CGAL::event  Predicate_leda_rat_bounded_side_2<K>::ev_leda_rat_triangle_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_are_ordered_along_line_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif  

  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_are_ordered_along_line_2::ev_leda_rat_point, p, q, r);
#endif    
      return Segment_2(p,r).contains(q);  
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_are_ordered_along_line_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_are_strictly_ordered_along_line_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif  

  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_are_strictly_ordered_along_line_2::ev_leda_rat_point, p, q, r);
#endif  
      return (Segment_2(p,r).contains(q) && ( q != p ) && ( q != r ));  
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_are_strictly_ordered_along_line_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_collinear_are_ordered_along_line_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;

#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif 

  // prec.: collinearity
  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Point_2&,const Point_2&,const Point_2&> \
       (Predicate_leda_rat_collinear_are_ordered_along_line_2::ev_leda_rat_point, p, q, r);
#endif  
      return Segment_2(p,r).contains(q);  
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_collinear_are_ordered_along_line_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_collinear_are_strictly_ordered_along_line_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_point;
#endif   

  // prec.: collinearity
  bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
  {
      return (Segment_2(p,r).contains(q) && ( q != p ) && ( q != r ));  
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_collinear_are_strictly_ordered_along_line_2<K>::ev_leda_rat_point;
#endif

// ------------------------------------------------------------------------------------------------------------

template<class K> 
class Predicate_leda_rat_counterclockwise_in_between_2 {

  typedef typename K::Direction_2   Direction_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef bool           result_type;
  
#if defined(CGAL_GEOMETRY_EVENTS)
  static CGAL::event ev_leda_rat_direction;
#endif   

  bool operator()(const Direction_2& d, 
                  const Direction_2& d1, 
		  const Direction_2& d2) const
  {
#if defined(CGAL_GEOMETRY_EVENTS)
     CGAL::occur<const Direction_2&,const Direction_2&,const Direction_2&>(Predicate_leda_rat_counterclockwise_in_between_2::ev_leda_rat_direction, d, d1, d2);
#endif   
    // we compare angles with the positive x- axis ...
    
    int w1 = LEDA_NAMESPACE_NAME::compare_by_angle(d1.get_vector(),d.get_vector());
    int w2 = LEDA_NAMESPACE_NAME::compare_by_angle(d.get_vector(),d2.get_vector());
    int w3 = LEDA_NAMESPACE_NAME::compare_by_angle(d1.get_vector(),d2.get_vector());
 
/*    
    std::cout << "d:" << d << "\n";
    std::cout << "d1:" << d1 << "\n";
    std::cout << "d2:" << d2 << "\n";       
    std::cout << w1 << " " << w2 << " " << w3 << "\n";
*/    
    // attention - LEDA compare_by_angle can return 2 !!!   
    
    // special cases:
    // cases that d==d1 || d==d2
    if (w1==0) return false;
    if (w2==0) return false;
    
    // d != d1 && d!= d2
    if (w3==0) return true;
    
    // end special cases ...
    
    if (w1 < 0) { // we meet first d1, then d ...
      if (w3 > 0) return true;// ccw: d1->d->d2 
      else {
        if (w2 < 0) return true;
	return false;
      }
    }
    else { // we meet first d, then d1 ...
      if (w3 < 0) return false;
      else { 
        if (w2 < 0) return true;
	return false;      
      }
    }
    
  }
};

#if defined(CGAL_GEOMETRY_EVENTS)
template<class K> CGAL::event  Predicate_leda_rat_counterclockwise_in_between_2<K>::ev_leda_rat_direction;
#endif

// ------------------------------------------------------------------------------------------------------------
// to do: add events ?

template<class K> 
class Predicate_leda_rat_do_intersect_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    
  typedef typename K::Line_2           Line_2;
  typedef typename K::Ray_2            Ray_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;   
  
  typedef typename K::Point_2          __My_Point_2;             
  
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;

  // points ...
 
  bool operator()(const Point_2& p, const Point_2& p2) const
  { return (__My_Point_2::cmp_xy(p,p2)==0); }
  
  bool operator()(const Point_2& p, const Line_2& l) const
  { return l.contains(p); }  
  
  bool operator()(const Line_2& l, const Point_2& p) const
  { return l.contains(p); }  
  
  bool operator()(const Point_2& p, const Ray_2& r) const
  { return r.contains(p); }  
  
  bool operator()(const Ray_2& r, const Point_2& p) const
  { return r.contains(p); }      
  
  bool operator()(const Point_2& p, const Segment_2& s) const
  { return s.contains(p); }  
  
  bool operator()(const Segment_2& s, const Point_2& p) const
  { return s.contains(p); }    

  bool operator()(const Point_2& p, const Triangle_2& t) const
  { return t.contains(p); }  
  
  bool operator()(const Triangle_2& t, const Point_2& p) const
  { return t.contains(p); }  
  
  bool operator()(const Point_2& p, const Iso_rectangle_2& r) const
  { return r.inside_or_contains(p); }  
  
  bool operator()(const Iso_rectangle_2& r, const Point_2& p) const
  { return r.inside_or_contains(p); }   
  
  // line ...
  
  bool operator()(const Line_2& l, const Line_2& l2) const
  {
    Point_2 inter;
    return l.intersection(l2, inter);
  }
  
  bool operator()(const Line_2& l, const Ray_2& r) const
  {
#if (__LEDA__ <= 420)
    int o1 = ::orientation(l, r.point1());
    int o2 = ::orientation(l, r.point2());
#else  
    int o1 = l.side_of(r.point1());
    int o2 = l.side_of(r.point2());
#endif    
    
    if (o1 != o2 || o1==0 || o2==0) return true;
    
    // dist. comparison ...
    // point2 must be nearer !!!
    int cmp = LEDA_NAMESPACE_NAME::cmp_signed_dist(l.point1(), l.point2(), r.point1(), r.point2());
    
    if (cmp == 1) return true;
    return false;
  } 
  
  bool operator()(const Ray_2& r, const Line_2& l) const        
  { return this->operator()(l,r); }
  
  bool operator()(const Line_2& l, const Segment_2& s) const
  { return l.intersection(s); }
  
  bool operator()(const Segment_2& s, const Line_2& l) const
  { return l.intersection(s); } 
  
  bool operator()(const Line_2& l, const Triangle_2& t) const
  { return t.intersection(l); }
  
  bool operator()(const Triangle_2& t, const Line_2& l) const
  { return t.intersection(l); } 
  
  bool operator()(const Line_2& l, const Iso_rectangle_2& r) const
  {
    leda_list<Point_2> LI = r.intersection(l);
    if (! LI.empty()) return true;
    return false;  
  }
  
  bool operator()(const Iso_rectangle_2& r, const Line_2& l) const
  { return this->operator()(l,r); } 
  
  // ray ...
  bool operator()(const Ray_2& r, const Ray_2& r2) const
  { 
    Point_2 inter;
    return r.intersection(r2,inter);
  }
  
  bool operator()(const Ray_2& r, const Segment_2& s) const
  { 
    Point_2 inter;
    return r.intersection(s,inter);    
  }
  
  bool operator()(const Segment_2& s, const Ray_2& r) const
  { 
    Point_2 inter;
    return r.intersection(s,inter);    
  }  

  bool operator()(const Ray_2& r, const Triangle_2& t) const
  { 
    // do we intersect one of the sides of t ?
    Point_2 inter;
    
    Segment_2 s1(t.point1(),t.point2());
    if (r.intersection(s1,inter)) return true;
    
    Segment_2 s2(t.point2(),t.point3());
    if (r.intersection(s2,inter)) return true;
    
    Segment_2 s3(t.point3(),t.point1());
    if (r.intersection(s3,inter)) return true;    
        
    return false;
  }  
  
  bool operator()(const Triangle_2& t, const Ray_2& r) const
  {    
    return this->operator()(r,t);
  } 
  
  bool operator()(const Ray_2& r, const Iso_rectangle_2& rect) const
  { 
    Line_2 l(r.point1(), r.point2());
    leda_list<Point_2> LI = rect.intersection(l);
    if (LI.empty()) return false;
    Point_2 p;
    forall(p,LI) if (r.contains(p)) return true;
    return false;
  }  
  
  bool operator()(const Iso_rectangle_2& rect, const Ray_2& r) const
  {    
    return this->operator()(r,rect);
  }
  
  // segment ... 
  bool operator()(const Segment_2& s, const Segment_2& s2) const
  { return s.intersection(s2); } 
  
  bool operator()(const Segment_2& s, const Triangle_2& t) const
  { return t.intersection(s);  }
  
  bool operator()(const Triangle_2& t, const Segment_2& s) const
  { return t.intersection(s);  } 
  
  bool operator()(const Segment_2& s, const Iso_rectangle_2& r) const  
  {
    leda_list<Point_2> LI = r.intersection(s);
    if (LI.empty()) return false;
    return true;
  }
  
  bool operator()(const Iso_rectangle_2& r, const Segment_2& s) const
  {
    leda_list<Point_2> LI = r.intersection(s);
    if (LI.empty()) return false;
    return true;  
  }   
  
  // spaeter effizienter implementieren ...
  // triangle ...
  bool operator()(const Triangle_2& a, const Triangle_2& b) const
  {
    Point_2 a1 = a.point1(), a2 = a.point2(), a3 = a.point3();
    Point_2 b1 = b.point1(), b2 = b.point2(), b3 = b.point3();  
    
    if (a.contains(b1) || a.contains(b2) || a.contains(b3)) return true;
    if (b.contains(a1) || b.contains(a2) || b.contains(a3)) return true;  
    
    // Segment - intersection tests ...
    Segment_2 s1(a1,a2), s2(a2,a3), s3(a3,a1);
    Segment_2 t1(b1,b2), t2(b2,b3), t3(b3,b1);    
    
    if (s1.intersection(t1) || s1.intersection(t2) || s1.intersection(t3)) return true;
    if (s2.intersection(t1) || s2.intersection(t2) || s2.intersection(t3)) return true; 
    // this is probably not needed ?
    if (s3.intersection(t1) || s3.intersection(t2) || s3.intersection(t3)) return true;        
    
    return false;   
  }
  
  bool operator()(const Triangle_2& t, const Iso_rectangle_2& r) const
  {
    // get the sides (segments) of t
    Point_2 a1 = t.point1(), a2 = t.point2(), a3 = t.point3();
    
    if (r.inside_or_contains(a1) || r.inside_or_contains(a2) || r.inside_or_contains(a3)) return true;
    
    Segment_2 s1(a1,a2), s2(a2,a3), s3(a3,a1);
    leda_list<Point_2> LI = r.intersection(s1);
    if (! LI.empty()) return true;
    LI = r.intersection(s2);
    if (! LI.empty()) return true;
    LI = r.intersection(s3);
    if (! LI.empty()) return true;
        
    return false;
  } 
  
  bool operator()(const Iso_rectangle_2& r, const Triangle_2& t) const
  {
    return this->operator()(t,r);
  }   
    
  // rectangle ...
  bool operator()(const Iso_rectangle_2& r, const Iso_rectangle_2& r2) const
  {
    leda_list<Iso_rectangle_2> LI = r.intersection(r2);
    
    return (! LI.empty()); 
  }        
};

// ------------------------------------------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif




