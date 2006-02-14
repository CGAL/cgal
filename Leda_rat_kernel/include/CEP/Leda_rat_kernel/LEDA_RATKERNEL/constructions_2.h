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
#include <LEDA/interval.h>


#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

/*
Todo:
provide complete bbox construction (switch to intervals)
*/

CGAL_BEGIN_NAMESPACE

template <class K>
class Construct_leda_rat_bbox
{
  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
  typedef typename K::Triangle_2       Triangle_2;
  typedef typename K::Circle_2         Circle_2;
  
public:
  typedef Bbox_2          result_type;
  typedef Arity_tag< 1 >   Arity;

  Bbox_2 get_point_bbox(const Point_2& p) const
  { 
      leda_rational x = p.xcoord();
      leda_rational y = p.ycoord();
      
      LEDA_NAMESPACE_NAME::interval xi(x);
      LEDA_NAMESPACE_NAME::interval yi(y);
      
      return Bbox_2(xi.lower_bound(),yi.lower_bound(),xi.upper_bound(),yi.upper_bound());  
  }
  
  Bbox_2
  operator()(const Point_2& p) const
  { return get_point_bbox(p); }  
  
  Bbox_2
  operator()(const Segment_2& s) const
  { return get_point_bbox(s.source()) + get_point_bbox(s.target()); }

  Bbox_2
  operator()(const Triangle_2& t) const
  { return get_point_bbox(t.point1()) + get_point_bbox(t.point2()) + get_point_bbox(t.point3()); }

  Bbox_2
  operator()(const Iso_rectangle_2& r) const
  { return get_point_bbox(r.upper_left()) + get_point_bbox(r.lower_right()); }

  Bbox_2
  operator()(const Circle_2& c) const
  {
    leda_rat_point p = c.center();
    leda_rational  r = c.sqr_radius();
    
    LEDA_NAMESPACE_NAME::interval i(r);
    i = LEDA_NAMESPACE_NAME::sqrt(i);
    
    // now build the bbox ...
    leda_rational x = p.xcoord();
    leda_rational y = p.ycoord();
      
    LEDA_NAMESPACE_NAME::interval xi(x);
    LEDA_NAMESPACE_NAME::interval yi(y);
    
    LEDA_NAMESPACE_NAME::interval xmin = xi - i;
    LEDA_NAMESPACE_NAME::interval ymin = yi - i;
    LEDA_NAMESPACE_NAME::interval xmax = xi + i;
    LEDA_NAMESPACE_NAME::interval ymax = yi + i;
      
    return Bbox_2(xmin.lower_bound(),ymin.lower_bound(),xmax.upper_bound(),ymax.upper_bound());     
  }
  
};


template<class K>
class Construct_leda_rat_point {

  typedef typename K::Point_2    Point_2;
  typedef typename K::FT         FT;
  typedef typename K::RT         RT;    

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_2           result_type;
  
  //----------------------------------------------------------------------------------
  //undocumented
  Point_2 operator()() const
  { 
   Point_2 p;
   return p;
  } 
  
  Point_2 operator()(const FT& x, const FT& y) const
  {     
    Point_2 p(x,y);
    return p;
  }  
  
  Point_2 operator()(const RT& x, const RT& y, const RT& w) const
  { 
    Point_2 p(x,y,w);
    return p;
  } 
  
  //---------------------------------------------------------------------------------- 
  Point_2 operator()(const CGAL::Origin& orig) const
  {  
    return Point_2(0,0,1); 
  }
};


template<class K>
class Construct_leda_rat_vector {

  typedef typename K::Point_2    Point_2;
  typedef typename K::Vector_2   Vector_2; 
  typedef typename K::FT         FT;
  typedef typename K::RT         RT;       

public:
  typedef Vector_2           result_type;
  typedef Arity_tag< 2 > Arity;

  //----------------------------------------------------------------------------------   
  //undocumented
  Vector_2 operator()() const
  { 
   Vector_2 v(2);
   return v;
  }  
  
  Vector_2 operator()(const FT& x, const FT& y) const
  {   
    Vector_2 v(x,y);
    return v;
  }  
  
  Vector_2 operator()(const RT& x, const RT& y, const RT& w) const
  {   
    Vector_2 v(x,y,w);
    return v;
  } 
  
  //----------------------------------------------------------------------------------   

  Vector_2 operator()(const Point_2& a, const Point_2& b) const
  {   
    return b-a; 
  }
   
  Vector_2 operator()(const CGAL::Null_vector& nv) const
  {     
    return Vector_2(2); 
  }
};


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
  
  //undocumented
  Direction_2 operator()() const
  { 
   Direction_2 d(2);
   return d;
  } 
  
  Direction_2 operator()(const FT& x, const FT& y) const
  { 
    Direction_2 d(x,y);
    return d;
  }        

  Direction_2 operator()(const Vector_2& v) const
  {
    return Direction_2(v); 
  }
    
  Direction_2 operator()(const Line_2& l) const
  {     
      Vector_2 v = l.point2()-l.point1();
      return Direction_2(v); 
  }
    
  Direction_2 operator()(const Ray_2& r) const
  {   
     Vector_2 v = r.point2()-r.point1(); 
     return Direction_2(v);
  }
    
  Direction_2 operator()(const Segment_2& s) const
  {  
    Vector_2 v = s.end()-s.start();
    return Direction_2(v); 
  }            
};

template<class K>
class Construct_leda_rat_segment {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Segment_2 Segment_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Segment_2      result_type;
  
  //undocumented
  Segment_2 operator()() const
  { 
   Segment_2 s;
   return s;
  }   

  Segment_2 operator()(const Point_2& p1, const Point_2& p2) const
  {   
    return Segment_2(p1,p2); 
  }
}; 

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
  

  //undocumented
  Line_2 operator()() const
  { 
   Line_2 l;
   return l;
  }    

  Line_2 operator()(const RT& a, const RT& b, const RT& c) const
  { 
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
    return Line_2(p1,p2); 
  }

  Line_2 operator()(const Point_2& p1, const Direction_2& dir) const
  { 
    return Line_2(p1, dir.get_vector()); 
  }

  Line_2 operator()(const Segment_2& s) const
  { 
    return Line_2(s.start(), s.end()); 
  }

  Line_2 operator()(const Ray_2& r) const
  {   
    return Line_2(r); 
  }
};

template<class K>
class Construct_leda_rat_ray {

  typedef typename K::Point_2   Point_2;
  typedef typename K::Direction_2  Direction_2;
  typedef typename K::Ray_2     Ray_2; 

public:
  typedef Arity_tag< 2 >  Arity;
  typedef Ray_2           result_type;
  
  //undocumented
  Ray_2 operator()() const
  { 
   Ray_2 r;
   return r;
  }    

  Ray_2 operator()(const Point_2& p1, const Point_2& p2) const
  { 
    return Ray_2(p1,p2); 
  }
    
  Ray_2 operator()(const Point_2& p1, const Direction_2& dir) const
  {  
    return Ray_2(p1, dir.get_vector()); 
  }    
};


#if defined(CGAL_COMPATIBLE_CIRCLES)
template<class K>
class Construct_leda_rat_circle {

  typedef typename K::Point_2   Point_2;
  typedef typename K::FT        FT;

public:
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle           result_type;
  typedef Arity_tag< 3 >          Arity;  

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
      return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,sqrad,ori);
    }    

    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, 
                                                    const Point_2& p2, 
						    const Point_2& p3) const
    {
     int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3);
     leda_rat_circle C(p1,p2,p3);
     Point_2  center = C.center();
     
     CGAL::Orientation cg_ori;
     switch(ori){
        case -1: { cg_ori = CGAL::RIGHT_TURN; break; }
	case  0: { cg_ori = CGAL::COLLINEAR; break; }
	case  1: { cg_ori = CGAL::LEFT_TURN; break; }
     }
    
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(center,p1,cg_ori);
    }
    
    // diameter version ...    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, const Point_2& p2, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {   
     Point_2 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(m,p1,ori);
    }
    
    LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const Point_2& p1, 
                                                    CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {  
     return LEDA_NAMESPACE_NAME::cgal_rat_circle(p1,ori);
    }        
};

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

  //undocumented
  leda_rat_circle operator()() const
  { 
   leda_rat_circle  c;
   return c;
  }   

    leda_rat_circle operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
    {
     return leda_rat_circle(p1,p2,p3);
    }
    
    leda_rat_circle operator()(const Point_2& p1, const Point_2& p2, CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
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
     return leda_rat_circle(p1,p1);
    }        
};

#endif

template<class K>
class Construct_leda_rat_triangle {

  typedef typename K::Point_2     Point_2;
  typedef typename K::Triangle_2  Triangle_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef Triangle_2     result_type;
  
  //undocumented
  Triangle_2 operator()() const
  { 
   Triangle_2  t;
   return t;
  }  

  Triangle_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
   return Triangle_2(p1,p2,p3);
  }
};

template<class K>
class Construct_leda_rat_rectangle {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;

public:
  typedef Arity_tag< 2 >    Arity;
  typedef Iso_rectangle_2   result_type;

  //undocumented
  Iso_rectangle_2 operator()() const
  { 
   Iso_rectangle_2 r;
   return r;
  }  

  Iso_rectangle_2 operator()(const Point_2& p1, const Point_2& p2) const
  {  
   return Iso_rectangle_2(p1,p2);
  } 
  
  // new
  Iso_rectangle_2 operator()(const Point_2& xmin, const Point_2& xmax,
                             const Point_2& ymin, const Point_2& ymax)
  {
   return Iso_rectangle_2(xmin.xcoord(),ymin.ycoord(),  xmax.xcoord(),ymax.ycoord());
  }   
};

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

  Vector_2 operator()(const Vector_2& v, const RT& scale) const
  {
   return scale * v;
  }      

  Vector_2 operator()(const Vector_2& v, const CGAL::Quotient<RT>& scale) const
  {
   FT fkt(scale.numerator(), scale.denominator());
   return fkt * v;
  } 
};

template<class K>
class Construct_leda_rat_translated_point_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Vector_2      Vector_2;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type;

  Point_2 operator()(const Point_2& p, const Vector_2& v) const
  {
    return p.translate(v);
  }    
};

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

  Point_2 operator()(const Line_2& l, int i = 0) const
  {
    Segment_2 s = l.seg();
    Vector_2  v = s.to_vector();
      
    return s.start() + (RT(i) * v);   
  }

  Point_2 operator()(const Ray_2& r, int i) const
  {
    if (i==0) return r.source();
    return r.point2(); // return a point different from the source ...
  }
    
  Point_2 operator()(const Segment_2& s, int i) const
  {
    i = i % 2;
    if (i==0) return s.start();
      
    return s.end();
  }        
    
};


template<class K> 
class Construct_leda_rat_projected_point_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Line_2        Line_2;
  typedef typename K::Segment_2     Segment_2; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type;

  // orthogonal projection of p onto l
  Point_2 operator()(const Line_2& l, const Point_2& p) const
  {
    Segment_2 s = l.perpendicular(p);
    return s.end();
  }
};


template<class K>
class Construct_leda_rat_vertex_2 {

  typedef typename K::Point_2       Point_2;
  typedef typename K::Triangle_2    Triangle_2;
  typedef typename K::Segment_2     Segment_2; 
  typedef typename K::Iso_rectangle_2  Iso_rectangle_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2        result_type; 

  Point_2 operator()(const Segment_2& s, int i) const
  { 
    i = i % 2;
      
    if (i==0) return s.start();
    return s.end();
  }
    
  Point_2 operator()(const Iso_rectangle_2& r, int i) const
  {
    i = i % 4;
    return r[i+1];
  }
    
  Point_2 operator()(const Triangle_2& t, int i) const
  {
    i = i % 3;
    return t[i+1];
  }    
};

template<class K>
class Construct_leda_rat_supporting_line_2 {

  typedef typename K::Line_2        Line_2;
  typedef typename K::Ray_2         Ray_2;
  typedef typename K::Segment_2     Segment_2; 

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_2  result_type;

  Line_2 operator()(const Ray_2& r) const
  {
     return Line_2(r.point1(), r.point2());
  }
    
  Line_2 operator()(const Segment_2& s) const
  {
     return Line_2(s.start(), s.end());
  }
};

template<class K>
class Construct_leda_rat_perpendicular_direction_2 {

  typedef typename K::Direction_2        Direction_2;
  typedef typename K::Vector_2           Vector_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Direction_2    result_type;

  Direction_2 operator()(const Direction_2& d, CGAL::Orientation ori) const
  {
      Vector_2 v = d.get_vector();
  
      if (ori == CGAL::COUNTERCLOCKWISE)
      {
#if (__LEDA__ <= 420)
        v = v.rotate90();
#else      
        v = v.rotate90(1);
#endif	
      }
      else { // clockwise ...
#if (__LEDA__ <= 420)
        Vector_2 hlp = v.rotate90();
	hlp = hlp.rotate90();
	hlp = hlp.rotate90();
	v = hlp;
#else      
        v = v.rotate90(-1);
#endif	
      }
      
      // collinear is not allowed
      return Direction_2(v);
  }
};

template<class K>
class Construct_leda_rat_perpendicular_vector_2 {

  typedef typename K::Vector_2           Vector_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_2       result_type;

  Vector_2 operator()(const Vector_2& v, CGAL::Orientation ori) const
  {

      if (ori == CGAL::COUNTERCLOCKWISE)
      {
#if (__LEDA__ <= 420)
        return v.rotate90();
#else      
        return v.rotate90(1);
#endif	
      }
      else { // clockwise ...
#if (__LEDA__ <= 420)
        Vector_2 hlp = v.rotate90();
	hlp = hlp.rotate90();
	hlp = hlp.rotate90();
	return hlp;
#else      
        return v.rotate90(-1);
#endif	
      }
      
      // collinear is not allowed
      return v;
  }
};

template<class K>
class Construct_leda_rat_perpendicular_line_2 {

  typedef typename K::Line_2           Line_2;
  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Line_2         result_type;

  Line_2 operator()(const Line_2& l, const Point_2& p) const
  { 
    // construct perp. line through p; rotation ccw by 90 degrees
#if (__LEDA__ <= 420)  
    int ori = ::orientation(l,p);
#else    
    int ori = l.side_of(p);
#endif    
      
    if (ori == 0) { // special case: collinear
#if (__LEDA__ <= 420)  
      return l.rotate90(p);
#else    
      return l.rotate90(p,1);
#endif      
    }
      
    Segment_2 s = l.perpendicular(p);
      
    if (ori == +1){
      return Line_2(s.reverse());
    }
      
    return Line_2(s);    
  }
};

template<class K>
class Construct_leda_rat_midpoint_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Segment_2        Segment_2;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2           result_type;

  Point_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
     return LEDA_NAMESPACE_NAME::midpoint(p1,p2);
  }
};    

template<class K>
class Construct_leda_rat_center_2 {

  typedef typename K::Point_2          Point_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_2        result_type;

#if defined(CGAL_COMPATIBLE_CIRCLES)
  Point_2 operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& C) const
  { 
     return C.center();
  }
#else
  Point_2 operator()(const leda_rat_circle& C) const
  {
     return C.center();
  }
#endif    
};  

template<class K>
class Construct_leda_rat_centroid_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::FT               FT;

public:
  typedef Point_2           result_type;
  typedef Arity_tag< 3 >    Arity;  

  Point_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
     // sum up coordinates, divide by 3
     FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() )/ FT(3.0);
     FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() )/ FT(3.0);
       
     return Point_2(x,y);
  }

  Point_2 operator()(const Point_2& p1, const Point_2& p2, 
                     const Point_2& p3, const Point_2& p4) const
  {
     // sum up coordinates, divide by 4
     FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() + p4.xcoord() )/ FT(4.0);
     FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() + p4.ycoord() )/ FT(4.0);
       
     return Point_2(x,y);       
  }
}; 

template<class K>
class Construct_leda_rat_circumcenter_2 {

  typedef typename K::Point_2          Point_2;

public:
  typedef Arity_tag< 3 > Arity;
  typedef Point_2        result_type;

  Point_2 operator()(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
    leda_rat_circle C(p1,p2,p3);
    return C.center();
  }
};

template<class K>
class Construct_leda_rat_bisector_2 {

  typedef typename K::Point_2          Point_2;
  typedef typename K::Line_2           Line_2;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Line_2         result_type;

  Line_2 operator()(const Point_2& p1, const Point_2& p2) const
  {
     return LEDA_NAMESPACE_NAME::p_bisector(p1,p2);
  }
};

template<class K>
class Construct_leda_rat_opposite_direction_2 {

  typedef typename K::Vector_2         Vector_2;
  typedef typename K::Direction_2      Direction_2;

public:

  typedef Arity_tag< 1 > Arity;
  typedef Direction_2    result_type;

  Direction_2 operator()(const Direction_2& dir) const
  {   
      Vector_2 v = dir.get_vector(); 
      return Direction_2(-v);
  }
};

template<class K>
class Construct_leda_rat_opposite_segment_2 {

  typedef typename K::Segment_2  Segment_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Segment_2      result_type;

  Segment_2 operator()(const Segment_2& s) const
  { 
    return s.reverse(); 
  }
};

template<class K>
class Construct_leda_rat_opposite_ray_2 {

  typedef typename K::Ray_2  Ray_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Ray_2          result_type;

  Ray_2 operator()(const Ray_2& r) const
  {  
    return r.reverse(); 
  }
};

template<class K>
class Construct_leda_rat_opposite_line_2 {

  typedef typename K::Line_2  Line_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_2         result_type;

  Line_2 operator()(const Line_2& l) const
  { 
    return l.reverse(); 
  }
};

template<class K>
class Construct_leda_rat_opposite_triangle_2 {

  typedef typename K::Triangle_2   Triangle_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Triangle_2     result_type;

  Triangle_2 operator()(const Triangle_2& t) const
  {
    return t.reverse(); 
  }
};

template<class K>
class Construct_leda_rat_opposite_circle_2 {

  typedef typename K::Point_2  Point_2;
  typedef typename K::FT       FT;

public:
  typedef Arity_tag< 1 > Arity;
  
#if defined(CGAL_COMPATIBLE_CIRCLES)
  typedef LEDA_NAMESPACE_NAME::cgal_rat_circle         result_type;

  LEDA_NAMESPACE_NAME::cgal_rat_circle operator()(const LEDA_NAMESPACE_NAME::cgal_rat_circle& c) const
  {  
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
    return c.reverse(); 
  }
#endif  
};

template<class K>
class Construct_leda_rat_opposite_vector_2 {

  typedef typename K::Vector_2 Vector_2;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Vector_2       result_type;

  Vector_2 operator()(const Vector_2& v) const
  { 
    return -v; 
  }
};

CGAL_END_NAMESPACE

#endif




