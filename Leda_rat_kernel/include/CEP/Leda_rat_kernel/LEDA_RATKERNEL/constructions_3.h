#ifndef CEP_LEDA_RAT_CONSTRUCTIONS_3_H
#define CEP_LEDA_RAT_CONSTRUCTIONS_3_H

// LEDA 3d rational kernel construction objects ...

// 3d constructions ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/d3_rat_support_functions.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/support_functions.h>

#include <LEDA/interval.h>

#if !defined(LEDA_NAMESPACE_NAME)
#define LEDA_NAMESPACE_NAME
#endif

/*
Todo:
provide complete bbox construction
*/

CGAL_BEGIN_NAMESPACE

template <class K>
class Construct_leda_d3_rat_bbox
{
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Sphere_3         Sphere_3;
public:
    typedef Bbox_3          result_type;
    typedef Arity_tag< 1 >   Arity;
    
    Bbox_3 get_point_bbox(const Point_3& p) const
    { 
      leda_rational x = p.xcoord();
      leda_rational y = p.ycoord();
      leda_rational z = p.zcoord();
      
      LEDA_NAMESPACE_NAME::interval xi(x);
      LEDA_NAMESPACE_NAME::interval yi(y);
      LEDA_NAMESPACE_NAME::interval zi(z);
      
      return Bbox_3(xi.lower_bound(),yi.lower_bound(),zi.lower_bound,
                    xi.upper_bound(),yi.upper_bound(),zi.upper_bound);    
    }    

    Bbox_3
    operator()(const Point_3& p) const
    { return get_point_bbox(p); }

    Bbox_3
    operator()(const Segment_3& s) const
    { return get_point_bbox(s.source()) + get_point_bbox(s.target()); }

    
    Bbox_3
    operator()(const Triangle_3& t) const
    { return get_point_bbox(t.point1()) + get_point_bbox(t.point2()) + get_point_bbox(t.point3()); } 

    Bbox_3
    operator()(const Iso_cuboid_3& r) const
    {
      leda_d3_rat_point min(r.xmin(), r.ymin(), r.zmin());
      leda_d3_rat_point max(r.xmax(), r.ymax(), r.zmax());
      
      return get_point_bbox(min) + get_point_bbox(max);      
    }

    Bbox_3
    operator()(const Tetrahedron_3& t) const
    {
      return get_point_bbox(t.point1()) + get_point_bbox(t.point2()) + 
             get_point_bbox(t.point3()) + get_point_bbox(t.point4()); 
    }

    Bbox_3
    operator()(const Sphere_3& s) const
    { 
    leda_d3_rat_point p = s.center();
    leda_rational     r = s.sqr_radius();
    
    LEDA_NAMESPACE_NAME::interval i(r);
    i = LEDA_NAMESPACE_NAME::sqrt(i);
    
    // now build the bbox ...
    leda_rational x = p.xcoord();
    leda_rational y = p.ycoord();
    leda_rational z = p.zcoord();
      
    LEDA_NAMESPACE_NAME::interval xi(x);
    LEDA_NAMESPACE_NAME::interval yi(y);
    LEDA_NAMESPACE_NAME::interval zi(z);
    
    LEDA_NAMESPACE_NAME::interval xmin = xi - i;
    LEDA_NAMESPACE_NAME::interval ymin = yi - i;
    LEDA_NAMESPACE_NAME::interval zmin = zi - i;
    LEDA_NAMESPACE_NAME::interval xmax = xi + i;
    LEDA_NAMESPACE_NAME::interval ymax = yi + i;
    LEDA_NAMESPACE_NAME::interval zmax = zi + i;
      
    return Bbox_2(xmin.lower_bound(),ymin.lower_bound(),zmin.lower_bound(),
                  xmax.upper_bound(),ymax.upper_bound(),zmax.upper_bound());        
    }
};



template<class K>
class Construct_leda_d3_rat_point {

  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_3           result_type;
  
  //----------------------------------------------------------------------------------
  //undocumented
  Point_3 operator()() const
  { 
   Point_3 p;
   return p;
  } 
  
  Point_3 operator()(const FT& x, const FT& y, const FT& z) const
  {
    Point_3 p(x,y,x);
    return p;
  }  
  
  Point_3 operator()(const RT& x, const RT& y, const RT& z, const RT& w) const
  {
    Point_3 p(x,y,z,w);
    return p;
  } 
  
  //----------------------------------------------------------------------------------
  Point_3 operator()(const CGAL::Origin& orig) const
  { return Point_3(0,0,0,1); }
};

template<class K>
class Construct_leda_d3_rat_vector {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Vector_3  Vector_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;

public:
  typedef Vector_3           result_type;
  typedef Arity_tag< 2 >            Arity;

  //----------------------------------------------------------------------------------   
  //undocumented
  Vector_3 operator()() const
  { 
   Vector_3 v(3);
   return v;
  }  
  
  Vector_3 operator()(const FT& x, const FT& y, const FT& z) const
  {
    Vector_3 v(x,y,z);
    return v;
  }  
  
  Vector_3 operator()(const RT& x, const RT& y, const RT& z, const RT& w) const
  {
    Vector_3 v(x,y,z,w);
    return v;
  } 
  
  //---------------------------------------------------------------------------------- 
  Vector_3 operator()(const Point_3& a, const Point_3& b) const
  { return b-a; }
   
  Vector_3 operator()(const CGAL::Null_vector& nv) const
  { return Vector_3(3); }
};

// Achtung : fuehre hier typedef aus fuer Construct_direction_of_line_3
// (in der Traits - Klasse)
// welches in 3d Delaunay gebraucht wird ...

template<class K>
class Construct_leda_d3_rat_direction {

  typedef typename K::Direction_3  Direction_3;  
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::Line_3    Line_3;
  typedef typename K::Ray_3     Ray_3;
  typedef typename K::Segment_3 Segment_3; 
  typedef typename K::FT        FT;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Direction_3           result_type;
  
  //undocumented
  Direction_3 operator()() const
  { 
   Direction_3 d(3);
   return d;
  } 
  
  Direction_3 operator()(const FT& x, const FT& y, const FT& z) const
  {
    Direction_3 d(x,y,z);
    return d;
  }          
  
  Direction_3 operator()(const Vector_3& v) const
  { return Direction_3(v); }
    
  Direction_3 operator()(const Line_3& l) const
  { Vector_3 v = l.point2()-l.point1();
    return Direction_3(v); 
  }
    
  Direction_3 operator()(const Ray_3& r) const
  {  Vector_3 v = r.point2()-r.point1(); 
     return Direction_3(v);
  }
    
  Direction_3 operator()(const Segment_3& s) const
  { Vector_3 v = s.target()-s.source();
    return Direction_3(v); 
  }            
};

template<class K>
class Construct_leda_d3_rat_plane {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3;
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::Direction_3  Direction_3;  
  typedef typename K::Line_3    Line_3;
  typedef typename K::Ray_3     Ray_3;
  typedef typename K::Segment_3 Segment_3;   
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;

public:
  typedef Plane_3           result_type;
  typedef Arity_tag< 3 >              Arity;
  
  Plane_3 operator()() const
  {
    Plane_3 p;
    return p;
  }  

  Plane_3 operator()(const RT& a, const RT& b, const RT& c, const RT& d) const
  {
      CGAL_precondition(! (a==0 && b==0 && c==0));

      RT x=0, y=0, z=0;
      FT other;
      RT n, den;
      Point_3 on_plane; // a point on the plane ...

      // compute point (x,y,z) on plane with
      // a*x + b*y + c*z + d = 0
      
      if (c!=0){
        other = FT(-a*x - b*y - d,c);
	n = other.numerator();
	den = other.denominator();
	on_plane = Point_3(x*den,y*den,n,den);
      }
      else {
        if (b!=0){
	  other = FT(-a*x - c*z -d,b);
	  n = other.numerator();
	  den = other.denominator();
	  on_plane = Point_3(x*den,n,z*den,den);	  
	}
	else {
	  other = FT(-b*y - c*z -d,a);
	  n = other.numerator();
	  den = other.denominator();
	  on_plane = Point_3(n,y*den,z*den,den);	  
	}
      }
      
      // normal vector ...
      Vector_3 nv(a,b,c,1);
      
      return Plane_3(on_plane, nv);
  }
    
  Plane_3 operator()(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  {
       return Plane_3(p1,p2,p3);
  }   
    
  Plane_3 operator()(const Point_3& p, const Direction_3& d) const
  {
       // construction with normal vector ...
       return Plane_3(p, d.get_vector());
  }
    
  Plane_3 operator()(const Line_3& l, const Point_3& p) const
  {
       return Plane_3(l.point1(), l.point2(), p);
  }
    
  Plane_3 operator()(const Ray_3& r, const Point_3& p) const
  {
       return Plane_3(r.point1(), r.point2(), p);    
  }
    
  Plane_3 operator()(const Segment_3& s, const Point_3& p) const
  {
       return Plane_3(s.source(), s.target(), p);        
  }            
    
};


template<class K>
class Construct_leda_d3_rat_iso_cuboid {

  typedef typename K::Point_3        Point_3;
  typedef typename K::Iso_cuboid_3   Iso_cuboid_3;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Iso_cuboid_3           result_type;

  Iso_cuboid_3 operator()() const
  {
    Iso_cuboid_3 i;
    return i;
  }

  Iso_cuboid_3 operator()(const Point_3& p, const Point_3& q) const
  {
     return Iso_cuboid_3(p,q);
  }
  
  Iso_cuboid_3 operator()(const Point_3& pxmin, const Point_3& pxmax,
                          const Point_3& pymin, const Point_3& pymax,
			  const Point_3& pzmin, const Point_3& pzmax)
  {
     return Iso_cuboid_3(pxmin,pxmax,pymin,pymax,pzmin,pzmax);
  }
};

template<class K>
class Construct_leda_d3_rat_line {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Line_3    Line_3;
  typedef typename K::Ray_3     Ray_3;
  typedef typename K::Segment_3 Segment_3;   
  typedef typename K::Direction_3 Direction_3; 
  typedef typename K::Vector_3  Vector_3;    

public:
  typedef Line_3           result_type;  
  typedef Arity_tag< 2 >   Arity;

  Line_3 operator()() const
  {
    Line_3 l;
    return l;
  }

  Line_3 operator()(const Point_3& p, const Point_3& q) const
  {
    return Line_3(p,q);
  }
     
  Line_3 operator()(const Point_3& p, const Direction_3& d) const
  {
    Vector_3 v = d.get_vector();
    Point_3  q = p + v;
	
    return Line_3(p,q);
  }
     
  Line_3 operator()(const Segment_3& s) const
  {
    return Line_3(s.source(), s.target());
  }  
     
  Line_3 operator()(const Ray_3& r) const
  {
    return Line_3(r.point1(), r.point2());
  }         
};

template<class K>
class Construct_leda_d3_rat_ray {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Direction_3 Direction_3;   
  typedef typename K::Vector_3  Vector_3;   
  typedef typename K::Ray_3     Ray_3;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Ray_3       result_type;

  Ray_3 operator()() const
  {
    Ray_3 r;
    return r;
  }

  Ray_3 operator()(const Point_3& p, const Point_3& q) const
  {
    return Ray_3(p,q);
  }
     
  Ray_3 operator()(const Point_3& p, const Direction_3& d) const
  {
    Vector_3 v = d.get_vector();
    Point_3  q = p + v;
	
    return Ray_3(p,q);
  }
};

#if defined(CGAL_COMPATIBLE_SPHERES)
template<class K>
class Construct_leda_d3_rat_sphere {

  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;  

public:
  typedef LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere       result_type;
  typedef Arity_tag< 4 >          Arity;
  
    LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const Point_3& center, 
                                                       const FT&  sqrad, 
						       CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center,sqrad,ori);
    } 
    
            
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const Point_3& p1, 
                                                        const Point_3& p2,
                                                        const Point_3& p3,
							const Point_3& p4) const
     {
      int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,p4);
      leda_d3_rat_sphere S(p1,p2,p3,p4);
      Point_3  center = S.center();
     
      CGAL::Orientation cg_ori;
      switch(ori){
        case -1: { cg_ori = CGAL::RIGHT_TURN; break; }
	case  0: { cg_ori = CGAL::COLLINEAR; break; }
	case  1: { cg_ori = CGAL::LEFT_TURN; break; }
      }
    
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center,p1,cg_ori);
     }

     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const Point_3& p1, 
                                                        const Point_3& p2,
                                                        const Point_3& p3, 
							CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
     {
       Point_3 center = leda_support::construct_circle_center_3(p1,p2,p3);
       FT     sq     = p1.sqr_dist(center);
       return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center, sq, ori);
     }   
     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const Point_3& p1, 
                                                        const Point_3& p2, 
                                                        CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
     {
      Point_3 m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(m,p1,ori);
     }       
     
     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const Point_3& p1, 
                                                        const CGAL::Orientation& ori = COUNTERCLOCKWISE) const
     {
        return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(p1,ori);
     }     
};
#else

template<class K>
class Construct_leda_d3_rat_sphere {

  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;  

public:
  typedef leda_d3_rat_sphere       result_type;
  typedef Arity_tag< 4 >           Arity;
  
  // squared radius version cannot be provided for LEDA rational kernel
     
  leda_d3_rat_sphere operator()(const Point_3& p1, const Point_3& p2,
                                const Point_3& p3, const Point_3& p4) const
  {
     return leda_d3_rat_sphere(p1,p2,p3,p4);
  }
     
     
  leda_d3_rat_sphere operator()(const Point_3& p1, const CGAL::Orientation& ori = COUNTERCLOCKWISE) const
  {
     return leda_d3_rat_sphere(p1,p1,p1,p1);
  }     
};
#endif

template<class K>
class Construct_leda_d3_rat_segment {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Segment_3 Segment_3;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Segment_3       result_type;

  Segment_3 operator()() const
  {
    Segment_3 s;
    return s;
  }

  Segment_3 operator()(const Point_3& p, const Point_3& q) const
  {
    return Segment_3(p,q);
  }
};

template<class K>
class Construct_leda_d3_rat_triangle {
  typedef typename K::Point_3     Point_3;
  typedef typename K::Triangle_3  Triangle_3;  

public:
  typedef Arity_tag< 3 > Arity;
  typedef Triangle_3       result_type;

  Triangle_3 operator()() const
  {
    Triangle_3 t;
    return t;
  }

  Triangle_3 operator()(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  {
    return Triangle_3(p1,p2,p3);
  }
};

template<class K>
class Construct_leda_d3_rat_simplex {

  typedef typename K::Point_3        Point_3;
  typedef typename K::Tetrahedron_3  Tetrahedron_3;  

public:
  typedef Arity_tag< 4 > Arity;
  typedef Tetrahedron_3       result_type;

  Tetrahedron_3 operator()() const
  {
    Tetrahedron_3 s;
    return s;
  }

  Tetrahedron_3 operator()(const Point_3& p1, const Point_3& p2,
                           const Point_3& p3, const Point_3& p4) const
  {
    return Tetrahedron_3(p1,p2,p3,p4);
  }
};

class Construct_leda_d3_rat_object {
public:
  typedef Arity_tag< 1 > Arity;
  typedef CGAL::Object       result_type;

  template<class T>
  CGAL::Object operator()(const T& obj) const
  { return CGAL::make_object(obj); }
};

template<class K>
class Construct_leda_d3_rat_scaled_vector {

  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_3       result_type;

  Vector_3 operator()(const Vector_3& v, const RT& scale) const
  {
   return scale * v;
  }      

  Vector_3 operator()(const Vector_3& v, const CGAL::Quotient<RT>& scale) const
  {
   FT fkt(scale.numerator(), scale.denominator());
   return fkt * v;
  } 
};

template<class K>
class Construct_leda_d3_rat_translated_point {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Vector_3  Vector_3; 

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_3       result_type;

  Point_3 operator()(const Point_3& p, const Vector_3& v) const
  {
   return p.translate(v);
  }         
};

template<class K>
class Construct_leda_d3_rat_point_on {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3;
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::Line_3    Line_3;
  typedef typename K::Ray_3     Ray_3;
  typedef typename K::Segment_3 Segment_3;   
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;  
  
public:
  typedef Point_3       result_type;
  typedef Arity_tag< 2 >          Arity;

  Point_3 operator()(const Line_3& l, int i = 0) const
  {
    Vector_3  v = l.to_vector();    
    return l.point1() + (RT(i) * v);     
  }
    
  Point_3 operator()(const Plane_3& pl) const
  { return pl.point1(); }
    
  Point_3 operator()(const Ray_3& r, int i = 0) const
  {
    if (i==0) return r.source(); 
    return r.point2(); // return a point different from the source ...
  }
    
  Point_3 operator()(const Segment_3& s, int i = 0) const
  {
    i = i % 2;
    if (i==0) return s.source();      
    return s.target();
  }               
};

template<class K>
class Construct_leda_d3_rat_projected_point {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3;
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::Line_3    Line_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_3       result_type;

  // orthogonal projection ...

  Point_3 operator()(const Line_3& l, const Point_3& p) const
  {
    // see CGAL implementation
    if (l.contains(p)) return p;

    // p not on l ...
    Point_3 pl = l.point1(); 
     
    Vector_3  v  = p  - pl;
    Vector_3  lv = l.to_vector();
    
    // homog. coordinates of vector from point on line to p ...
    const RT& vx = v.hcoord(0);
    const RT& vy = v.hcoord(1);
    const RT& vz = v.hcoord(2);
    const RT& vw = v.hcoord(3);
    
    // homog. coordinates of direction vector of the line
    const RT& lvx = lv.hcoord(0);
    const RT& lvy = lv.hcoord(1);
    const RT& lvz = lv.hcoord(2);
    const RT& lvw = lv.hcoord(3);

    RT numerator = (vx*lvx + vy*lvy + vz*lvz)*lvw; 
    RT denom = (lvx*lvx + lvy*lvy + lvz*lvz)*vw; 
    
    // add to pl ...
    return pl + ((numerator*lv)/denom);      
  }

    
  Point_3 operator()(const Plane_3& pl, const Point_3& p) const
  {
       Vector_3 v = pl.normal_project(p);
       return p + v;
  }    
};

template<class K>
class Construct_leda_d3_rat_lifted_point {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Point_2   Point_2; 
  typedef typename K::Vector_3  Vector_3;   
  typedef typename K::Plane_3   Plane_3;
  typedef typename K::FT        FT;
  typedef typename K::RT        RT;  

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_3       result_type;

    Point_3 operator()(const Plane_3& pl, const Point_2& p) const
    {
      // precondition ??? -> C must be nonzero !!!!
      // other variant: compute intersection with a line   
      
      //std::cout << pl << "\n";
#if (__LEDA__ >= 440)           
      FT A(pl.A());
      FT B(pl.B());
      FT C(pl.C());
      FT D(pl.D());
#else
      // no A(), B(), C(), D() for plane ...
      Point_3 pt = pl.point1();
      Vector_3 n = pl.normal();
      RT xi = n.X();
      RT yi = n.Y();
      RT zi = n.Z();
      
      FT A(xi*pt.W());
      FT B(yi*pt.W());
      FT C(zi*pt.W());
      FT D(-xi*pt.X()-yi*pt.Y()-zi*pt.Z());      
#endif      
      
      // compute z value
      FT z;
      
      if (C != 0){
        z = (- A*p.xcoord() - B*p.ycoord() - D)/C;  
      }
      else {
        z = 0;
	// should we call an error handler because the
	// lifting is not possible ???
      }
      
      return Point_3(p.xcoord(),p.ycoord(),z);    
    }
};

// this is a rather strange functor from the 2d kernel that deals 
// with projections of 3d points into the xy - plane ...

template<class K>
class Construct_leda_rat_projected_xy_point {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Point_2   Point_2; 
  typedef typename K::Plane_3   Plane_3;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_2       result_type;

  Point_2 operator()(const Plane_3& pl, const Point_3& p) const
  {
       // we need a projection that "maps" the plane into the
       // xy plane;
       // it seems to be forbidden to have a result of the plane
       // projection that is a line (?)
       
       // so first project the point1/2/3 of pl into the xy plane and see 
       // if it is a line !!!!
       
       Point_2 p1 = pl.point1().project_xy();
       Point_2 p2 = pl.point2().project_xy();
       Point_2 p3 = pl.point3().project_xy();
       
       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3) ){
          // projecting into xy - plane by removing z is OK !!!
          return Point_2(p.xcoord(),p.ycoord());
       }
       // projection doesn't map the 3d plane into the xy plane ...
       // try other projections ...

       p1 = pl.point1().project_yz();
       p2 = pl.point2().project_yz();
       p3 = pl.point3().project_yz();       

       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3) ){
          return Point_2(p.ycoord(),p.zcoord());
       }
       
       return Point_2(p.xcoord(),p.zcoord());
  }
};

template<class K>
class Construct_leda_d3_rat_vertex {

  typedef typename K::Point_3       Point_3;
  typedef typename K::Segment_3     Segment_3;   
  typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename K::Triangle_3    Triangle_3;
  typedef typename K::Tetrahedron_3 Tetrahedron_3;      

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_3       result_type;

  Point_3 operator()(const Segment_3& s, int i) const
  {
    i = i % 2;
      
    if (i==0) return s.source();
    return s.target();
  }
    
  Point_3 operator()(const Iso_cuboid_3& c, int i) const
  {
    i = i % 8;
    return c.vertex(i);
  }
    
  Point_3 operator()(const Triangle_3& t, int i) const
  {
    i = i % 3;
    switch(i){
       case  0: { return t.point1(); }
       case  1: { return t.point2(); }
       default: { return t.point3(); }
    }
  }
    
  Point_3 operator()(const Tetrahedron_3& s, int i) const
  {
    i = i % 4;
    switch(i){
       case  0: { return s.point1(); }
       case  1: { return s.point2(); }
       case  2: { return s.point3(); }
       default: { return s.point4(); }      
    }
  }
};

template<class K>
class Construct_leda_d3_rat_supporting_line {

  typedef typename K::Line_3    Line_3;
  typedef typename K::Ray_3     Ray_3;
  typedef typename K::Segment_3 Segment_3;   

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_3       result_type;

  Line_3 operator()(const Ray_3& r) const
  { return Line_3(r.point1(), r.point2()); }
    
  Line_3 operator()(const Segment_3& s) const
  { return Line_3(s); }     
};

template<class K>
class Construct_leda_d3_rat_supporting_plane {

  typedef typename K::Plane_3      Plane_3;  
  typedef typename K::Triangle_3   Triangle_3; 
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;    

public:
  typedef Arity_tag< 1 > Arity;
  typedef Plane_3       result_type;
   
  Plane_3 operator()(const Triangle_3& t) const
  {
      Point_3 p1 = t.point1();
      Point_3 p2 = t.point2();
      Point_3 p3 = t.point3();
      
      if (LEDA_NAMESPACE_NAME::collinear(p1,p2,p3)){ // special case ...
         if (p1==p2 && p2==p3) { // the triangle is a point ...
	   Vector_3 v1(1,0,0,1);
	   Vector_3 v2(0,1,0,1);
	   p2 = p1 + v1;
	   p3 = p1 + v2;	   
	 }
	 else { // triangle is degenerated to a segment ...
	   if (p1==p2) p2=p3;
           // how we have 2 points on the segment; get another point
	   // not on the segment ....
	   Point_3 other(0,0,0,1);
	   if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	   else {
	     other = Point_3(1,0,0,1);
	     if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	     else {
	       other = Point_3(0,1,0,1);
	       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	       else p3 = Point_3(0,0,1,1);
	     }
	   } 
	 }
      }     
      return Plane_3(p1,p2,p3);            
  }
};

template<class K>
class Construct_leda_d3_rat_orthogonal_vector {

  typedef typename K::Plane_3      Plane_3; 
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;    


public:
  typedef Arity_tag< 1 > Arity;
  typedef Vector_3       result_type;

  Vector_3 operator()(const Plane_3& pl) const
  {
       // returned vector must be oriented from the negative to the positive
       // side of pl ...
       
       Vector_3 n = pl.normal();
       
       Point_3 p1 = pl.point1();
       Point_3 p2 = p1 + n;
       
       // in what direction points the normal vector ???
       int ori = LEDA_NAMESPACE_NAME::orientation(pl, p2);
       
       if (ori == +1) return n;
       
       return -n;
    }
};

template<class K>
class Construct_leda_d3_rat_base_vector {

  typedef typename K::Point_3      Point_3;
  typedef typename K::Plane_3      Plane_3;  
  typedef typename K::Vector_3     Vector_3;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_3       result_type;

  // construct orthogonal vectors with some precondition for pl ...

  Vector_3 operator()(const Plane_3& pl, int index) const
  {
       Point_3 p1 = pl.point1();
       Point_3 p2 = pl.point2();
       Point_3 p3 = pl.point3();
       
       if (index==1){ // return vector b1 orthogonal to normal ...
         Vector_3 b1 = p3 - p2;
	 return b1;
       }
       else { // index == 2
          // return vector orthogonal to normal and b1; also some orientation prec.
	  // must be fulfilled ...
	  Vector_3 n = pl.normal();
	  Plane_3 other(p2,p3,p2 + n);
	  
	  // get the normal projection vector for p1 and other plane ...
	  Vector_3 b2 = other.normal_project(p1);
	  
	  // this vector is orthogonal, but is the orientation correct ???
	  int res = LEDA_NAMESPACE_NAME::orientation(p2,p3,p2+b2,p2+n);
	  
	  if (res == 1) return b2;
	  
	  return -b2;
       }
  }
};

template<class K>
class Construct_leda_d3_rat_perpendicular_plane {

  typedef typename K::Plane_3      Plane_3;  
  typedef typename K::Point_3      Point_3;
  typedef typename K::Line_3       Line_3;    

public:
  typedef Arity_tag< 2 > Arity;
  typedef Plane_3       result_type;

  Plane_3 operator()(const Line_3& l, const Point_3& p) const
  { return Plane_3(p, l.to_vector()); }
};

template<class K>
class Construct_leda_d3_rat_perpendicular_line {

  typedef typename K::Plane_3      Plane_3;  
  typedef typename K::Point_3      Point_3;
  typedef typename K::Line_3       Line_3;
  typedef typename K::Vector_3     Vector_3;      

public:
  typedef Arity_tag< 2 > Arity;
  typedef Line_3       result_type;

  Line_3 operator()(const Plane_3& pl, const Point_3& p1) const
  {
       // returned line must be oriented from the negative to the positive
       // side of pl ...
       
       Vector_3 n = pl.normal();
       
       Point_3 p2 = p1 + n;
       
       // in what direction points the normal vector ???
       int ori = LEDA_NAMESPACE_NAME::orientation(pl, pl.point1() + n);
       
       if (ori == +1) return Line_3(p1,p2);
       
       return Line_3(p2,p1);
  }
};

template<class K>
class Construct_leda_d3_rat_midpoint {

  typedef typename K::Point_3      Point_3;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Point_3       result_type;

  Point_3 operator()(const Point_3& p1, const Point_3& p2) const
  { return LEDA_NAMESPACE_NAME::midpoint(p1,p2); }    
};

template<class K>
class Construct_leda_d3_rat_center {

  typedef typename K::Point_3      Point_3;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Point_3       result_type;

#if defined(CGAL_COMPATIBLE_SPHERES)
  Point_3 operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& S) const
  {
       return S.center();
  }
#else  
  Point_3 operator()(const leda_d3_rat_sphere& S) const
  {
       return S.center();
  }
#endif  
};

template<class K>
class Construct_leda_d3_rat_centroid {

  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT;  

public:
  typedef Point_3          result_type;
  typedef Arity_tag< 4 >   Arity;

  Point_3 operator()(const Point_3& p1, const Point_3& p2, const Point_3& p3) const
  {
       // sum up coordinates, divide by 3
       FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() )/ FT(3.0);
       FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() )/ FT(3.0);
       FT z = (p1.zcoord() + p2.zcoord() + p3.zcoord() )/ FT(3.0);
       
       return Point_3(x,y,z);
  }

  Point_3 operator()(const Point_3& p1, const Point_3& p2, 
                     const Point_3& p3, const Point_3& p4) const
  {
       // sum up coordinates, divide by 4
       FT x = (p1.xcoord() + p2.xcoord() + p3.xcoord() + p4.xcoord() )/ FT(4.0);
       FT y = (p1.ycoord() + p2.ycoord() + p3.ycoord() + p4.ycoord() )/ FT(4.0);
       FT z = (p1.zcoord() + p2.zcoord() + p3.zcoord() + p4.zcoord() )/ FT(4.0);
       
       return Point_3(x,y,z);       
  }
};    

template<class K>
class Construct_leda_d3_rat_circumcenter {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3;  
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::FT        FT;  

public:
  typedef Point_3          result_type;
  typedef Arity_tag< 4 >   Arity;

   Point_3 operator()(const Point_3& p1, const Point_3& p2,
                      const Point_3& p3, const Point_3& p4) const
   {
     leda_d3_rat_sphere S(p1,p2,p3,p4);     
     return S.center();
   }  
  
   //variant circle on a plane:
   //taken from GEOMLEP
   //can we improve this ?????
   
   FT r_det(FT a11,FT a12,FT a13,
            FT a21,FT a22,FT a23,
            FT a31,FT a32,FT a33) const
   {
     return a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a11*a23*a32-a12*a21*a33;
   }    
   
   
   Point_3 intersec(Vector_3 n1,Vector_3 n2,Vector_3 n3,
                    FT D1,FT D2,FT D3) const
   {
     FT xw,yw,zw;
     FT det,xdet,ydet,zdet;
     FT A1,B1,C1,A2,B2,C2,A3,B3,C3;

     A1=n1.xcoord(); B1=n1.ycoord(); C1=n1.zcoord();
     A2=n2.xcoord(); B2=n2.ycoord(); C2=n2.zcoord();
     A3=n3.xcoord(); B3=n3.ycoord(); C3=n3.zcoord();

     det= r_det(A1,B1,C1,A2,B2,C2,A3,B3,C3);
     xdet=r_det(D1,B1,C1,D2,B2,C2,D3,B3,C3);
     ydet=r_det(A1,D1,C1,A2,D2,C2,A3,D3,C3);
     zdet=r_det(A1,B1,D1,A2,B2,D2,A3,B3,D3);

     /*if (det==0) error_handler(1,"intersec: determinant==0 .");*/

     xw=-xdet/det; yw=-ydet/det; zw=-zdet/det;
     xw=xw.normalize();
     yw=yw.normalize();
     zw=zw.normalize();

     return Point_3(xw,yw,zw);
   }
   
   
   Point_3 operator()(const Point_3& a, const Point_3& b, const Point_3& c) const
   {
     Point_3 pneu;
     Point_3 m1,m2;
     Plane_3 pl1(a,b,c);
     Point_3 h=c + pl1.normal();
     Plane_3 plh(a,b,h); 
     Plane_3 plh2(b,c,h);
  
     m1= LEDA_NAMESPACE_NAME::midpoint(a,b);
     m2= LEDA_NAMESPACE_NAME::midpoint(b,c);
 
     Vector_3 r1=pl1.normal();
     Vector_3 r2=plh.normal();
     Vector_3 r3=plh2.normal();
     Plane_3 pl2(m1,m1+r1,m1+r2),pl3(m2,m2+r1,m2+r3);

     Vector_3 n1,n2,n3;
     n1=pl1.normal(); n2=pl2.normal(); n3=pl3.normal();

     FT D1,D2,D3;
     D1=-(n1*m1.to_vector()); D2=-(n2*m1.to_vector()); D3=-(n3*m2.to_vector());

     //compute intersection of the three planes
     pneu=intersec(n1,n2,n3,D1,D2,D3);
     return pneu;        
   }    
};

template<class K>
class Construct_leda_d3_rat_cross_product_vector {

  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::FT        FT;

public:
  typedef Arity_tag< 2 > Arity;
  typedef Vector_3       result_type;
 
  Vector_3 operator()(const Vector_3& a, const Vector_3& b) const
  {
    CGAL_precondition(a.dim() == 3);
    CGAL_precondition(b.dim() == 3);    
    
    FT a1 = a.xcoord();
    FT a2 = a.ycoord();
    FT a3 = a.zcoord();
    FT b1 = b.xcoord();
    FT b2 = b.ycoord();
    FT b3 = b.zcoord();    
    
    FT v1 = a2*b3-a3*b2;
    FT v2 = a3*b1-a1*b3;
    FT v3 = a1*b2-a2*b1;
    
    return Vector_3(v1,v2,v3);
  }
};

template<class K>
class Construct_leda_d3_rat_opposite_direction {

  typedef typename K::Vector_3    Vector_3;
  typedef typename K::Direction_3 Direction_3;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Direction_3       result_type;

  Direction_3 operator()(const Direction_3& dir) const
  { Vector_3 v = dir.get_vector(); 
    return Direction_3(-v);
  }
};

template<class K>
class Construct_leda_d3_rat_opposite_segment {

  typedef typename K::Segment_3 Segment_3;   

public:
  typedef Arity_tag< 1 > Arity;
  typedef Segment_3       result_type;

  Segment_3 operator()(const Segment_3& s) const
  { return s.reverse(); }
};

template<class K>
class Construct_leda_d3_rat_opposite_ray {

  typedef typename K::Ray_3     Ray_3;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Ray_3       result_type;

  Ray_3 operator()(const Ray_3& r) const
  { return r.reverse(); }
};

template<class K>
class Construct_leda_d3_rat_opposite_line {

  typedef typename K::Line_3    Line_3;

public:
  typedef Arity_tag< 1 > Arity;
  typedef Line_3       result_type;

  Line_3 operator()(const Line_3& l) const
  { return l.reverse(); }
};

template<class K>
class Construct_leda_d3_rat_opposite_plane {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3; 

public:
  typedef Arity_tag< 1 > Arity;
  typedef Plane_3       result_type;

  Plane_3 operator()(const Plane_3& pl) const
  {
    Point_3 p1 = pl.point1();
    Point_3 p2 = pl.point2();
    Point_3 p3 = pl.point3();

    return Plane_3(p3,p2,p1);
  }
};

template<class K>
class Construct_leda_d3_rat_opposite_sphere {

  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT; 

public:
  typedef Arity_tag< 1 > Arity;
  
#if defined(CGAL_COMPATIBLE_SPHERES)  
  typedef LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere       result_type;

  LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s) const
  {
    Point_3 p  = s.center();
    FT  sq = s.sqr_radius();
    CGAL::Orientation ori = s.orientation();
    
    ori = reverse_orientation(ori);
  
    return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(p,sq,ori);   
  }
#else
  typedef leda_d3_rat_sphere       result_type;

  leda_d3_rat_sphere operator()(const leda_d3_rat_sphere& s) const
  {
    Point_3 p1 = s.point1();
    Point_3 p2 = s.point2();
    Point_3 p3 = s.point3();
    Point_3 p4 = s.point4();   
    
    return leda_d3_rat_sphere(p4,p3,p2,p1);
  }
#endif   
};

template<class K>
class Construct_leda_d3_rat_opposite_vector {

  typedef typename K::Vector_3  Vector_3; 

public:
  typedef Arity_tag< 1 > Arity;
  typedef Vector_3       result_type;
  
  Vector_3 operator()(const Vector_3& v) const
  { return -v; }
};


CGAL_END_NAMESPACE

#endif




