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

CGAL_BEGIN_NAMESPACE

class Construct_leda_d3_rat_point {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_point           result_type;
  
  //----------------------------------------------------------------------------------
  //undocumented
  leda_d3_rat_point operator()() const
  { 
   leda_d3_rat_point p;
   return p;
  } 
  
  leda_d3_rat_point operator()(const leda_rational& x, 
                               const leda_rational& y, const leda_rational& z) const
  {
    leda_d3_rat_point p(x,y,x);
    return p;
  }  
  
  leda_d3_rat_point operator()(const leda_integer& x, const leda_integer& y, 
                            const leda_integer& z, const leda_integer& w) const
  {
    leda_d3_rat_point p(x,y,z,w);
    return p;
  } 
  
  //----------------------------------------------------------------------------------

  leda_d3_rat_point operator()(const CGAL::Origin& orig) const
  { return leda_d3_rat_point(0,0,0,1); }
};

class Construct_leda_d3_rat_vector {
public:
  typedef leda_rat_vector           result_type;
  typedef Arity_tag< 2 >            Arity;

  //----------------------------------------------------------------------------------   
  //undocumented
  leda_rat_vector operator()() const
  { 
   leda_rat_vector v(2);
   return v;
  }  
  
  leda_rat_vector operator()(const leda_rational& x, const leda_rational& y, 
                             const leda_rational& z) const
  {
    leda_rat_vector v(x,y,z);
    return v;
  }  
  
  leda_rat_vector operator()(const leda_integer& x, const leda_integer& y, 
                             const leda_integer& z, const leda_integer& w) const
  {
    leda_rat_vector v(x,y,z,w);
    return v;
  } 
  
  //---------------------------------------------------------------------------------- 

   leda_rat_vector operator()(const leda_d3_rat_point& a, const leda_d3_rat_point& b) const
   { return b-a; }
   
   leda_rat_vector operator()(const CGAL::Null_vector& nv) const
   { return leda_rat_vector(3); }
};

// Achtung : fuehre hier typedef aus fuer Construct_direction_of_line_3
// (in der Traits - Klasse)
// welches in 3d Delaunay gebraucht wird ...

class Construct_leda_d3_rat_direction {
public:
  typedef Arity_tag< 1 > Arity;
  typedef LEDA_NAMESPACE_NAME::rat_direction           result_type;
  
  //undocumented
  LEDA_NAMESPACE_NAME::rat_direction operator()() const
  { 
   LEDA_NAMESPACE_NAME::rat_direction d(3);
   return d;
  } 
  
  LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rational& x, const leda_rational& y,
                                                const leda_rational& z) const
  {
    LEDA_NAMESPACE_NAME::rat_direction d(x,y,z);
    return d;
  }          
  
    LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_rat_vector& v) const
    { return LEDA_NAMESPACE_NAME::rat_direction(v); }
    
    LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_d3_rat_line& l) const
    { leda_rat_vector v = l.point2()-l.point1();
      return LEDA_NAMESPACE_NAME::rat_direction(v); 
    }
    
    LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_d3_rat_ray& r) const
    {  leda_rat_vector v = r.point2()-r.point1(); 
       return LEDA_NAMESPACE_NAME::rat_direction(v);
    }
    
    LEDA_NAMESPACE_NAME::rat_direction operator()(const leda_d3_rat_segment& s) const
    { leda_rat_vector v = s.target()-s.source();
      return LEDA_NAMESPACE_NAME::rat_direction(v); 
    }            
};

class Construct_leda_d3_rat_plane {
public:
  typedef leda_d3_rat_plane           result_type;
  typedef Arity_tag< 3 >              Arity;
  
  leda_d3_rat_plane operator()() const
  {
    leda_d3_rat_plane p;
    return p;
  }  

  leda_d3_rat_plane operator()(const leda_integer& a, const leda_integer& b,
                               const leda_integer& c, const leda_integer& d) const
  {
      CGAL_precondition(! (a==0 && b==0 && c==0));

      leda_integer x=0, y=0, z=0;
      leda_rational other;
      leda_integer n, den;
      leda_d3_rat_point on_plane; // a point on the plane ...

      // compute point (x,y,z) on plane with
      // a*x + b*y + c*z + d = 0
      
      if (c!=0){
        other = leda_rational(-a*x - b*y - d,c);
	n = other.numerator();
	den = other.denominator();
	on_plane = leda_d3_rat_point(x*den,y*den,n,den);
      }
      else {
        if (b!=0){
	  other = leda_rational(-a*x - c*z -d,b);
	  n = other.numerator();
	  den = other.denominator();
	  on_plane = leda_d3_rat_point(x*den,n,z*den,den);	  
	}
	else {
	  other = leda_rational(-b*y - c*z -d,a);
	  n = other.numerator();
	  den = other.denominator();
	  on_plane = leda_d3_rat_point(n,y*den,z*den,den);	  
	}
      }
      
      // normal vector ...
      leda_rat_vector nv(a,b,c,1);
      
      return leda_d3_rat_plane(on_plane, nv);
  }
    
  leda_d3_rat_plane operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                               const leda_d3_rat_point& p3) const
  {
       return leda_d3_rat_plane(p1,p2,p3);
  }   
    
  leda_d3_rat_plane operator()(const leda_d3_rat_point& p, const LEDA_NAMESPACE_NAME::rat_direction& d) const
  {
       // construction with normal vector ...
       return leda_d3_rat_plane(p, d.get_vector());
  }
    
  leda_d3_rat_plane operator()(const leda_d3_rat_line& l, const leda_d3_rat_point& p) const
  {
       return leda_d3_rat_plane(l.point1(), l.point2(), p);
  }
    
  leda_d3_rat_plane operator()(const leda_d3_rat_ray& r, const leda_d3_rat_point& p) const
  {
       return leda_d3_rat_plane(r.point1(), r.point2(), p);    
  }
    
  leda_d3_rat_plane operator()(const leda_d3_rat_segment& s, const leda_d3_rat_point& p) const
  {
       return leda_d3_rat_plane(s.source(), s.target(), p);        
  }            
    
};

class Construct_leda_d3_rat_iso_cuboid {
public:
  typedef Arity_tag< 2 > Arity;
  typedef LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid           result_type;

  LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid operator()() const
  {
    LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid i;
    return i;
  }

  LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
  {
     return LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid(p,q);
  }
};

class Construct_leda_d3_rat_line {
public:
  typedef leda_d3_rat_line           result_type;  
  typedef Arity_tag< 2 >             Arity;

  leda_d3_rat_line operator()() const
  {
    leda_d3_rat_line l;
    return l;
  }

  leda_d3_rat_line operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
  {
    return leda_d3_rat_line(p,q);
  }
     
  leda_d3_rat_line operator()(const leda_d3_rat_point& p, const LEDA_NAMESPACE_NAME::rat_direction& d) const
  {
    leda_rat_vector v = d.get_vector();
    leda_d3_rat_point  q = p + v;
	
    return leda_d3_rat_line(p,q);
  }
     
  leda_d3_rat_line operator()(const leda_d3_rat_segment& s) const
  {
    return leda_d3_rat_line(s.source(), s.target());
  }  
     
  leda_d3_rat_line operator()(const leda_d3_rat_ray& r) const
  {
    return leda_d3_rat_line(r.point1(), r.point2());
  }         
};

class Construct_leda_d3_rat_ray {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_ray       result_type;

  leda_d3_rat_ray operator()() const
  {
    leda_d3_rat_ray r;
    return r;
  }

  leda_d3_rat_ray operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
  {
    return leda_d3_rat_ray(p,q);
  }
     
  leda_d3_rat_ray operator()(const leda_d3_rat_point& p, const LEDA_NAMESPACE_NAME::rat_direction& d) const
  {
    leda_rat_vector v = d.get_vector();
    leda_d3_rat_point  q = p + v;
	
    return leda_d3_rat_ray(p,q);
  }
};

#if defined(CGAL_COMPATIBLE_SPHERES)
class Construct_leda_d3_rat_sphere {
public:
  typedef LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere       result_type;
  typedef Arity_tag< 4 >          Arity;
  
    LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const leda_d3_rat_point& center, 
                                                       const leda_rational&  sqrad, 
						       CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
    {
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center,sqrad,ori);
    } 
    
            
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const leda_d3_rat_point& p1, 
                                                        const leda_d3_rat_point& p2,
                                                        const leda_d3_rat_point& p3,
							const leda_d3_rat_point& p4) const
     {
      int ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,p4);
      leda_d3_rat_sphere S(p1,p2,p3,p4);
      leda_d3_rat_point  center = S.center();
     
      CGAL::Orientation cg_ori;
      switch(ori){
        case -1: { cg_ori = CGAL::RIGHTTURN; break; }
	case  0: { cg_ori = CGAL::COLLINEAR; break; }
	case  1: { cg_ori = CGAL::LEFTTURN; break; }
      }
    
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center,p1,cg_ori);
     }

     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const leda_d3_rat_point& p1, 
                                                        const leda_d3_rat_point& p2,
                                                        const leda_d3_rat_point& p3, 
							CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
     {
       leda_d3_rat_point center = leda_support::construct_circle_center_3(p1,p2,p3);
       leda_rational     sq     = p1.sqr_dist(center);
       return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center, sq, ori);
     }   
     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const leda_d3_rat_point& p1, 
                                                        const leda_d3_rat_point& p2, 
                                                        CGAL::Orientation ori = CGAL::COUNTERCLOCKWISE) const
     {
      leda_d3_rat_point m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     
      return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(m,p1,ori);
     }       
     
     
     LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const leda_d3_rat_point& p1, 
                                                        const CGAL::Orientation& ori = COUNTERCLOCKWISE) const
     {
        return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(p1,ori);
     }     
};
#else
class Construct_leda_d3_rat_sphere {
public:
  typedef leda_d3_rat_sphere       result_type;
  typedef Arity_tag< 4 >           Arity;
  
     // squared radius version cannot be provided for LEDA rational kernel
     
     leda_d3_rat_sphere operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                                   const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
     {
        return leda_d3_rat_sphere(p1,p2,p3,p4);
     }
     
     
     leda_d3_rat_sphere operator()(const leda_d3_rat_point& p1, const CGAL::Orientation& ori = COUNTERCLOCKWISE) const
     {
        return leda_d3_rat_sphere(p1,p1,p1,p1);
     }     
};
#endif

class Construct_leda_d3_rat_segment {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_segment       result_type;

  leda_d3_rat_segment operator()() const
  {
    leda_d3_rat_segment s;
    return s;
  }

  leda_d3_rat_segment operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
  {
    return leda_d3_rat_segment(p,q);
  }
};

class Construct_leda_d3_rat_triangle {
public:
  typedef Arity_tag< 3 > Arity;
  typedef leda_d3_rat_triangle       result_type;

  leda_d3_rat_triangle operator()() const
  {
    leda_d3_rat_triangle t;
    return t;
  }

  leda_d3_rat_triangle operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, const leda_d3_rat_point& p3) const
  {
    return leda_d3_rat_triangle(p1,p2,p3);
  }
};

class Construct_leda_d3_rat_simplex {
public:
  typedef Arity_tag< 4 > Arity;
  typedef LEDA_NAMESPACE_NAME::d3_rat_simplex       result_type;

  LEDA_NAMESPACE_NAME::d3_rat_simplex operator()() const
  {
    LEDA_NAMESPACE_NAME::d3_rat_simplex s;
    return s;
  }

  LEDA_NAMESPACE_NAME::d3_rat_simplex operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                                                 const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
  {
    return LEDA_NAMESPACE_NAME::d3_rat_simplex(p1,p2,p3,p4);
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

class Construct_leda_d3_rat_scaled_vector {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_vector       result_type;

    leda_rat_vector operator()(const leda_rat_vector& v, const leda_integer& scale) const
    {
     return scale * v;
    }      

    leda_rat_vector operator()(const leda_rat_vector& v, const CGAL::Quotient<leda_integer>& scale) const
    {
     leda_rational fkt(scale.numerator(), scale.denominator());
     return fkt * v;
    } 
};

class Construct_leda_d3_rat_translated_point {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_point       result_type;

    leda_d3_rat_point operator()(const leda_d3_rat_point& p, const leda_rat_vector& v) const
    {
     return p.translate(v);
    }         
};

class Construct_leda_d3_rat_point_on {
public:
  typedef leda_d3_rat_point       result_type;
  typedef Arity_tag< 2 >          Arity;

    leda_d3_rat_point operator()(const leda_d3_rat_line& l, int i = 0) const
    {
      leda_rat_vector  v = l.to_vector();
      
      return l.point1() + (leda_integer(i) * v);     
    }
    
    leda_d3_rat_point operator()(const leda_d3_rat_plane& pl) const
    {
      return pl.point1();
    }
    
    leda_d3_rat_point operator()(const leda_d3_rat_ray& r, int i = 0) const
    {
      if (i==0) return r.source(); 
      return r.point2(); // return a point different from the source ...
    }
    
    leda_d3_rat_point operator()(const leda_d3_rat_segment& s, int i = 0) const
    {
      i = i % 2;
      if (i==0) return s.source();
      
      return s.target();
    }               
};

class Construct_leda_d3_rat_projected_point {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_point       result_type;

  // orthogonal projection ...

  leda_d3_rat_point operator()(const leda_d3_rat_line& l, const leda_d3_rat_point& p) const
  {
    // see CGAL implementation
    if (l.contains(p)) return p;

    // p not on l ...
    leda_d3_rat_point pl = l.point1(); 
     
    leda_rat_vector  v  = p  - pl;
    leda_rat_vector  lv = l.to_vector();
    
    // homog. coordinates of vector from point on line to p ...
    const leda_integer& vx = v.hcoord(0);
    const leda_integer& vy = v.hcoord(1);
    const leda_integer& vz = v.hcoord(2);
    const leda_integer& vw = v.hcoord(3);
    
    // homog. coordinates of direction vector of the line
    const leda_integer& lvx = lv.hcoord(0);
    const leda_integer& lvy = lv.hcoord(1);
    const leda_integer& lvz = lv.hcoord(2);
    const leda_integer& lvw = lv.hcoord(3);

    leda_integer numerator = (vx*lvx + vy*lvy + vz*lvz)*lvw; 
    leda_integer denom = (lvx*lvx + lvy*lvy + lvz*lvz)*vw; 
    
    // add to pl ...
    return pl + ((numerator*lv)/denom);      
  }

    
  leda_d3_rat_point operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
  {
       leda_rat_vector v = pl.normal_project(p);
       return p + v;
  }    
};

class Construct_leda_d3_rat_lifted_point {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_point       result_type;

    leda_d3_rat_point operator()(const leda_d3_rat_plane& pl, const leda_rat_point& p) const
    {
      // precondition ??? -> C must be nonzero !!!!
      // other variant: compute intersection with a line   
      
      //std::cout << pl << "\n";
#if (__LEDA__ >= 440)           
      leda_rational A(pl.A());
      leda_rational B(pl.B());
      leda_rational C(pl.C());
      leda_rational D(pl.D());
#else
      // no A(), B(), C(), D() for plane ...
      leda_d3_rat_point pt = pl.point1();
      leda_rat_vector n = pl.normal();
      leda_integer xi = n.X();
      leda_integer yi = n.Y();
      leda_integer zi = n.Z();
      
      leda_rational A(xi*pt.W());
      leda_rational B(yi*pt.W());
      leda_rational C(zi*pt.W());
      leda_rational D(-xi*pt.X()-yi*pt.Y()-zi*pt.Z());      
#endif      
      
      // compute z value
      leda_rational z;
      
      if (C != 0){
        z = (- A*p.xcoord() - B*p.ycoord() - D)/C;  
      }
      else {
        z = 0;
	// should we call an error handler because the
	// lifting is not possible ???
      }
      
      return leda_d3_rat_point(p.xcoord(),p.ycoord(),z);    
    }
};

// this is a rather strange predicate from the 2d kernel that deals 
// with projections of 3d points into the xy - plane ...

class Construct_leda_rat_projected_xy_point {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_point       result_type;

    leda_rat_point operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
    {
       // we need a projection that "maps" the plane into the
       // xy plane;
       // it seems to be forbidden to have a result of the plane
       // projection that is a line (?)
       
       // so first project the point1/2/3 of pl into the xy plane and see 
       // if it is a line !!!!
       
       leda_rat_point p1 = pl.point1().project_xy();
       leda_rat_point p2 = pl.point2().project_xy();
       leda_rat_point p3 = pl.point3().project_xy();
       
       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3) ){
          // projecting into xy - plane by removing z is OK !!!
          return leda_rat_point(p.xcoord(),p.ycoord());
       }
       // projection doesn't map the 3d plane into the xy plane ...
       // try other projections ...

       p1 = pl.point1().project_yz();
       p2 = pl.point2().project_yz();
       p3 = pl.point3().project_yz();       

       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3) ){
          return leda_rat_point(p.ycoord(),p.zcoord());
       }
       
       return leda_rat_point(p.xcoord(),p.zcoord());
    }
};


class Construct_leda_d3_rat_vertex {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_point       result_type;

    leda_d3_rat_point operator()(const leda_d3_rat_segment& s, int i) const
    {
      i = i % 2;
      
      if (i==0) return s.source();
      return s.target();
    }
    
    leda_d3_rat_point operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& c, int i) const
    {
      i = i % 8;
      return c.vertex(i);
    }
    
    leda_d3_rat_point operator()(const leda_d3_rat_triangle& t, int i) const
    {
      i = i % 3;
      switch(i){
       case  0: { return t.point1(); }
       case  1: { return t.point2(); }
       default: { return t.point3(); }
      }
    }
    
    leda_d3_rat_point operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, int i) const
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

class Construct_leda_d3_rat_supporting_line {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_line       result_type;

    leda_d3_rat_line operator()(const leda_d3_rat_ray& r) const
    {
      return leda_d3_rat_line(r.point1(), r.point2());
    }
    
    leda_d3_rat_line operator()(const leda_d3_rat_segment& s) const
    {
      return leda_d3_rat_line(s);    
    }     
};

class Construct_leda_d3_rat_supporting_plane {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_plane       result_type;

    leda_d3_rat_plane operator()(const leda_d3_rat_triangle& t) const
    {
      leda_d3_rat_point p1 = t.point1();
      leda_d3_rat_point p2 = t.point2();
      leda_d3_rat_point p3 = t.point3();
      
      if (LEDA_NAMESPACE_NAME::collinear(p1,p2,p3)){ // special case ...
         if (p1==p2 && p2==p3) { // the triangle is a point ...
	   leda_rat_vector v1(1,0,0,1);
	   leda_rat_vector v2(0,1,0,1);
	   p2 = p1 + v1;
	   p3 = p1 + v2;	   
	 }
	 else { // triangle is degenerated to a segment ...
	   if (p1==p2) p2=p3;
           // how we have 2 points on the segment; get another point
	   // not on the segment ....
	   leda_d3_rat_point other(0,0,0,1);
	   if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	   else {
	     other = leda_d3_rat_point(1,0,0,1);
	     if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	     else {
	       other = leda_d3_rat_point(0,1,0,1);
	       if (! LEDA_NAMESPACE_NAME::collinear(p1,p2,other)) p3 = other;
	       else p3 = leda_d3_rat_point(0,0,1,1);
	     }
	   } 
	 }
      }
      
      return leda_d3_rat_plane(p1,p2,p3);            
    }
};

class Construct_leda_d3_rat_orthogonal_vector {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_vector       result_type;

    leda_rat_vector operator()(const leda_d3_rat_plane& pl) const
    {
       // returned vector must be oriented from the negative to the positive
       // side of pl ...
       
       leda_rat_vector n = pl.normal();
       
       leda_d3_rat_point p1 = pl.point1();
       leda_d3_rat_point p2 = p1 + n;
       
       // in what direction points the normal vector ???
       int ori = LEDA_NAMESPACE_NAME::orientation(pl, p2);
       
       if (ori == +1) return n;
       
       return -n;
    }
};

class Construct_leda_d3_rat_base_vector {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_vector       result_type;

    // construct orthogonal vectors with some precondition for pl ...

    leda_rat_vector operator()(const leda_d3_rat_plane& pl, int index) const
    {
       leda_d3_rat_point p1 = pl.point1();
       leda_d3_rat_point p2 = pl.point2();
       leda_d3_rat_point p3 = pl.point3();
       
       if (index==1){ // return vector b1 orthogonal to normal ...
         leda_rat_vector b1 = p3 - p2;
	 return b1;
       }
       else { // index == 2
          // return vector orthogonal to normal and b1; also some orientation prec.
	  // must be fulfilled ...
	  leda_rat_vector n = pl.normal();
	  leda_d3_rat_plane other(p2,p3,p2 + n);
	  
	  // get the normal projection vector for p1 and other plane ...
	  leda_rat_vector b2 = other.normal_project(p1);
	  
	  // this vector is orthogonal, but is the orientation correct ???
	  int res = LEDA_NAMESPACE_NAME::orientation(p2,p3,p2+b2,p2+n);
	  
	  if (res == 1) return b2;
	  
	  return -b2;
       }
    }
};

class Construct_leda_d3_rat_perpendicular_plane {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_plane       result_type;

    leda_d3_rat_plane operator()(const leda_d3_rat_line& l, const leda_d3_rat_point& p) const
    {
        return leda_d3_rat_plane(p, l.to_vector()); 
    }
};


class Construct_leda_d3_rat_perpendicular_line {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_line       result_type;

    leda_d3_rat_line operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p1) const
    {
       // returned line must be oriented from the negative to the positive
       // side of pl ...
       
       leda_rat_vector n = pl.normal();
       
       leda_d3_rat_point p2 = p1 + n;
       
       // in what direction points the normal vector ???
       int ori = LEDA_NAMESPACE_NAME::orientation(pl, pl.point1() + n);
       
       if (ori == +1) return leda_d3_rat_line(p1,p2);
       
       return leda_d3_rat_line(p2,p1);
    }
};

class Construct_leda_d3_rat_midpoint {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_d3_rat_point       result_type;

    leda_d3_rat_point operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
    {
       return LEDA_NAMESPACE_NAME::midpoint(p1,p2);
    }    
};

class Construct_leda_d3_rat_center {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_point       result_type;

#if defined(CGAL_COMPATIBLE_SPHERES)
  leda_d3_rat_point operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& S) const
  {
       return S.center();
  }
#else  
  leda_d3_rat_point operator()(const leda_d3_rat_sphere& S) const
  {
       return S.center();
  }
#endif  
};

class Construct_leda_d3_rat_centroid {
public:
  typedef leda_d3_rat_point       result_type;
  typedef Arity_tag< 4 >          Arity;

    leda_d3_rat_point operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, const leda_d3_rat_point& p3) const
    {
       // sum up coordinates, divide by 3
       leda_rational x = (p1.xcoord() + p2.xcoord() + p3.xcoord() )/ leda_rational(3.0);
       leda_rational y = (p1.ycoord() + p2.ycoord() + p3.ycoord() )/ leda_rational(3.0);
       leda_rational z = (p1.zcoord() + p2.zcoord() + p3.zcoord() )/ leda_rational(3.0);
       
       return leda_d3_rat_point(x,y,z);
    }

    leda_d3_rat_point operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                                 const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
    {
       // sum up coordinates, divide by 4
       leda_rational x = (p1.xcoord() + p2.xcoord() + p3.xcoord() + p4.xcoord() )/ leda_rational(4.0);
       leda_rational y = (p1.ycoord() + p2.ycoord() + p3.ycoord() + p4.ycoord() )/ leda_rational(4.0);
       leda_rational z = (p1.zcoord() + p2.zcoord() + p3.zcoord() + p4.zcoord() )/ leda_rational(4.0);
       
       return leda_d3_rat_point(x,y,z);       
    }
};    

class Construct_leda_d3_rat_circumcenter {
public:
  typedef leda_d3_rat_point       result_type;
  typedef Arity_tag< 4 >          Arity;

   leda_d3_rat_point operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                                const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
   {
     leda_d3_rat_sphere S(p1,p2,p3,p4);
     
     return S.center();
   }  
  
   //variant circle on a plane:
   // taken from GEOMLEP
   // can we improve this ?????
   
   leda_rational r_det(leda_rational a11,leda_rational a12,leda_rational a13,
                       leda_rational a21,leda_rational a22,leda_rational a23,
                       leda_rational a31,leda_rational a32,leda_rational a33) const
   {
     return a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a11*a23*a32-a12*a21*a33;
   }    
   
   
   leda_d3_rat_point intersec(leda_rat_vector n1,leda_rat_vector n2,leda_rat_vector n3,
                         leda_rational D1,leda_rational D2,leda_rational D3) const
   {
     leda_rational xw,yw,zw;
     leda_rational det,xdet,ydet,zdet;
     leda_rational A1,B1,C1,A2,B2,C2,A3,B3,C3;

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

     return leda_d3_rat_point(xw,yw,zw);
   }
   
   
   leda_d3_rat_point operator()(const leda_d3_rat_point& a, const leda_d3_rat_point& b,
                                const leda_d3_rat_point& c) const
   {
     leda_d3_rat_point pneu;
     leda_d3_rat_point m1,m2;
     leda_d3_rat_plane pl1(a,b,c);
     leda_d3_rat_point h=c + pl1.normal();
     leda_d3_rat_plane plh(a,b,h); 
     leda_d3_rat_plane plh2(b,c,h);
  
     m1= LEDA_NAMESPACE_NAME::midpoint(a,b);
     m2= LEDA_NAMESPACE_NAME::midpoint(b,c);
 
     leda_rat_vector r1=pl1.normal();
     leda_rat_vector r2=plh.normal();
     leda_rat_vector r3=plh2.normal();
     leda_d3_rat_plane pl2(m1,m1+r1,m1+r2),pl3(m2,m2+r1,m2+r3);

     leda_rat_vector n1,n2,n3;
     n1=pl1.normal(); n2=pl2.normal(); n3=pl3.normal();

     leda_rational D1,D2,D3;
     D1=-(n1*m1.to_vector()); D2=-(n2*m1.to_vector()); D3=-(n3*m2.to_vector());

     //compute intersection of the three planes
     pneu=intersec(n1,n2,n3,D1,D2,D3);
     return pneu;        
   }  
  
};

class Construct_leda_d3_rat_cross_product_vector {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rat_vector       result_type;
 
  leda_rat_vector operator()(const leda_rat_vector& a, const leda_rat_vector& b) const
  {
    CGAL_precondition(a.dim() == 3);
    CGAL_precondition(b.dim() == 3);    
    
    leda_rational a1 = a.xcoord();
    leda_rational a2 = a.ycoord();
    leda_rational a3 = a.zcoord();
    leda_rational b1 = b.xcoord();
    leda_rational b2 = b.ycoord();
    leda_rational b3 = b.zcoord();    
    
    leda_rational v1 = a2*b3-a3*b2;
    leda_rational v2 = a3*b1-a1*b3;
    leda_rational v3 = a1*b2-a2*b1;
    
    return leda_rat_vector(v1,v2,v3);
  }
};

class Construct_leda_d3_rat_opposite_direction {
public:
  typedef Arity_tag< 1 > Arity;
  typedef LEDA_NAMESPACE_NAME::rat_direction       result_type;

    LEDA_NAMESPACE_NAME::rat_direction operator()(const LEDA_NAMESPACE_NAME::rat_direction& dir) const
    { leda_rat_vector v = dir.get_vector(); 
      return LEDA_NAMESPACE_NAME::rat_direction(-v);
    }
};

class Construct_leda_d3_rat_opposite_segment {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_segment       result_type;

    leda_d3_rat_segment operator()(const leda_d3_rat_segment& s) const
    { return s.reverse(); }
};


class Construct_leda_d3_rat_opposite_ray {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_ray       result_type;

    leda_d3_rat_ray operator()(const leda_d3_rat_ray& r) const
    { return r.reverse(); }
};

class Construct_leda_d3_rat_opposite_line {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_line       result_type;

    leda_d3_rat_line operator()(const leda_d3_rat_line& l) const
    { return l.reverse(); }
};

class Construct_leda_d3_rat_opposite_plane {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_d3_rat_plane       result_type;

  leda_d3_rat_plane operator()(const leda_d3_rat_plane& pl) const
  {
    leda_d3_rat_point p1 = pl.point1();
    leda_d3_rat_point p2 = pl.point2();
    leda_d3_rat_point p3 = pl.point3();

    return leda_d3_rat_plane(p3,p2,p1);
  }
};

class Construct_leda_d3_rat_opposite_sphere {
public:
  typedef Arity_tag< 1 > Arity;
  
#if defined(CGAL_COMPATIBLE_SPHERES)  
  typedef LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere       result_type;

   LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s) const
   {
    leda_d3_rat_point p  = s.center();
    leda_rational  sq = s.sqr_radius();
    CGAL::Orientation ori = s.orientation();
    
    ori = reverse_orientation(ori);
  
    return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(p,sq,ori);   
   }
#else
  typedef leda_d3_rat_sphere       result_type;

   leda_d3_rat_sphere operator()(const leda_d3_rat_sphere& s) const
   {
    leda_d3_rat_point p1 = s.point1();
    leda_d3_rat_point p2 = s.point2();
    leda_d3_rat_point p3 = s.point3();
    leda_d3_rat_point p4 = s.point4();   
    
    return leda_d3_rat_sphere(p4,p3,p2,p1);
   }
#endif   
};

class Construct_leda_d3_rat_opposite_vector {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rat_vector       result_type;
  
    leda_rat_vector operator()(const leda_rat_vector& v) const
    { return -v; }
};


CGAL_END_NAMESPACE

#endif




