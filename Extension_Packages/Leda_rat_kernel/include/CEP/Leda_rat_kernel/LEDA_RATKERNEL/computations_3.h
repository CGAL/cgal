#ifndef CEP_LEDA_RAT_COMPUTATIONS_3_H
#define CEP_LEDA_RAT_COMPUTATIONS_3_H

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

// LEDA 3d rational kernel computation objects ...

// 3d computations ...


CGAL_BEGIN_NAMESPACE

template<class HELP_KERNEL>
class CGAL_compute_leda_d3_rat_squared_distance {
public:
  typedef Arity_tag< 2 > Arity;
  typedef leda_rational           result_type;
  
  typedef typename HELP_KERNEL::Compute_squared_distance_3  Compute_squared_distance_3;

  template<class T1, class T2>
  leda_rational operator()(const T1& obj1, const T2& obj2) const
  {
     leda_to_cgal_3 conv;
     Compute_squared_distance_3 compute;
     
     CGAL::Quotient<leda_integer> result = compute(conv(obj1), conv(obj2));
     return leda_rational(result.numerator(), result.denominator());  
  }
};


class Compute_leda_d3_rat_squared_area {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rational  result_type;

  leda_rational operator()(const leda_d3_rat_triangle& t) const
  {
    leda_d3_rat_point p1 = t.point1();
    leda_d3_rat_point p2 = t.point2();
    leda_d3_rat_point p3 = t.point3();
    
    leda_rat_vector a = p2 - p1;
    leda_rat_vector b = p3 - p1;

    // compute the cross product ...
    leda_rational a1 = a.xcoord();
    leda_rational a2 = a.ycoord();
    leda_rational a3 = a.zcoord();
    leda_rational b1 = b.xcoord();
    leda_rational b2 = b.ycoord();
    leda_rational b3 = b.zcoord();    
    
    leda_rational vx = a2*b3-a3*b2;
    leda_rational vy = a3*b1-a1*b3;
    leda_rational vz = a1*b2-a2*b1;
    
    leda_rat_vector cp = leda_rat_vector(vx,vy,vz);    
    leda_rational   scal_prod = cp*cp;
    
    // divide by four (by multiplying denominator)...
    // why division by 4 ?  --> we want squared value 
    // to  compute the size of the parallelogram spanned by a,b we would not divide,
    // for the triangle we would divide by 2 ...
    
    leda_integer n = scal_prod.numerator();
    leda_integer d = scal_prod.denominator() * leda_integer(4);
    return leda_rational(n,d);     
  }
};

/*
class Compute_leda_d3_rat_squared_distance {
public:
   typedef Arity_tag< 2 > Arity;
   typedef leda_rational  result_type;

   leda_rational operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
      return p1.sqr_dist(p2);
   }

   leda_rational operator()(const leda_d3_rat_point& p1, const leda_d3_rat_line& l) const
   { return l.sqr_dist(p1); }
   
   leda_rational operator()(const leda_d3_rat_line& l, const leda_d3_rat_point& p1) const
   { return l.sqr_dist(p1); }  
   
   leda_rational operator()(const leda_d3_rat_point& p, const leda_d3_rat_segment& s) const
   { 
     leda_d3_rat_point p1 = s.source();
     leda_d3_rat_point p2 = s.target();
     
     if (p1==p2){
       return p.sqr_dist(p1);
     }
     // a real segment ...
     leda_d3_rat_line l(p1,p2);
     
     leda_rational d1 = l.sqr_dist(p);
     leda_rational d2 = p.sqr_dist(p1);
     leda_rational d3 = p.sqr_dist(p2);
     
     if (d2 < d1) d1=d2;
     if (d3 < d1) d1=d3;
     
     return d1;
   }
   
   leda_rational operator()(const leda_d3_rat_segment& s, const leda_d3_rat_point& p1) const
   { return this->operator()(p1,s); }       

   leda_rational operator()(const leda_d3_rat_point& p, const leda_d3_rat_ray& r) const
   { 
     leda_d3_rat_point p1 = r.source();
     
     leda_d3_rat_line l(p1,r.point2());
     
     leda_rational d1 = l.sqr_dist(p);
     leda_rational d2 = p.sqr_dist(p1);
     
     if (d2 < d1) d1=d2;
     
     return d1;
   }
   
   leda_rational operator()(const leda_d3_rat_ray& r, const leda_d3_rat_point& p1) const
   { return this->operator()(p1,r); }       
   
   leda_rational operator()(const leda_d3_rat_point& p1, const leda_d3_rat_plane& pl) const
   {
      return pl.sqr_dist(p1);
   }

   leda_rational operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p1) const
   {
      return pl.sqr_dist(p1);
   }     
   
   // line
   
   // ray
   
   // segment
   
   // plane
   leda_rational operator()(const leda_d3_rat_plane& pl1, const leda_d3_rat_plane& pl2) const    
   {
      if (! pl1.parallel(pl2)) return leda_rational(leda_integer(0), leda_integer(1));
      
      // the planes are parallel ...
      leda_d3_rat_point p1 = pl1.point1();
      
      // normal projection ...
      leda_rat_vector vec = pl2.normal_project(p1);
      
      return vec.sqr_length();
   }   
};
*/


class Compute_leda_d3_rat_squared_length {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rational  result_type;

  leda_rational operator()(const leda_d3_rat_segment& s) const
  {
     return s.sqr_length();
  }   
};

class Compute_leda_d3_rat_squared_radius {
public:
  typedef leda_rational           result_type;
  typedef Arity_tag< 1 > Arity;

#if defined(CGAL_COMPATIBLE_SPHERES)
  leda_rational operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s) const
  {
    return s.sqr_radius();
  }
#else
   leda_rational operator()(const leda_d3_rat_sphere& s) const
   {
    return s.sqr_radius();
   }
#endif
  
   leda_rational operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                            const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
   {
    CGAL_precondition( ! LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p4));
    leda_d3_rat_sphere s(p1,p2,p3,p4);
    return s.sqr_radius();
   }

   // we have similar stuff in Construct_leda_d3_rat_circumcenter ...

   // special version: 3 points, center has to be in the plane defined  by them
   // no precondition documented, but probably collinearity is forbidden ??
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
  
  
  leda_rational operator()(const leda_d3_rat_point& a, const leda_d3_rat_point& b,
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
     // this is the center of the circle defined by these points ...
     pneu=intersec(n1,n2,n3,D1,D2,D3);
     
     return a.sqr_dist(pneu);     
  }  
    
};

class Compute_leda_d3_rat_volume {
public:
  typedef Arity_tag< 1 > Arity;
  typedef leda_rational  result_type;

   leda_rational operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& c) const
   {
     leda_rational xd = c.xmax() - c.xmin();
     leda_rational yd = c.ymax() - c.ymin();
     leda_rational zd = c.zmax() - c.zmin();
     
     return xd*yd*zd;
   }
   
   leda_rational operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s) const
   {
     return s.vol();
   }   
};

CGAL_END_NAMESPACE

#endif




