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
  typedef leda_rational  result_type;
  
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

template<class K>
class Compute_leda_d3_rat_squared_area {

  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3; 
  typedef typename K::Triangle_3   Triangle_3;
  typedef typename K::FT           FT;
  typedef typename K::RT           RT;  

public:
  typedef Arity_tag< 1 > Arity;
  typedef FT  result_type;

  FT operator()(const Triangle_3& t) const
  {
    Point_3 p1 = t.point1();
    Point_3 p2 = t.point2();
    Point_3 p3 = t.point3();
    
    Vector_3 a = p2 - p1;
    Vector_3 b = p3 - p1;

    // compute the cross product ...
    FT a1 = a.xcoord();
    FT a2 = a.ycoord();
    FT a3 = a.zcoord();
    FT b1 = b.xcoord();
    FT b2 = b.ycoord();
    FT b3 = b.zcoord();    
    
    FT vx = a2*b3-a3*b2;
    FT vy = a3*b1-a1*b3;
    FT vz = a1*b2-a2*b1;
    
    Vector_3 cp = Vector_3(vx,vy,vz);    
    FT   scal_prod = cp*cp;
    
    // divide by four (by multiplying denominator)...
    // why division by 4 ?  --> we want squared value 
    // to  compute the size of the parallelogram spanned by a,b we would not divide,
    // for the triangle we would divide by 2 ...
    
    RT n = scal_prod.numerator();
    RT d = scal_prod.denominator() * RT(4);
    return FT(n,d);     
  }
};

/*
class Compute_leda_d3_rat_squared_distance {
public:
   typedef Arity_tag< 2 > Arity;
   typedef leda_rational  result_type;

   leda_rational operator()(const Point_3& p1, const Point_3& p2) const
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

template<class K>
class Compute_leda_d3_rat_squared_length {

  typedef typename K::FT        FT;
  typedef typename K::Segment_3 Segment_3;

public:
  typedef Arity_tag< 1 > Arity;
  typedef FT  result_type;

  FT operator()(const Segment_3& s) const
  { return s.sqr_length(); }   
};

template<class K>
class Compute_leda_d3_rat_squared_radius {

  typedef typename K::Point_3   Point_3;
  typedef typename K::Plane_3   Plane_3;  
  typedef typename K::Vector_3  Vector_3; 
  typedef typename K::FT        FT;  

public:
  typedef FT           result_type;
  typedef Arity_tag< 1 > Arity;

#if defined(CGAL_COMPATIBLE_SPHERES)
  FT operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s) const
  {
    return s.sqr_radius();
  }
#else
   FT operator()(const leda_d3_rat_sphere& s) const
   {
    return s.sqr_radius();
   }
#endif
  
   FT operator()(const Point_3& p1, const Point_3& p2,
                 const Point_3& p3, const Point_3& p4) const
   {
    CGAL_precondition( ! LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p4));
    leda_d3_rat_sphere s(p1,p2,p3,p4);
    return s.sqr_radius();
   }

   // we have similar stuff in Construct_leda_d3_rat_circumcenter ...

   // special version: 3 points, center has to be in the plane defined  by them
   // no precondition documented, but probably collinearity is forbidden ??
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
  
  
  FT operator()(const Point_3& a, const Point_3& b,const Point_3& c) const
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
     // this is the center of the circle defined by these points ...
     pneu=intersec(n1,n2,n3,D1,D2,D3);
     
     return a.sqr_dist(pneu);     
  }  
    
};

template<class K>
class Compute_leda_d3_rat_volume {

  typedef typename K::FT            FT;
  typedef typename K::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename K::Tetrahedron_3 Tetrahedron_3;   

public:
  typedef Arity_tag< 1 > Arity;
  typedef FT  result_type;

  FT operator()(const Iso_cuboid_3& c) const
  {
     FT xd = c.xmax() - c.xmin();
     FT yd = c.ymax() - c.ymin();
     FT zd = c.zmax() - c.zmin();
     
     return xd*yd*zd;
  }
   
  FT operator()(const Tetrahedron_3& s) const
  { return s.vol(); }   
};

CGAL_END_NAMESPACE

#endif




