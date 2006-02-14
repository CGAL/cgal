#ifndef CEP_LEDA_RAT_SUPPORT_FUNCTIONS_3_H
#define CEP_LEDA_RAT_SUPPORT_FUNCTIONS_3_H

#include <CGAL/basic.h>

#include <LEDA/d3_rat_point.h>
#include <LEDA/d3_rat_sphere.h>
#include <LEDA/d3_rat_plane.h>
#include <LEDA/d3_rat_triangle.h>
#include <LEDA/d3_rat_simplex.h>
#include <LEDA/rat_triangle.h>

// some 3d support functions  ...

namespace leda_support {

// ------------------------------------------------------------------------------------------
// bounded side computation for d3_rat_simplex ...
// ------------------------------------------------------------------------------------------


inline CGAL::Bounded_side bounded_side(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, 
                         const leda_d3_rat_point& p)
{
     // inside returns true if p is inside or on s ...
     
     bool inside = s.in_simplex(p);
     if (! inside) return CGAL::ON_UNBOUNDED_SIDE;
     
     // is the point ON the simplex ???
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();   
     
     if (LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p) ||  LEDA_NAMESPACE_NAME::coplanar(p1,p2,p4,p) ||
         LEDA_NAMESPACE_NAME::coplanar(p2,p3,p4,p) ||  LEDA_NAMESPACE_NAME::coplanar(p3,p1,p4,p)) return CGAL::ON_BOUNDARY; 
     
     return CGAL::ON_BOUNDED_SIDE;
}    

// ------------------------------------------------------------------------------------------
// construction of a 3d circle center ...
// ------------------------------------------------------------------------------------------

inline leda_rational r_det(leda_rational a11,leda_rational a12,leda_rational a13,
                    leda_rational a21,leda_rational a22,leda_rational a23,
                    leda_rational a31,leda_rational a32,leda_rational a33)
{
     return a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a11*a23*a32-a12*a21*a33;
}    
   
   
inline leda_d3_rat_point intersection_of_planes(leda_rat_vector n1,leda_rat_vector n2,leda_rat_vector n3,
                                         leda_rational D1,leda_rational D2,leda_rational D3)
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


inline leda_d3_rat_point construct_circle_center_3(const leda_d3_rat_point& a,const leda_d3_rat_point& b,
                                            const leda_d3_rat_point& c)
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
     pneu=intersection_of_planes(n1,n2,n3,D1,D2,D3);
     return pneu;         
}

// ------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------
// projections of 3d triangle ...
// ------------------------------------------------------------------------------------------

inline leda_rat_triangle project_xy(const leda_d3_rat_triangle& t)
{
    leda_rat_point p1 = t.point1().project_xy();
    leda_rat_point p2 = t.point2().project_xy();    
    leda_rat_point p3 = t.point2().project_xy();        
    
    return leda_rat_triangle(p1,p2,p3);
}

inline leda_rat_triangle project_xz(const leda_d3_rat_triangle& t)
{
    leda_rat_point p1 = t.point1().project_xz();
    leda_rat_point p2 = t.point2().project_xz();    
    leda_rat_point p3 = t.point2().project_xz();        
    
    return leda_rat_triangle(p1,p2,p3);
}

inline leda_rat_triangle project_yz(const leda_d3_rat_triangle& t)
{
    leda_rat_point p1 = t.point1().project_yz();
    leda_rat_point p2 = t.point2().project_yz();    
    leda_rat_point p3 = t.point2().project_yz();        
    
    return leda_rat_triangle(p1,p2,p3);
}

// ------------------------------------------------------------------------------------------


}


#endif



