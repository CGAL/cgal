// conversion from CGAL to LEDA kernel (3d) and other direction ...

#ifndef CEP_LEDA_RAT_CONVERSIONS_3_H
#define CEP_LEDA_RAT_CONVERSIONS_3_H

#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>

// LEDA types ...
#include <LEDA/d3_rat_point.h>
#include <LEDA/d3_rat_segment.h>
#include <LEDA/d3_rat_line.h>
#include <LEDA/d3_rat_sphere.h>
#include <LEDA/d3_rat_ray.h>
#include <LEDA/d3_rat_triangle.h>
#include <LEDA/d3_rat_simplex.h>
#include <LEDA/d3_rat_plane.h>
#include <LEDA/rat_vector.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/rat_direction.h>
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/d3_rat_iso_cuboid.h>

#if defined(CGAL_COMPATIBLE_SPHERES)
#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/cgal_d3_rat_sphere.h>
#endif


CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------------------
//  ... to CGAL:
//-----------------------------------------------------------------------------------


struct leda_to_cgal_3 {

  typedef CGAL::Homogeneous<leda_integer> HELP_KERNEL_3;
 

  HELP_KERNEL_3::Point_3 operator()(const leda_d3_rat_point& p) const
  { return HELP_KERNEL_3::Point_3(p.X(), p.Y(), p.Z(), p.W()); }

  HELP_KERNEL_3::Segment_3 operator()(const leda_d3_rat_segment& s)  const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(s.source());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(s.target());
    return HELP_KERNEL_3::Segment_3(p1,p2); 
  }

  HELP_KERNEL_3::Line_3 operator()(const leda_d3_rat_line& l) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(l.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(l.point2());
    return HELP_KERNEL_3::Line_3(p1,p2); 
  }

  HELP_KERNEL_3::Ray_3 operator()(const leda_d3_rat_ray& r) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(r.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(r.point2());
    return HELP_KERNEL_3::Ray_3(p1,p2); 
  }

  HELP_KERNEL_3::Plane_3 operator()(const leda_d3_rat_plane& pl) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(pl.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(pl.point2());
    HELP_KERNEL_3::Point_3 p3 = this->operator()(pl.point3());
    
    return HELP_KERNEL_3::Plane_3(p1,p2,p3); 
  }

  HELP_KERNEL_3::Triangle_3 operator()(const leda_d3_rat_triangle& t) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(t.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(t.point2());
    HELP_KERNEL_3::Point_3 p3 = this->operator()(t.point3());
    
    return HELP_KERNEL_3::Triangle_3(p1,p2,p3); 
  }

  HELP_KERNEL_3::Tetrahedron_3 operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(s.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(s.point2());
    HELP_KERNEL_3::Point_3 p3 = this->operator()(s.point3());
    HELP_KERNEL_3::Point_3 p4 = this->operator()(s.point4());
    
    return HELP_KERNEL_3::Tetrahedron_3(p1,p2,p3,p4); 
  }

  HELP_KERNEL_3::Vector_3 operator()(const leda_rat_vector& v) const
  {     
    return HELP_KERNEL_3::Vector_3(v.X(),v.Y(),v.Z(),v.W()); 
  }

  HELP_KERNEL_3::Direction_3 operator()(const LEDA_NAMESPACE_NAME::rat_direction& d) const
  {     
    HELP_KERNEL_3::Vector_3 v = this->operator()(d.get_vector());
  
    return HELP_KERNEL_3::Direction_3(v); 
  }

  HELP_KERNEL_3::Iso_cuboid_3 operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& i) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(i.vertex(0));
    HELP_KERNEL_3::Point_3 p2 = this->operator()(i.vertex(7));
    
    return HELP_KERNEL_3::Iso_cuboid_3(p1,p2); 
  }

  HELP_KERNEL_3::Sphere_3 operator()(const leda_d3_rat_sphere& s) const
  { 
    HELP_KERNEL_3::Point_3 p1 = this->operator()(s.point1());
    HELP_KERNEL_3::Point_3 p2 = this->operator()(s.point2());
    HELP_KERNEL_3::Point_3 p3 = this->operator()(s.point3());
    HELP_KERNEL_3::Point_3 p4 = this->operator()(s.point4());    
    
    return HELP_KERNEL_3::Sphere_3(p1,p2,p3,p4); 
  }


#if defined(CGAL_COMPATIBLE_SPHERES)
  HELP_KERNEL_3::Sphere_3 operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s) const
  { 
    HELP_KERNEL_3::Point_3 c = this->operator()(s.center());
    leda_rational        rad = s.sqr_radius();
    
    return HELP_KERNEL_3::Sphere_3(c, CGAL::Quotient<leda_integer>(rad.numerator(),rad.denominator())); 
  }
#endif

  // object conversion; the Object stores a leda 3d object ...
  CGAL::Object  operator()(const CGAL::Object& leda_obj) const
  {
    leda_d3_rat_point p;
     
    if (CGAL::assign(p, leda_obj)){
      return CGAL::make_object(this->operator()(p));
    }     
     
    leda_d3_rat_segment s;
     
    if (CGAL::assign(s, leda_obj)){
      return CGAL::make_object(this->operator()(s));
    }     
     
    leda_d3_rat_line l;
     
    if (CGAL::assign(l, leda_obj)){
      return CGAL::make_object(this->operator()(l));
    }     
     
    leda_d3_rat_ray r;
     
    if (CGAL::assign(r, leda_obj)){
      return CGAL::make_object(this->operator()(r));
    }     
     
    leda_d3_rat_plane pl;
    
    if (CGAL::assign(pl, leda_obj)){
      return CGAL::make_object(this->operator()(pl));
    }    
     
    leda_d3_rat_triangle t;
    
    if (CGAL::assign(t, leda_obj)){
      return CGAL::make_object(this->operator()(t));
    }    
     
    LEDA_NAMESPACE_NAME::d3_rat_simplex sim;
    
    if (CGAL::assign(sim, leda_obj)){
      return CGAL::make_object(this->operator()(sim));
    }    
     
    leda_rat_vector v;
    
    if (CGAL::assign(v, leda_obj)){
      return CGAL::make_object(this->operator()(v));
    }    
     
    LEDA_NAMESPACE_NAME::rat_direction dir;
    
    if (CGAL::assign(dir, leda_obj)){
      return CGAL::make_object(this->operator()(dir));
    }    
     
    LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid i;
    
    if (CGAL::assign(i, leda_obj)){
      return CGAL::make_object(this->operator()(i));
    }    
     
    leda_d3_rat_sphere sph;
    
    if (CGAL::assign(sph, leda_obj)){
      return CGAL::make_object(this->operator()(sph));
    }    

#if defined(CGAL_COMPATIBLE_SPHERES)     
    LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere  cg_sph;
    
    if (CGAL::assign(cg_sph, leda_obj)){
      return CGAL::make_object(this->operator()(cg_sph));
    }    
#endif
  
    return  leda_obj;    
  }

};

//-----------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------
//  ... to LEDA:
//-----------------------------------------------------------------------------------


struct cgal_to_leda_3 {

  typedef CGAL::Homogeneous<leda_integer> HELP_KERNEL_3;

  leda_d3_rat_point operator()(const HELP_KERNEL_3::Point_3& p) const
  { return leda_d3_rat_point(p.hx(), p.hy(), p.hz(), p.hw()); }

  leda_d3_rat_segment operator()(const HELP_KERNEL_3::Segment_3& s) const
  {    
    leda_d3_rat_point p1 = this->operator()(s.source());
    leda_d3_rat_point p2 = this->operator()(s.target());
    return leda_d3_rat_segment(p1,p2); 
  }

  leda_d3_rat_line operator()(const HELP_KERNEL_3::Line_3& l) const
  { 
    leda_d3_rat_point p1 = this->operator()(l.point(0));
    leda_d3_rat_point p2 = this->operator()(l.point(1));
    return leda_d3_rat_line(p1,p2); 
  }

  leda_d3_rat_ray operator()(const HELP_KERNEL_3::Ray_3& r) const
  {    
    leda_d3_rat_point p1 = this->operator()(r.point(0));
    leda_d3_rat_point p2 = this->operator()(r.point(1));
    return leda_d3_rat_ray(p1,p2); 
  }

  leda_rat_vector operator()(const HELP_KERNEL_3::Vector_3& v) const
  {     
    return leda_rat_vector(v.hx(),v.hy(),v.hz(),v.hw()); 
  }

  leda_d3_rat_plane operator()(const HELP_KERNEL_3::Plane_3& pl) const
  {   
    // get a point on pl ...
    leda_d3_rat_point p = this->operator()(pl.point());
    leda_rat_vector   v = this->operator()(pl.orthogonal_vector());
    
    return leda_d3_rat_plane(p,v); 
  }

  leda_d3_rat_triangle operator()(const HELP_KERNEL_3::Triangle_3& t) const
  {    
    leda_d3_rat_point p1 = this->operator()(t.vertex(0));
    leda_d3_rat_point p2 = this->operator()(t.vertex(1));
    leda_d3_rat_point p3 = this->operator()(t.vertex(2));
    
    return leda_d3_rat_triangle(p1,p2,p3); 
  }

  LEDA_NAMESPACE_NAME::d3_rat_simplex operator()(const HELP_KERNEL_3::Tetrahedron_3& t) const
  { 
    leda_d3_rat_point p1 = this->operator()(t.vertex(0));
    leda_d3_rat_point p2 = this->operator()(t.vertex(1));
    leda_d3_rat_point p3 = this->operator()(t.vertex(2));
    leda_d3_rat_point p4 = this->operator()(t.vertex(3));
    
    return LEDA_NAMESPACE_NAME::d3_rat_simplex(p1,p2,p3,p4); 
  }

  LEDA_NAMESPACE_NAME::rat_direction operator()(const HELP_KERNEL_3::Direction_3& d) const
  {     
    leda_rat_vector v = this->operator()(d.vector());
  
    return LEDA_NAMESPACE_NAME::rat_direction(v); 
  }

  LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid operator()(const HELP_KERNEL_3::Iso_cuboid_3& i) const
  { 
    leda_d3_rat_point p1 = this->operator()(i.min());
    leda_d3_rat_point p2 = this->operator()(i.max());
    
    return LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid(p1,p2); 
  }

  // conversion to leda spheres does not work ...

#if defined(CGAL_COMPATIBLE_SPHERES)
  LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere operator()(const HELP_KERNEL_3::Sphere_3& s) const
  { 
    leda_d3_rat_point center = this->operator()(s.center());
    CGAL::Quotient<leda_integer>  q = s.squared_radius();
    leda_rational  rad(q.numerator(), q.denominator() );
    
    return LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere(center,rad); 
  }
#endif

  // object conversion; the Object stores a cgal 3d object ...
  CGAL::Object  operator()(const CGAL::Object& cgal_obj) const
  {
    HELP_KERNEL_3::Point_3 p;
    
    if (CGAL::assign(p, cgal_obj)){
      return CGAL::make_object(this->operator()(p));
    }
    
    HELP_KERNEL_3::Segment_3 s;
    
    if (CGAL::assign(s, cgal_obj)){
      return CGAL::make_object(this->operator()(s));
    }    
    
    HELP_KERNEL_3::Line_3 l;
    
    if (CGAL::assign(l, cgal_obj)){
      return CGAL::make_object(this->operator()(l));
    }    
    
    HELP_KERNEL_3::Ray_3 r;
    
    if (CGAL::assign(r, cgal_obj)){
      return CGAL::make_object(this->operator()(r));
    }    
    
    HELP_KERNEL_3::Triangle_3 t;
    
    if (CGAL::assign(t, cgal_obj)){
      return CGAL::make_object(this->operator()(t));
    } 
    
    HELP_KERNEL_3::Vector_3 vec;
    
    if (CGAL::assign(vec, cgal_obj)){
      return CGAL::make_object(this->operator()(vec));
    }     
    
    HELP_KERNEL_3::Plane_3 pl;
    
    if (CGAL::assign(pl, cgal_obj)){
      return CGAL::make_object(this->operator()(pl));
    }     
    
    HELP_KERNEL_3::Triangle_3 tr;
    
    if (CGAL::assign(tr, cgal_obj)){
      return CGAL::make_object(this->operator()(tr));
    }     
    
    HELP_KERNEL_3::Tetrahedron_3  tetr;
    
    if (CGAL::assign(tetr, cgal_obj)){
      return CGAL::make_object(this->operator()(tetr));
    }     
    
    HELP_KERNEL_3::Direction_3  dir;
    
    if (CGAL::assign(dir, cgal_obj)){
      return CGAL::make_object(this->operator()(dir));
    }     
    
    HELP_KERNEL_3::Iso_cuboid_3  i;
    
    if (CGAL::assign(i, cgal_obj)){
      return CGAL::make_object(this->operator()(i));
    }     
    
#if defined(CGAL_COMPATIBLE_SPHERES)
    HELP_KERNEL_3::Sphere_3  sph;
    
    if (CGAL::assign(sph, cgal_obj)){
      return CGAL::make_object(this->operator()(sph));
    }     
#endif

    return cgal_obj;    
  }     

};


CGAL_END_NAMESPACE

#endif


