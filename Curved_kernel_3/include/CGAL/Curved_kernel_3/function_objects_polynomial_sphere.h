// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Julien Hazebrouck
//           Damien Leroy
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H
#define CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H

#include <CGAL/kernel_basic.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_circular_arc_point_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_line_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Curved_kernel_3/internal_functions_on_circle_3.h>
#include <CGAL/Curved_kernel_3/internal_function_has_on_spherical_kernel.h>
#include <CGAL/Object.h>


namespace CGAL {
namespace SphericalFunctors {

  template < class SK >
  class Compare_x_3
    : public SK::Linear_kernel::Compare_x_3
  {
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Point_3 Point_3;

  public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>             Arity;

    using SK::Linear_kernel::Compare_x_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_x<SK>(p0, p1);}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Point_3 &p1) const
    { return compare_x<SK>(p0, p1);}

    result_type
    operator() (const Point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_x<SK>(p0, p1);}

  };

  template < class SK >
  class Compare_y_3
    : public SK::Linear_kernel::Compare_y_3
  {
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Point_3 Point_3;

  public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>             Arity;

    using SK::Linear_kernel::Compare_x_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_y<SK>(p0, p1);}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Point_3 &p1) const
    { return compare_y<SK>(p0, p1);}

    result_type
    operator() (const Point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_y<SK>(p0, p1);}

  };

  template < class SK >
  class Compare_z_3
    : public SK::Linear_kernel::Compare_z_3
  {
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Point_3 Point_3;

  public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>             Arity;

    using SK::Linear_kernel::Compare_x_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_z<SK>(p0, p1);}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Point_3 &p1) const
    { return compare_z<SK>(p0, p1);}

    result_type
    operator() (const Point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_z<SK>(p0, p1);}

  };

  template < class SK >
  class Compare_xy_3
    : public SK::Linear_kernel::Compare_xy_3
  {
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Point_3 Point_3;

  public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>             Arity;

    using SK::Linear_kernel::Compare_xy_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_xy<SK>(p0, p1);}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Point_3 &p1) const
    { return compare_xy<SK>(p0, p1);}

    result_type
    operator() (const Point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_xy<SK>(p0, p1);}

  };

template < class SK >
  class Compare_xyz_3
    : public SK::Linear_kernel::Compare_xyz_3
  {
    typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
    typedef typename SK::Point_3 Point_3;

  public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>             Arity;

    using SK::Linear_kernel::Compare_xyz_3::operator();

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_xyz<SK>(p0, p1);}

    result_type
    operator() (const Circular_arc_point_3 &p0,
                const Point_3 &p1) const
    { return compare_xyz<SK>(p0, p1);}

    result_type
    operator() (const Point_3 &p0,
                const Circular_arc_point_3 &p1) const
    { return compare_xyz<SK>(p0, p1);}

  };

 template <class SK>
  class Compute_circular_x_3: Has_qrt
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_3 & a) const
    {
      return (a.rep().x());
    }
  };

template <class SK>
  class Compute_circular_y_3: Has_qrt
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_3 & a) const
    {
      return (a.rep().y());
    }
  };

template <class SK>
  class Compute_circular_z_3: Has_qrt
  {
    typedef typename SK::Circular_arc_point_3   Circular_arc_point_3;
    typedef typename SK::Root_of_2                 Root_of_2;

  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    typedef Arity_tag<1>             Arity;

    qualified_result_type operator() (const Circular_arc_point_3 & a) const
    {
      return (a.rep().z());
    }
  };

template < class SK >
  class Equal_3
    : public SK::Linear_kernel::Equal_3
  {
    typedef typename SK::Circular_arc_point_3     Circular_arc_point_3;
    typedef typename SK::Circle_3                 Circle_3;
  public:
    typedef bool result_type;
    typedef Arity_tag<2>             Arity;
    
    using SK::Linear_kernel::Equal_3::operator();

    // Our Circle_3 dont have orientation
    result_type
    operator() (const Circle_3 &c0,
                const Circle_3 &c1) const
    { return equal<SK>(c0, c1); }

    result_type
    operator() (const Circular_arc_point_3 &c0,
                const Circular_arc_point_3 &c1) const
    { return equal<SK>(c0, c1); }

  };

template < class SK >
  class Construct_circular_arc_point_3
  {
    typedef typename SK::Point_3                            Point_3;
    typedef typename SK::Plane_3                            Plane_3;
    typedef typename SK::Line_3                             Line_3;
    typedef typename SK::Circle_3                           Circle_3;
    typedef typename SK::Sphere_3                           Sphere_3;
    typedef typename SK::Circular_arc_point_3               Circular_arc_point_3;
    typedef typename SK::Kernel_base::Circular_arc_point_3  RCircular_arc_point_3;
    typedef typename SK::Root_of_2                          Root_of_2;
    typedef typename Circular_arc_point_3::Rep              Rep;
    typedef typename Circular_arc_point_3::Root_for_spheres_2_3  Root_for_spheres_2_3;

  public:
    typedef  Circular_arc_point_3 result_type;
    typedef Arity_tag<1>             Arity;


    result_type
    operator()(void) 
    { return Rep(); }
    
    result_type
      operator()(const Root_of_2 & x,
		 const Root_of_2 & y,
		 const Root_of_2 & z
		 ) const
    { return Rep(x,y,z); }
    
    result_type
    operator()(const Root_for_spheres_2_3 & np) const
    { return Rep(np); }

    result_type
    operator()(const Point_3 & p) const
    { return Rep(p); }

    result_type
    operator()(const Sphere_3 & s1,
               const Sphere_3 & s2,
               const Sphere_3 & s3,
               const bool less_xyz = true) const
    { return Rep(s1,s2,s3,less_xyz); }

    result_type
    operator()(const Plane_3 & p,
               const Sphere_3 & s1,
               const Sphere_3 & s2,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    result_type
    operator()(const Sphere_3 & s1,
               const Plane_3 & p,
               const Sphere_3 & s2,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    result_type
    operator()(const Sphere_3 & s1,
               const Sphere_3 & s2,
               const Plane_3 & p,
               const bool less_xyz = true) const
    { return Rep(p,s1,s2,less_xyz); }

    result_type
    operator()(const Plane_3 & p1,
               const Plane_3 & p2,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    result_type
    operator()(const Plane_3 & p1,
               const Sphere_3 & s,
               const Plane_3 & p2,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    result_type
    operator()(const Sphere_3 & s,
               const Plane_3 & p1,
               const Plane_3 & p2,
               const bool less_xyz = true) const
    { return Rep(p1,p2,s,less_xyz); }

    result_type
    operator()(const Line_3 & l,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(l,s,less_xyz); }

    result_type
    operator()(const Sphere_3 & s,
               const Line_3 & l,
               const bool less_xyz = true) const
    { return Rep(l,s,less_xyz); }

    result_type
    operator()(const Circle_3 & c,
               const Sphere_3 & s,
               const bool less_xyz = true) const
    { return Rep(c,s,less_xyz); }

    result_type
    operator()(const Sphere_3 & s,
               const Circle_3 & c,
               const bool less_xyz = true) const
    { return Rep(c,s,less_xyz); }

    result_type
    operator()(const Circle_3 & c,
               const Plane_3 & p,
               const bool less_xyz = true) const
    { return Rep(c,p,less_xyz); }

    result_type
    operator()(const Plane_3 & p,
               const Circle_3 & c,
               const bool less_xyz = true) const
    { return Rep(c,p,less_xyz); }

  };

template < class SK >
  class Construct_sphere_3 : public  SK::Linear_kernel::Construct_sphere_3
  {
  public:
    
    typedef typename SK::Sphere_3 result_type;
    typedef Arity_tag<1>          Arity; 

    using SK::Linear_kernel::Construct_sphere_3::operator();

    result_type
    operator() ( const typename SK::Polynomial_for_spheres_2_3 &eq )
    {
      return construct_sphere_3<SK>(eq);
    }
  };

template < class SK >
  class Construct_plane_3 : public  SK::Linear_kernel::Construct_plane_3
  {
  public:
    
    typedef typename SK::Plane_3 result_type;
    typedef Arity_tag<1>          Arity; 

    using SK::Linear_kernel::Construct_plane_3::operator();

    result_type
    operator() ( const typename SK::Polynomial_1_3 &eq )
    {
      return construct_plane_3<SK>(eq);
    }
  };

template < class SK >
  class Construct_line_3 : public  SK::Linear_kernel::Construct_line_3
  {
  public:
    
    typedef typename SK::Line_3 result_type;
    typedef Arity_tag<1>          Arity; 

    using SK::Linear_kernel::Construct_line_3::operator();

    result_type
    operator() ( const typename SK::Polynomials_for_line_3 &eq )
    {
      return construct_line_3<SK>(eq);
    }
  };

template < class SK >
  class Construct_circle_3
  {
  public:

    typedef typename SK::FT FT;
    typedef typename SK::Point_3 Point_3;
    typedef typename SK::Plane_3 Plane_3;
    typedef typename SK::Sphere_3 Sphere_3;
    typedef typename SK::Circle_3 Circle_3;
    typedef typename SK::Kernel_base::Circle_3  RCircle_3;
    typedef typename Circle_3::Rep              Rep;

    typedef Circle_3 result_type;
    typedef Arity_tag<1>          Arity; 

    result_type
    operator()(void) 
    { return Rep(); }

    result_type
    operator()(const Circle_3 &c) const
    { return Rep(c); }

    result_type
    operator() ( const typename SK::Polynomials_for_circle_3 &eq )
    { return Rep(construct_circle_3<SK>(eq)); }

    result_type
    operator() (const Point_3& p, const FT& sr, const Plane_3& plane)
    { return Rep(p,sr,plane); }

    result_type
    operator() (const Sphere_3& s1, const Sphere_3& s2)
    { return Rep(s1,s2); }

    result_type
    operator() (const Plane_3& p, const Sphere_3& s)
    { return Rep(p,s); }

    result_type
    operator() (const Sphere_3& s, const Plane_3& p)
    { return Rep(p,s); }
  };

template <class SK>
  class Construct_supporting_plane_3//: Has_qrt
  {
    typedef typename SK::Plane_3      Plane_3;
    typedef typename SK::Circle_3      Circle_3;

  public:

    typedef Plane_3        result_type;
    typedef const Plane_3& qualified_result_type;
    typedef Arity_tag<1>   Arity;
    
    qualified_result_type operator() (const Circle_3 & c) const
    {
      return c.rep().supporting_plane();
    }
  };

template <class SK>
  class Construct_diametral_sphere_3//: Has_qrt
  {
    typedef typename SK::Sphere_3      Sphere_3;
    typedef typename SK::Circle_3      Circle_3;

  public:

    typedef Sphere_3        result_type;
    typedef const Sphere_3& qualified_result_type;
    typedef Arity_tag<1>    Arity;
    
    qualified_result_type operator() (const Circle_3 & c) const
    {
      return c.rep().diametral_sphere();
    }
  };

template < class SK >
  class Has_on_3
  {
    typedef typename SK::Point_3                 Point_3;
    typedef typename SK::Sphere_3                Sphere_3;
    typedef typename SK::Plane_3                 Plane_3;
    typedef typename SK::Line_3                  Line_3;
    typedef typename SK::Circular_arc_point_3    Circular_arc_point_3;
    typedef typename SK::Circle_3                Circle_3;
    

  public:
    typedef bool result_type;
    typedef Arity_tag<2>             Arity;
    
    result_type
    operator()(const Sphere_3 &a, const Point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Point_3 &p, const Sphere_3 &a) const
    { return false; }

    result_type
    operator()(const Sphere_3 &a, const Circular_arc_point_3 &p) const
    { 
      return has_on<SK>(a, p); 
    }

    result_type
    operator()(const Circular_arc_point_3 &p, const Sphere_3 &a) const
    { return false; }

    result_type
    operator()(const Plane_3 &a, const Point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Point_3 &p, const Plane_3 &a) const
    { return false; }

    result_type
    operator()(const Plane_3 &a, const Circular_arc_point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Circular_arc_point_3 &p, const Plane_3 &a) const
    { return false; }

    result_type
    operator()(const Line_3 &a, const Point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Point_3 &p, const Line_3 &a) const
    { return false; }

    result_type
    operator()(const Line_3 &a, const Circular_arc_point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Circular_arc_point_3 &p, const Line_3 &a) const
    { return false; }

    result_type
    operator()(const Circle_3 &a, const Point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Point_3 &p, const Circle_3 &a) const
    { return false; }

    result_type
    operator()(const Circle_3 &a, const Circular_arc_point_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Circular_arc_point_3 &p, const Circle_3 &a) const
    { return false; }

    result_type
    operator()(const Sphere_3 &a, const Circle_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Circle_3 &p, const Sphere_3 &a) const
    { return false; }

    result_type
    operator()(const Plane_3 &a, const Line_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Line_3 &a, const Plane_3 &p) const
    { return false; }

    result_type
    operator()(const Plane_3 &a, const Circle_3 &p) const
    { return has_on<SK>(a, p); }

    result_type
    operator()(const Circle_3 &p, const Plane_3 &a) const
    { return false; }

    // We assume Spheres cannot be points, should we?
    result_type
    operator()(const Sphere_3 &a, const Plane_3 &p) const
    { return false; }

    result_type
    operator()(const Plane_3 &a, const Sphere_3 &p) const
    { return false; }

    result_type
    operator()(const Sphere_3 &a, const Line_3 &p) const
    { return false; }

  };

template < class SK >
  class Intersect_3
    : public SK::Linear_kernel::Intersect_3
  {
  
    typedef typename SK::Sphere_3                 Sphere_3;
    typedef typename SK::Line_3                   Line_3;
    typedef typename SK::Plane_3                  Plane_3;
    typedef typename SK::Circle_3                 Circle_3;
    
    public:

    typedef void result_type;   
    //typedef Arity_tag<2> Arity; // The Arity can be 2 and 3
                                  // Is there some solution for this problem??
    
    using SK::Linear_kernel::Intersect_3::operator();

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s1, const Sphere_3 & s2, 
	       OutputIterator res) const
      { return intersect_3<SK> (s1,s2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Sphere_3 & s, 
	       OutputIterator res) const
      { return intersect_3<SK> (p,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Plane_3 & p, 
	       OutputIterator res) const
      { return intersect_3<SK> (p,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Line_3 & l, 
	       OutputIterator res) const
      { return intersect_3<SK> (s,l,res); }
    
     template < class OutputIterator >
    OutputIterator
      operator()(const Line_3 & l,const Sphere_3 & s, 
	       OutputIterator res) const
      { return intersect_3<SK> (s,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s1, const Sphere_3 & s2, 
	       const Sphere_3 & s3, OutputIterator res) const
      { return intersect_3<SK> (s1,s2,s3,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s1, const Sphere_3 & s2, 
	       const Plane_3 & p, OutputIterator res) const
      { return intersect_3<SK> (p,s1,s2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Sphere_3 & s1,  
	       const Sphere_3 & s2, OutputIterator res) const
      { return intersect_3<SK> (p,s1,s2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p1, const Plane_3 & p2,  
	       const Sphere_3 & s, OutputIterator res) const
      { return intersect_3<SK> (p1,p2,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Plane_3 & p1, 
	       const Plane_3 & p2, OutputIterator res) const
      { return intersect_3<SK> (p1,p2,s,res); }

    // INTERSECTIONS WITH CIRCLE
    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Plane_3 & p, 
	       OutputIterator res) const
      { return intersect_3<SK> (c,p,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Plane_3 & p, const Circle_3 & c, 
	       OutputIterator res) const
      { return intersect_3<SK> (c,p,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Sphere_3 & s, 
	       OutputIterator res) const
      { return intersect_3<SK> (c,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Sphere_3 & s, const Circle_3 & c, 
	       OutputIterator res) const
      { return intersect_3<SK> (c,s,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c1, const Circle_3 & c2, 
	       OutputIterator res) const
      { return intersect_3<SK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle_3 & c, const Line_3 & l, 
	       OutputIterator res) const
      { return intersect_3<SK> (c,l,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_3 & l, const Circle_3 & c,
	       OutputIterator res) const
      { return intersect_3<SK> (c,l,res); }

  };

} // namespace SphericalFunctors
} // namespace CGAL

#endif // CGAL_SPHERICAL_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_SPHERE_H
