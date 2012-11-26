// Copyright (c) 1997   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://mcaroli@scm.gforge.inria.fr/svn/cgal/trunk/Triangulation_2/include/CGAL/Regular_triangulation_euclidean_traits_2.h $
// $Id: Regular_triangulation_euclidean_traits_2.h 46206 2008-10-11 20:21:08Z spion $
// 
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion
//                 Manuel Caroli

#ifndef CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H

#include <CGAL/Weighted_point.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Triangulation_sphere_traits_2.h>

#include <CGAL/algebraic_kernel_for_spheres_extensions_ntC3.h>

namespace CGAL { 

// TODO: Constructions to be adapted.

// constructions for DUALITY: weighted_circumcenter and radical 
//  axis

// template < class Bare_point, class We >
// inline
// Bare_point
// weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
// 		      const Weighted_point< Bare_point,We >& q,
// 		      const Weighted_point< Bare_point,We >& r,
// 		      Cartesian_tag )
// {
//   typename Kernel_traits<Bare_point>::Kernel::RT x,y;
//   weighted_circumcenterC2(p.x(),p.y(),p.weight(),
// 			  q.x(),q.y(),q.weight(),
// 			  r.x(),r.y(),r.weight(),x,y);
//   return Bare_point(x,y);
// }

// template < class Bare_point, class We >
// inline
// Bare_point
// weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
// 		      const Weighted_point< Bare_point,We >& q,
// 		      const Weighted_point< Bare_point,We >& r,
// 		      Homogeneous_tag )
// {
//   typename Kernel_traits<Bare_point>::Kernel::RT x,y,w;
//   weighted_circumcenterH2(p.hx(),p.hy(),p.hw(),p.weight(),
// 			  q.hx(),q.hy(),q.hw(),q.weight(),
// 			  r.hx(),r.hy(),r.hw(),r.weight(),
// 			  x,y,w);
//   return Bare_point(x,y,w);
// }


// template < class Bare_point, class We >
// inline
// Bare_point
// weighted_circumcenter(const Weighted_point< Bare_point,We >& p,
// 		      const Weighted_point< Bare_point,We >& q,
// 		      const Weighted_point< Bare_point,We >& r)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
//   return weighted_circumcenter(p, q, r, Tag()); 
// }


template < typename K >
class Dummy_constructor
{
public:
  typedef typename K::Weighted_point_2         Weighted_point_2;
  typedef typename K::Bare_point_2               Bare_point;

  typedef Bare_point       result_type;

  Bare_point operator() ( const Weighted_point_2 & p,
		          const Weighted_point_2 & ,
		          const Weighted_point_2 & ,
		          const Weighted_point_2 & ) const
    { return p; }

  Bare_point operator() ( const Weighted_point_2 & p,
		          const Weighted_point_2 & ,
		          const Weighted_point_2 & ) const
    { return p; }
};

// template < typename K >
// class Construct_weighted_circumcenter_2
// {
// public:
//   typedef typename K::Weighted_point_2         Weighted_point_2;
//   typedef typename K::Bare_point               Bare_point;

//   typedef Bare_point       result_type;

//   Bare_point operator() ( const Weighted_point_2 & p,
// 		          const Weighted_point_2 & q,
// 		          const Weighted_point_2 & r) const
//     {
//       CGAL_triangulation_precondition( ! collinear(p, q, r) );
//       return CGAL::weighted_circumcenter(p,q,r);
//     }
// };
 


// template < class Bare_point, class We >
// inline
// Line_2<typename Kernel_traits<Bare_point>::Kernel>
// radical_axis(const Weighted_point< Bare_point,We >& p,
// 	     const Weighted_point< Bare_point,We >& q,
// 	     Cartesian_tag )
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::RT RT;
//   typedef typename Kernel_traits<Bare_point>::Kernel     Rep;
//   RT a,b,c;
//   radical_axisC2(p.x(),p.y(),p.weight(),q.x(),q.y(),q.weight(),a,b,c);
//   return Line_2<Rep>(a,b,c);
// }

// template < class Bare_point, class We >
// inline
// Line_2<typename Kernel_traits<Bare_point>::Kernel>
// radical_axis(const Weighted_point< Bare_point,We >& p,
// 	     const Weighted_point< Bare_point,We >& q,
// 	      Homogeneous_tag)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::RT RT;
//   typedef typename Kernel_traits<Bare_point>::Kernel     Rep;
//   RT a,b,c;
//   radical_axisH2(p.hx(),p.hy(), p.hw(), p.weight(),
// 		 q.hx(),q.hy(), q.hw(), q.weight(),a,b,c);
//   return Line_2<Rep>(a,b,c);
// }

// template < class Bare_point, class We >
// inline
// Line_2<typename Kernel_traits<Bare_point>::Kernel>
// radical_axis(const Weighted_point< Bare_point,We >& p,
// 	     const Weighted_point< Bare_point,We >& q)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
//   return radical_axis(p, q, Tag()); 
// }


// template < typename K >
// class Construct_radical_axis_2
// {
// public:
//   typedef typename K::Weighted_point_2                Weighted_point_2;
//   typedef typename K::Line_2                          Line_2;

//   typedef Line_2           result_type;

//   Line_2
//   operator() ( const Weighted_point_2 & p, const Weighted_point_2 & q) const
//   {
//     return CGAL::radical_axis(p,q);
//   }
// };

//-----------------------------------------------------------
template < typename K >
class Orientation_on_sphere_3
{
public:
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Comparison_result        Comparison_result;

  typedef Comparison_result   result_type;

  Comparison_result operator()(const Point_2& p,
			       const Point_2& q,
			       const Point_2& t) const
  {
		return orientation(p,q,Point_2(0,0,0),t);	
  }
};

template < typename K >
class In_cone_3
{
public:
  typedef typename K::Point_2                Point_2;
  typedef typename K::Orientation            Orientation;
  typedef typename K::FT                     FT;

  typedef Orientation   result_type;

  result_type
  operator()( const Point_2& p, const Point_2& q,
        const Point_2& r, const Point_2& s) const
  {
		const FT px = p.x();
		const FT py = p.y();
		const FT pz = p.z();
		const FT qx = q.x();
		const FT qy = q.y();
		const FT qz = q.z();
		const FT rx = r.x();
		const FT ry = r.y();
		const FT rz = r.z();
		const FT sx = s.x();
		const FT sy = s.y();
		const FT sz = s.z();
		return in_coneC3(
       px, py, pz,
		   qx, qy, qz,
		   rx, ry, rz,
		   sx, sy, sz
		);		
  }

  Oriented_side operator() (const Point_2& p,
			    const Point_2& q,
			    const Point_2& r) const
  {
	  static Point_2 _sphere = Point_2(0,0,0);
    return -coplanar_orientation(p,q,_sphere,r);
  }

   Oriented_side operator() (const Point_2& p,
			     const Point_2& q) const
  {
	  static Point_2 _sphere = Point_2(0,0,0);
		Comparison_result pq=compare_xyz(p,q);
      
    if(pq==EQUAL){
      return ON_ORIENTED_BOUNDARY;
    }
    Comparison_result sq=compare_xyz(_sphere,q);
    if(pq==sq){
      return ON_POSITIVE_SIDE;
    }
    return ON_NEGATIVE_SIDE;
  }


};


template < typename K >
class Compare_power_distance_sphere_3
{
public:
    typedef typename K::Point_3            Point_3;
    typedef typename K::FT                 FT;
  public:
    typedef typename K::Comparison_result  result_type;

    // the projected distance is invariant to radius
    // and the distance itself can be replace by the inverse of the dot product
    // to compare
    // supposant (0,0,0)
    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& s) const
    {
	    // TODO, precondition: pq, ps, qs < pi and center != p, q or s (p, q, s != (0,0,0))
      return cmp_dist_to_point_projected_on_sphereC3(
        p.x(), p.y(), p.z(),
        q.x(), q.y(), q.z(),
        s.x(), s.y(), s.z()
      );
    }
};

// template < class Bare_point, class Weight >
// inline
// typename Kernel_traits<Bare_point>::Kernel::Comparison_result
// compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
// 		       const Weighted_point<Bare_point, Weight>& q,
// 		       const Bare_point& r, Cartesian_tag)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::FT  FT;
//   return compare_power_distanceC2(p.x(), p.y(), FT(p.weight()),
// 				  q.x(), q.y(), FT(q.weight()),
// 				  r.x(), r.y());
// }

// template < class Bare_point, class Weight >
// inline
// typename Kernel_traits<Bare_point>::Kernel::Comparison_result
// compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
// 		       const Weighted_point<Bare_point, Weight>& q,
// 		       const Bare_point& r, Homogeneous_tag)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::RT  RT;
//   return compare_power_distanceH2(p.hx(), p.hy(), p.hw(), FT(p.weight()),
// 				  q.hx(), q.hy(), q.hw(), FT(q.weight()),
// 				  r.hx(), r.hy(), r.hw());
// }

// template < class Bare_point, class Weight >
// inline
// typename Kernel_traits<Bare_point>::Kernel::Comparison_result
// compare_power_distance(const Weighted_point<Bare_point, Weight>& p,
// 		       const Weighted_point<Bare_point, Weight>& q,
// 		       const Bare_point& r)
// {
//   typedef typename Kernel_traits<Bare_point>::Kernel::Rep_tag Tag;
//   return compare_power_distance(p, q, r, Tag());
// }

// template < typename K >
// class Compare_power_distance_2
// {
// public:
//   typedef typename K::Weighted_point_2         Weighted_point_2;
//   typedef typename K::Point_2                  Point_2;
//   typedef typename K::Comparison_result        Comparison_result;

//   typedef Comparison_result   result_type;

//   Comparison_result operator()(const Point_2& p,
// 			       const Weighted_point_2& q,
// 			       const Weighted_point_2& r) const
//   {
//     return CGAL::compare_power_distance(q, r, p);
//   }
// };

template < class R, class W = typename R::RT>
class Delaunay_triangulation_sphere_traits_base_2
  : public R
{
public:
  typedef R                              Kernel;
  typedef R                              Rep;
  typedef W                              Weight;
  typedef R                              Traits;
  typedef typename Traits::Point_3       Linear_point_3;
  typedef typename Traits::Point_3       Point_2;
  typedef typename Traits::Point_3       Bare_point_2;
  typedef typename Traits::Point_3       Weighted_point_2;
  // This is required for the point() function of vertex base class
  // to be correctly return a weighted_point;
  // patch 27/11/00
  // 18/03/03 I put now the same typedef in Regulat_triangulation_2
  // for the need of hierarchy
  // don't know if this is definitive
  //typedef Weighted_point                        Point_2;

  typedef Delaunay_triangulation_sphere_traits_base_2<R, W>   Self;

  typedef CGAL::In_cone_3<Self>                       Power_test_2;

  typedef CGAL::Orientation_sphere_2<Self>            Orientation_2;
  typedef CGAL::Compare_power_distance_sphere_3<Self> Compare_power_distance_2;

  // bidouille pour utiliser la meme traits pour Delaunay et Regular
  typedef CGAL::In_cone_3<Self>                       Side_of_oriented_circle_2;
  typedef CGAL::Compare_power_distance_sphere_3<Self> Compare_distance_2;

  typedef CGAL::In_cone_3<Self>                       In_cone_3;

  typedef CGAL::Coradial_sphere_2<Self>       Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Self>           Inside_cone_2;
  typedef CGAL::Orientation_sphere_1<Self>    Orientation_1;

  // construction objects
  typedef Dummy_constructor<Self> Construct_circumcenter_3;
  //typedef Construct_weighted_circumcenter_de_Pedro_2<Self>
  //                                          Construct_weighted_circumcenter_2;
  //typedef Construct_radical_axis_de_Pedro_2<Self>  Construct_radical_axis_2;

  // Compare_x / Less_x etc. not needed!
  
  Orientation_2
  orientation_2_object() const {
		static Point_2 p = Point_2(0,0,0);
    return Orientation_2(p);
  }

  Power_test_2 
  power_test_2_object() const
    {  return Power_test_2();}

  Compare_power_distance_2
  compare_power_distance_2_object() const {
    return Compare_power_distance_2();
  }

  Compare_distance_2
  compare_distance_2_object() const {
    return Compare_distance_2();
  }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2();
  }

  In_cone_3
  in_cone_3_object() const {
    return In_cone_3();
  }

  //TODO: uncomment once Pedro's types have been inserted.
  //constructions for dual:
  Construct_circumcenter_3
  construct_circumcenter_3_object() const
    {return Construct_circumcenter_3();}

  // new predicates from Olivier Rouiller
  Orientation_1
  orientation_1_object() const {
		static Point_2 p = Point_2(0,0,0);
    return Orientation_1(p);
  }

  Coradial_sphere_2
  coradial_sphere_2_object() const
  {
	  static Point_2 p = Point_2(0,0,0);
	  return Coradial_sphere_2(p);
	}

  Inside_cone_2
  inside_cone_2_object() const 
  {
	 	static Point_2 p = Point_2(0,0,0);
    return Inside_cone_2(p);
  }


  //Construct_weighted_circumcenter_2
  //construct_weighted_circumcenter_2_object() const
  //  {return Construct_weighted_circumcenter_2();}
  
  //Construct_radical_axis_2
  //construct_radical_axis_2_object() const
  //  {return Construct_radical_axis_2();}
};
 
// We use a base class here to have the specialization below to work.
// Otherwise there is a cycle in the derivation.
template < class R, class W = typename R::RT>
class Delaunay_triangulation_sphere_traits_2
  : public Delaunay_triangulation_sphere_traits_base_2<R, W>
{};

} //namespace CGAL

// Now specialize for Filtered_kernel<CK>, to get
// the filtered traits automatically.

#include <CGAL/Delaunay_triangulation_sphere_filtered_traits_2.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/static_in_cone_ntC3.h>

namespace CGAL {

// This declaration is needed to break the cyclic dependency.
 template < typename K >
 class Delaunay_triangulation_sphere_filtered_traits_2;

#ifndef CGAL_NO_STATIC_FILTERS

// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K >
class Delaunay_triangulation_sphere_static_filtered_traits_2
  : public Delaunay_triangulation_sphere_filtered_traits_2<K>
{

public:

  typedef K               Kernel;

  typedef CGAL::SF_In_cone_3< Delaunay_triangulation_sphere_filtered_traits_2<K> > Power_test_2;
  typedef CGAL::SF_In_cone_3< Delaunay_triangulation_sphere_filtered_traits_2<K> > Side_of_oriented_circle_2;
  typedef CGAL::SF_In_cone_3< Delaunay_triangulation_sphere_filtered_traits_2<K> > In_cone_3;

  Power_test_2 power_test_2_object() const
	{ return Power_test_2();}

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2();
  }

  In_cone_3
  in_cone_3_object() const {
    return In_cone_3();
  }

  // The following are inherited since they are constructions :
  // Construct_weighted_circumcenter_2
  // Construct_radical_axis_2
};

 template < typename K >
 class Delaunay_triangulation_sphere_static_filtered_traits_2;

 template < typename CK >
 class Delaunay_triangulation_sphere_traits_2 < Filtered_kernel<CK> >
   : public Delaunay_triangulation_sphere_static_filtered_traits_2 < Filtered_kernel<CK> >
 {
 public:
   typedef Filtered_kernel<CK>   Kernel;
 };
#else

 template < typename CK >
 class Delaunay_triangulation_sphere_traits_2 < Filtered_kernel<CK> >
   : public Delaunay_triangulation_sphere_filtered_traits_2 < Filtered_kernel<CK> >
 {
 public:
   typedef Filtered_kernel<CK>   Kernel;
 };

#endif

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
