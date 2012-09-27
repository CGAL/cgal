
#ifndef CGAL_REGULAR_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_SPHERE_TRAITS_2_H

//#include <CGAL/Weighted_point.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Triangulation_sphere_traits_2.h>

namespace CGAL { 

template <typename K >
class Power_test_2
{
 public:
  typedef typename K::Point_2   Point_2;
  typedef typename K::Oriented_side     Oriented_side;
  typedef typename K::Comparison_result Comparison_result;

  Power_test_2(const Point_2& sphere);

  Oriented_side operator() (const Point_2& p,
			    const Point_2& q,
			    const Point_2& r,
			    const Point_2& s) const
  {
    return orientation(p,q,r,s);
  }

  Oriented_side operator() (const Point_2& p,
			    const Point_2& q,
			    const Point_2& r) const
  {
    return -coplanar_orientation(p,q,_sphere,r);
  }

   Oriented_side operator() (const Point_2& p,
			     const Point_2& q) const
  {
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

 protected:
   Point_2 _sphere;
};

template < typename K >
Power_test_2<K>::
Power_test_2(const Point_2& sphere)
: _sphere(sphere)
{}


template < class R >
class Regular_triangulation_sphere_traits_2
  : public R
{
public:
  typedef Triangulation_sphere_traits_2<R>                     Base;
  typedef typename R::Point_3                               Point_2; 
  typedef typename R::Point_3                      Weighted_point_2;
                      

  typedef Regular_triangulation_sphere_traits_2<R>   Self;
  typedef CGAL::Power_test_2<Self>            Power_test_2;
  typedef CGAL::Orientation_sphere_2<Self>    Orientation_2;
  typedef CGAL::Coradial_sphere_2<Self>       Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Self>           Inside_cone_2;
  typedef CGAL::Orientation_sphere_1<Self>    Orientation_1;

  Regular_triangulation_sphere_traits_2(const Point_2& sphere=Point_2(0,0,0));

  Orientation_2
  orientation_2_object()const
  {return Orientation_2(_sphere);}

  Orientation_1
  orientation_1_object() const {
    return Orientation_1(_sphere);
  }

  Power_test_2 
  power_test_2_object() const
    {  return Power_test_2(_sphere);}

  Coradial_sphere_2
  coradial_sphere_2_object() const
  {return Coradial_sphere_2(_sphere);}

  Inside_cone_2
  inside_cone_2_object() const {
    return Inside_cone_2(_sphere);
  }

protected :
  Point_2 _sphere;

};

template < class R >
Regular_triangulation_sphere_traits_2<R> ::
Regular_triangulation_sphere_traits_2(const Point_2& sphere)
: _sphere(sphere)
{}
 


} //namespace CGAL

#endif // CGAL_Reg_TRIANGULATION_SPHERE_TRAITS_2_H
