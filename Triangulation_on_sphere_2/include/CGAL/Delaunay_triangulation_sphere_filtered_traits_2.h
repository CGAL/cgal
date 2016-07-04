// Author(s)     : 

#ifndef CGAL_DELAUNAY_TRIANGULATION_SPHERE_FILTERED_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_SPHERE_FILTERED_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Filtered_predicate.h>
//#include <CGAL/static_in_cone_ntC3.h>

namespace CGAL { 

// The Weighted_converter is parametrized by a usual kernel converter,
// and adds the conversions for the Weighted_point.
template < typename Converter >
struct Delaunay_weighted_converter_2
  : Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Delaunay_triangulation_sphere_traits_2<Source_kernel>
                   ::Weighted_point_2  Source_wp;

  typedef typename Delaunay_triangulation_sphere_traits_2<Target_kernel>
                   ::Weighted_point_2  Target_wp;

  /*
  We dont have weight, so our converter wont converter a weighted point.
  We even not need such converter roughly speaking
  Target_wp
  operator()(const Source_wp &wp) const
  {
    return Target_wp(Converter::operator()(wp.point()),
                     Converter::operator()(wp.weight()));
  }*/
	Target_wp
	operator()(const Source_wp &wp) const
	{
	  return Converter::operator()(wp);
	}
};

	/*
	
// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K >
class Delaunay_triangulation_sphere_filtered_traits_2
  : public Delaunay_triangulation_sphere_traits_base_2<K>
{
  // Exact traits is based on the exact kernel.
  typedef Delaunay_triangulation_sphere_traits_2<typename K::Exact_kernel>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Delaunay_triangulation_sphere_traits_2<typename K::Approximate_kernel>
                                                   Filtering_traits;

  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:

  typedef K               Kernel;

  typedef Filtered_predicate<
            typename Exact_traits::Power_test_2,
            typename Filtering_traits::Power_test_2,
            Delaunay_weighted_converter_2<C2E>,
            Delaunay_weighted_converter_2<C2F> >  Power_test_2;

  typedef Filtered_predicate<
    typename Exact_traits::Compare_power_distance_2,
    typename Filtering_traits::Compare_power_distance_2,
    Delaunay_weighted_converter_2<C2E>,
    Delaunay_weighted_converter_2<C2F> >  Compare_power_distance_2;


  // bidouille pour marcher avec DT et RT au meme temp (la traits)
  // je pense que c'est mieu d'en avoir une pour RT et une pour DT

  typedef Filtered_predicate<
            typename Exact_traits::Side_of_oriented_circle_2,
            typename Filtering_traits::Side_of_oriented_circle_2,
            Delaunay_weighted_converter_2<C2E>,
            Delaunay_weighted_converter_2<C2F> >  Side_of_oriented_circle_2;

  typedef Filtered_predicate<
    typename Exact_traits::Compare_distance_2,
    typename Filtering_traits::Compare_distance_2,
    Delaunay_weighted_converter_2<C2E>,
    Delaunay_weighted_converter_2<C2F> >  Compare_distance_2;

  // pour faire jolie pour le static

  typedef Filtered_predicate<
            typename Exact_traits::In_cone_3,
            typename Filtering_traits::In_cone_3,
            Delaunay_weighted_converter_2<C2E>,
            Delaunay_weighted_converter_2<C2F> >  In_cone_3;

  Power_test_2 power_test_2_object() const
	{ return Power_test_2();}

  Compare_power_distance_2 compare_power_distance_2_object() const
  { return Compare_power_distance_2(); }

  // bidouille pour marcher avec DT et RT au meme temp

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


  // The following are inherited since they are constructions :
  // Construct_weighted_circumcenter_2
  // Construct_radical_axis_2
};*/

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_SPHERE_FILTERED_TRAITS_2_H
