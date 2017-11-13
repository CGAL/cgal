// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_CARTESIAN_LA_BASE_H
#define CGAL_KERNEL_D_CARTESIAN_LA_BASE_H

#include <CGAL/basic.h>
#include <CGAL/Origin.h>
#include <boost/type_traits/integral_constant.hpp>
#include <CGAL/representation_tags.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/typeset.h>
#include <CGAL/NewKernel_d/Dimension_base.h>
#include <CGAL/NewKernel_d/Cartesian_LA_functors.h>
#include <CGAL/NewKernel_d/Vector/array.h>
#include <CGAL/NewKernel_d/Vector/vector.h>
#include <CGAL/NewKernel_d/Vector/mix.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/NewKernel_d/LA_eigen/LA.h>
#else
#error Eigen3 is required
#endif

namespace CGAL {

template < typename FT_, typename Dim_,
#if 1
	 typename Vec_=Mix_vector<Array_vector<FT_, Dim_>,
				  Vector_vector<FT_, Dim_>,
				  FT_, Dim_>,
#elif 0
	 typename Vec_=Array_vector<FT_, Dim_>,
#elif 0
	 typename Vec_=Vector_vector<FT_, Dim_>,
#else
	 // Dangerous because of alignment. Ok on x86_64 without AVX.
	 typename Vec_=LA_eigen<FT_, Dim_>,
#endif
	 typename LA_=LA_eigen<FT_,Dim_> >
  /* Default LA to Vec or to LA_eigen? */
struct Cartesian_LA_base_d : public Dimension_base<Dim_>
{
    typedef Cartesian_LA_base_d<FT_,Dim_>               Self;
    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;
    typedef Dim_              Default_ambient_dimension;
    typedef Dim_              Max_ambient_dimension;
    typedef Dim_              Dimension;
    typedef LA_               LA;
    template <class> struct Ambient_dimension { typedef Dim_ type; };

    typedef Vec_     LA_vector;
    typedef typename LA_vector::Vector Point;
    typedef typename LA_vector::Vector Vector;
    typedef typename LA_vector::Vector Vector_;
    typedef typename LA_vector::Construct_vector Constructor;
    typedef typename LA_vector::Vector_const_iterator Point_cartesian_const_iterator;
    typedef typename LA_vector::Vector_const_iterator Vector_cartesian_const_iterator;

    template<class, class=void> struct Type {};
    template<class D> struct Type< Point_tag, D> { typedef Vector_ type; };
    template<class D> struct Type<Vector_tag, D> { typedef Vector_ type; };
    template<class D> struct Type<    FT_tag, D> { typedef     FT_ type; };
    template<class D> struct Type<    RT_tag, D> { typedef     FT_ type; };

    typedef typeset<Point_tag>
      ::add<Vector_tag>::type
    // FIXME: These have nothing to do here.
      ::add<Segment_tag>::type
      ::add<Hyperplane_tag>::type
      ::add<Sphere_tag>::type
      ::add<Weighted_point_tag>::type
      Object_list;

    typedef typeset< Point_cartesian_const_iterator_tag>::type
      ::add<Vector_cartesian_const_iterator_tag>::type
      Iterator_list;

    template<class, class=void, class=boost::integral_constant<int,0> > struct Functor {
	    typedef Null_functor type;
    };
    template<class D> struct Functor<Construct_ttag<Vector_tag>,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Null_vector> type;
    };
    template<class D> struct Functor<Construct_ttag<Point_tag>,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Origin> type;
    };
    template<class D> struct Functor<Construct_ttag<Point_cartesian_const_iterator_tag>,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
    template<class D> struct Functor<Construct_ttag<Vector_cartesian_const_iterator_tag>,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
    template<class D> struct Functor<Sum_of_vectors_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_vector_plus_minus_tag>::value> > {
	    typedef CartesianDVectorBase::Sum_of_vectors<Self> type;
    };
    template<class D> struct Functor<Difference_of_vectors_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_vector_plus_minus_tag>::value> > {
	    typedef CartesianDVectorBase::Difference_of_vectors<Self> type;
    };
    template<class D> struct Functor<Opposite_vector_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_vector_plus_minus_tag>::value> > {
	    typedef CartesianDVectorBase::Opposite_vector<Self> type;
    };
    template<class D> struct Functor<Midpoint_tag,D,
      boost::integral_constant<int,
	   !LA_vector::template Property<Has_vector_plus_minus_tag>::value
	|| !LA_vector::template Property<Has_vector_scalar_ops_tag>::value> > {
	    typedef CartesianDVectorBase::Midpoint<Self> type;
    };
    template<class D> struct Functor<Compute_point_cartesian_coordinate_tag,D> {
	    typedef CartesianDVectorBase::Compute_cartesian_coordinate<Self> type;
    };
    template<class D> struct Functor<Compute_vector_cartesian_coordinate_tag,D> {
	    typedef CartesianDVectorBase::Compute_cartesian_coordinate<Self> type;
    };
    template<class D> struct Functor<Point_dimension_tag,D> {
	    typedef CartesianDVectorBase::PV_dimension<Self> type;
    };
    template<class D> struct Functor<Vector_dimension_tag,D> {
	    typedef CartesianDVectorBase::PV_dimension<Self> type;
    };
    template<class D> struct Functor<Orientation_of_vectors_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_determinant_of_iterator_to_vectors_tag>::value> > {
	    typedef CartesianDVectorBase::Orientation_of_vectors<Self> type;
    };
    template<class D> struct Functor<Orientation_of_points_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_determinant_of_iterator_to_points_tag>::value> > {
	    typedef CartesianDVectorBase::Orientation_of_points<Self> type;
    };
    template<class D> struct Functor<Scalar_product_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Has_dot_product_tag>::value> > {
	    typedef CartesianDVectorBase::Scalar_product<Self> type;
    };
    template<class D> struct Functor<Squared_distance_to_origin_tag,D,
      boost::integral_constant<int,!LA_vector::template Property<Stores_squared_norm_tag>::value> > {
	    typedef CartesianDVectorBase::Squared_distance_to_origin_stored<Self> type;
    };
    // Use integral_constant<int,2> in case of failure, to distinguish from the previous one.
    template<class D> struct Functor<Squared_distance_to_origin_tag,D,
      boost::integral_constant<int,
	(LA_vector::template Property<Stores_squared_norm_tag>::value
	 || !LA_vector::template Property<Has_dot_product_tag>::value)*2> > {
	    typedef CartesianDVectorBase::Squared_distance_to_origin_via_dotprod<Self> type;
    };
    template<class D> struct Functor<Point_to_vector_tag,D> {
	    typedef CartesianDVectorBase::Identity_functor<Self> type;
    };
    template<class D> struct Functor<Vector_to_point_tag,D> {
	    typedef CartesianDVectorBase::Identity_functor<Self> type;
    };

    CGAL_CONSTEXPR Cartesian_LA_base_d(){}
    CGAL_CONSTEXPR Cartesian_LA_base_d(int d):Dimension_base<Dim_>(d){}
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_LA_BASE_H
