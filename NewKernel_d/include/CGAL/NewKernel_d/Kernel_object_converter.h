// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_KO_CONVERTER_H
#define CGAL_KD_KO_CONVERTER_H
#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Kernel/mpl.h> // First_if_different
#include <CGAL/Dimension.h>
namespace CGAL {
template <class Tag_, class K1, class K2> struct KO_converter;
//TODO: It would probably be better if this was a Misc Functor in K1.
// This way K1 could chose how it wants to present its points (sparse
// iterator?) and derived classes would inherit it.

namespace internal {
template <class D /*=Dynamic_dimension_tag*/, class K1, class K2>
struct Point_converter_help {
        typedef typename Get_type<K1, Point_tag>::type        argument_type;
        typedef typename Get_type<K2, Point_tag>::type        result_type;
        template <class C>
        result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
                typename Get_functor<K1, Construct_ttag<Point_cartesian_const_iterator_tag> >::type i(k1);
                typename Get_functor<K2, Construct_ttag<Point_tag> >::type cp(k2);
                return cp(conv(i(p,Begin_tag())),conv(i(p,End_tag())));
        }
};

// This doesn't seem so useful, the compiler should be able to handle
// the iterators just as efficiently.
template <int d, class K1, class K2>
struct Point_converter_help<Dimension_tag<d>,K1,K2> {
        typedef typename Get_type<K1, Point_tag>::type        argument_type;
        typedef typename Get_type<K2, Point_tag>::type        result_type;
        template <class C,std::size_t...I>
        result_type help(std::index_sequence<I...>, K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
                typename Get_functor<K1, Compute_point_cartesian_coordinate_tag>::type cc(k1);
                typename Get_functor<K2, Construct_ttag<Point_tag> >::type cp(k2);
                return cp(conv(cc(p,I))...);
        }
        template <class C>
        result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
                return help(std::make_index_sequence<d>(),k1,k2,conv,p);
        }
};
}
template <class K1, class K2> struct KO_converter<Point_tag,K1,K2>
: internal::Point_converter_help<typename K1::Default_ambient_dimension,K1,K2>
{};

template <class K1, class K2> struct KO_converter<Vector_tag,K1,K2>{
        typedef typename Get_type<K1, Vector_tag>::type        K1_Vector;

        // Disabling is now done in KernelD_converter
        // // can't use vector without at least a placeholder point because of this
        // typedef typename K1:: Point K1_Point;
        // typedef typename First_if_different<K1_Vector,K1_Point>::Type argument_type;

        typedef K1_Vector argument_type;
        typedef typename Get_type<K2, Vector_tag>::type        result_type;
        template <class C>
        result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& v) const {
                typename Get_functor<K1, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type i(k1);
                typename Get_functor<K2, Construct_ttag<Vector_tag> >::type cp(k2);
                return cp(conv(i(v,Begin_tag())),conv(i(v,End_tag())));
        }
};

template <class K1, class K2> struct KO_converter<Segment_tag,K1,K2>{
        typedef typename Get_type<K1, Segment_tag>::type        argument_type;
        typedef typename Get_type<K2, Segment_tag>::type        result_type;
        template <class C>
        result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
                typename Get_functor<K1, Segment_extremity_tag>::type f(k1);
                typename Get_functor<K2, Construct_ttag<Segment_tag> >::type cs(k2);
                return cs(conv(f(s,0)),conv(f(s,1)));
        }
};

template <class K1, class K2> struct KO_converter<Iso_box_tag,K1,K2>{
        typedef typename Get_type<K1, Iso_box_tag>::type        argument_type;
        typedef typename Get_type<K2, Iso_box_tag>::type        result_type;
        template <class C>
        result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
                typename Get_functor<K1, Construct_min_vertex_tag>::type f(k1);
                typename Get_functor<K1, Construct_max_vertex_tag>::type g(k1);
                typename Get_functor<K2, Construct_ttag<Iso_box_tag> >::type cib(k2);
                return cib(conv(f(s)),conv(g(s)));
        }
};

template <class K1, class K2> struct KO_converter<Hyperplane_tag,K1,K2>{
  typedef typename Get_type<K1, Hyperplane_tag>::type        argument_type;
  typedef typename Get_type<K2, Hyperplane_tag>::type        result_type;
  template <class C>
    result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& h) const {
      typename Get_functor<K1, Orthogonal_vector_tag>::type ov(k1);
      typename Get_functor<K1, Hyperplane_translation_tag>::type ht(k1);
      typename Get_functor<K2, Construct_ttag<Hyperplane_tag> >::type ch(k2);
      return ch(conv(ov(h)),conv(ht(h)));
    }
};

template <class K1, class K2> struct KO_converter<Sphere_tag,K1,K2>{
  typedef typename Get_type<K1, Sphere_tag>::type        argument_type;
  typedef typename Get_type<K2, Sphere_tag>::type        result_type;
  template <class C>
    result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
      typename Get_functor<K1, Center_of_sphere_tag>::type cos(k1);
      typename Get_functor<K1, Squared_radius_tag>::type sr(k1);
      typename Get_functor<K2, Construct_ttag<Sphere_tag> >::type cs(k2);
      return cs(conv(cos(s)),conv(sr(s)));
    }
};

template <class K1, class K2> struct KO_converter<Weighted_point_tag,K1,K2>{
  typedef typename Get_type<K1, Weighted_point_tag>::type        argument_type;
  typedef typename Get_type<K2, Weighted_point_tag>::type        result_type;
  template <class C>
    result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
      typename Get_functor<K1, Point_drop_weight_tag>::type pdw(k1);
      typename Get_functor<K1, Point_weight_tag>::type pw(k1);
      typename Get_functor<K2, Construct_ttag<Weighted_point_tag> >::type cwp(k2);
      return cwp(conv(pdw(s)),conv(pw(s)));
    }
};

}
#endif
