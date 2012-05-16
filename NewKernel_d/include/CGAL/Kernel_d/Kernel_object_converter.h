#ifndef CGAL_KD_KO_CONVERTER_H
#define CGAL_KD_KO_CONVERTER_H
#include <CGAL/marcutils.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h> // First_if_different
namespace CGAL {
template <class Tag_, class K1, class K2> struct KO_converter;
//TODO: It would probably be better if this was a Misc Functor in K1.
// This way K1 could chose how it wants to present its points (sparse
// iterator?) and derived classes would inherit it.

namespace internal {
template <class D /*=Dynamic_dimension_tag*/, class K1, class K2>
struct Point_converter_help {
	typedef typename K1::Point argument_type;
	typedef typename K2::Point result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
		typename K1::template Functor<Construct_ttag<Point_cartesian_const_iterator_tag> >::type i(k1);
		typename K2::template Functor<Construct_ttag<Point_tag> >::type cp(k2);
		return cp(conv(i(p,Begin_tag())),conv(i(p,End_tag())));
	}
};
#ifdef CGAL_CXX0X
// This doesn't seem so useful, the compiler should be able to handle
// the iterators just as efficiently.
template <int d, class K1, class K2>
struct Point_converter_help<Dimension_tag<d>,K1,K2> {
	typedef typename K1::Point argument_type;
	typedef typename K2::Point result_type;
	template <class C,int...I>
	result_type help(Indices<I...>, K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
		typename K1::template Functor<Compute_point_cartesian_coordinate_tag>::type cc(k1);
		typename K2::template Functor<Construct_ttag<Point_tag> >::type cp(k2);
		return cp(conv(cc(p,I))...);
	}
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
		return help(typename N_increasing_indices<d>::type(),k1,k2,conv,p);
	}
};
#endif
}
template <class K1, class K2> struct KO_converter<Point_tag,K1,K2>
: internal::Point_converter_help<typename K1::Default_ambient_dimension,K1,K2>
{};

template <class K1, class K2> struct KO_converter<Vector_tag,K1,K2>{
	typedef typename K1::Vector K1_Vector;
	
	// Disabling is now done in KernelD_converter
	// // can't use vector without at least a placeholder point because of this
	// typedef typename K1:: Point K1_Point;
	// typedef typename First_if_different<K1_Vector,K1_Point>::Type argument_type;

	typedef K1_Vector argument_type;
	typedef typename K2::Vector result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& v) const {
		typename K1::template Functor<Construct_ttag<Vector_cartesian_const_iterator_tag> >::type i(k1);
		typename K2::template Functor<Construct_ttag<Vector_tag> >::type cp(k2);
		return cp(conv(i(v,Begin_tag())),conv(i(v,End_tag())));
	}
};

template <class K1, class K2> struct KO_converter<Segment_tag,K1,K2>{
	typedef typename K1::Segment argument_type;
	typedef typename K2::Segment result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
		typename K1::template Functor<Segment_extremity_tag>::type f(k1);
		typename K2::template Functor<Construct_ttag<Segment_tag> >::type cs(k2);
		return cs(conv(f(s,0)),conv(f(s,1)));
	}
};

}
#endif
