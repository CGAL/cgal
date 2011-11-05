#ifndef CGAL_KD_KO_CONVERTER_H
#define CGAL_KD_KO_CONVERTER_H
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h> // First_if_different
namespace CGAL {
template <class Tag_, class K1, class K2> struct KO_converter;
//TODO: It would probably be better if this was a Misc Functor in K1.

template <class K1, class K2> struct KO_converter<Point_tag,K1,K2>{
	typedef typename K1::template Type<Point_tag>::type argument_type;
	typedef typename K2::template Type<Point_tag>::type result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& p) const {
		typename K1::template Functor<Construct_point_cartesian_const_iterator_tag>::type i(k1);
		typename K2::template Functor<Construct_ttag<Point_tag> >::type cp(k2);
		return cp(conv(i(p,Begin_tag())),conv(i(p,End_tag())));
	}
};

template <class K1, class K2> struct KO_converter<Vector_tag,K1,K2>{
	typedef typename K1::template Type<Vector_tag>::type K1_Vector;
	typedef typename K1::template Type< Point_tag>::type K1_Point;
	// can't use vector without at least a placeholder point because of this
	typedef typename First_if_different<K1_Vector,K1_Point>::Type argument_type;
	typedef typename K2::template Type<Vector_tag>::type result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& v) const {
		typename K1::template Functor<Construct_vector_cartesian_const_iterator_tag>::type i(k1);
		typename K2::template Functor<Construct_ttag<Vector_tag> >::type cp(k2);
		return cp(conv(i(v,Begin_tag())),conv(i(v,End_tag())));
	}
};

template <class K1, class K2> struct KO_converter<Segment_tag,K1,K2>{
	typedef typename K1::template Type<Segment_tag>::type argument_type;
	typedef typename K2::template Type<Segment_tag>::type result_type;
	template <class C>
	result_type operator()(K1 const& k1, K2 const& k2, C const& conv, argument_type const& s) const {
		typename K1::template Functor<Construct_segment_extremity_tag>::type f(k1);
		typename K2::template Functor<Construct_ttag<Segment_tag> >::type cs(k2);
		return cs(conv(f(s,0)),conv(f(s,1)));
	}
};

}
#endif
