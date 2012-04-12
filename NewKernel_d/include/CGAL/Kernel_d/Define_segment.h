#ifndef CGAL_KD_DEFINE_SEGMENT_H
#define CGAL_KD_DEFINE_SEGMENT_H
#include <utility>
#include <CGAL/marcutils.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel_d/Segmentd.h>

namespace CGAL {
namespace CartesianDKernelFunctors {

template<class R_> struct Construct_segment {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_segment)
	typedef R_ R;
	typedef typename R_::Point Point;
	typedef typename R_::Segment Segment;
	typedef typename R_::template Functor<Construct_ttag<Point_tag> >::type CP;
	typedef Segment result_type;
	result_type operator()(Point const&a, Point const&b)const{
		return result_type(a,b);
	}
	// T should only be std::piecewise_construct_t, but we shouldn't fail if it doesn't exist.
	template<class T,class U,class V>
	result_type operator()(CGAL_FORWARDABLE(T),CGAL_FORWARDABLE(U) u,CGAL_FORWARDABLE(V) v)const{
		CP cp(this->kernel());
		result_type r = {{
			call_on_tuple_elements<Point>(cp, CGAL_FORWARD(U,u)),
			call_on_tuple_elements<Point>(cp, CGAL_FORWARD(V,v)) }};
		return r;
	}
};

// This should be part of Construct_point, according to Kernel_23 conventions
template<class R_> struct Construct_segment_extremity {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_segment_extremity)
	typedef R_ R;
	typedef typename R_::Point Point;
	typedef typename R_::Segment Segment;
	typedef Point result_type;
	result_type operator()(Segment const&s, int i)const{
		if(i==0) return s.source();
		CGAL_assertion(i==1);
		return s.target();
	}
#ifdef CGAL_CXX0X
	result_type operator()(Segment &&s, int i)const{
		if(i==0) return std::move(s.source());
		CGAL_assertion(i==1);
		return std::move(s.target());
	}
#endif
};
} // CartesianDKernelFunctors

template<class Base_, class Derived_=Default>
struct Define_segment : public Base_ {
	typedef Base_ Base;
	typedef Define_segment<Base_,Derived_> Self;
	typedef typename Default::Get<Derived_,Self>::type Derived;

	typedef SegmentCd<Derived> Segment;
	typedef typename Base::Object_list::template add<Segment_tag>::type Object_list;

	// TODO: forward the second Functor argument (like fast, no_filter)
	template<class T,class=void> struct Functor : Base_::template Functor<T> {};

	template<class D> struct Functor<Construct_ttag<Segment_tag>,D> {
		typedef CartesianDKernelFunctors::Construct_segment<Derived> type;
	};
	template<class D> struct Functor<Construct_segment_extremity_tag,D> {
		typedef CartesianDKernelFunctors::Construct_segment_extremity<Derived> type;
	};
};
}
#endif
