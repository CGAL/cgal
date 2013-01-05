#ifndef CGAL_KD_DEFINE_SEGMENT_H
#define CGAL_KD_DEFINE_SEGMENT_H
#include <utility>
#include <CGAL/marcutils.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel_d/Types/Segmentd.h>

namespace CGAL {
namespace CartesianDKernelFunctors {

template<class R_> struct Construct_segment {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_segment)
	typedef R_ R;
	typedef typename Get_type<R_, Point_tag>::type	Point;
	typedef typename Get_type<R_, Segment_tag>::type	Segment;
	typedef typename Get_functor<R_, Construct_ttag<Point_tag> >::type CP;
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
template<class R_> struct Segment_extremity {
	CGAL_FUNCTOR_INIT_IGNORE(Segment_extremity)
	typedef R_ R;
	typedef typename Get_type<R_, Point_tag>::type	Point;
	typedef typename Get_type<R_, Segment_tag>::type	Segment;
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

	typedef CGAL::Segment<Derived> Segment;
	typedef typename Base::Object_list::template add<Segment_tag>::type Object_list;
	template<class T,class=void> struct Type : Inherit_type<Base_, T> {};
	template<class D> struct Type<Segment_tag, D> {
	  typedef CGAL::Segment<Derived> type;
	};

	// TODO: forward the second Functor argument (like fast, no_filter)
	template<class T,class=void> struct Functor : Inherit_functor<Base_, T> {};

	template<class D> struct Functor<Construct_ttag<Segment_tag>,D> {
		typedef CartesianDKernelFunctors::Construct_segment<Derived> type;
	};
	template<class D> struct Functor<Segment_extremity_tag,D> {
		typedef CartesianDKernelFunctors::Segment_extremity<Derived> type;
	};
};
}
#endif
