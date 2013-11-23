#ifndef CGAL_KERNELD_SEGMENTD_H
#define CGAL_KERNELD_SEGMENTD_H
#include <utility>
#include <CGAL/functor_tags.h>
namespace CGAL {
template <class R_> class Segment {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	//typedef typename R_::Vector Vector_;
	//typedef typename Get_functor<R_, Construct_ttag<Vector_tag> >::type Cv_;
//	typedef typename R_::Squared_distance Csd_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	//typedef Segmentd<R_> Segment;
#ifdef CGAL_CXX0X
	//FIXME: don't forward directly, piecewise_constuct should call the point construction functor (I guess? or is it unnecessary?)
	template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Segment>>::value>::type>
	Segment(U&&...u):data(std::forward<U>(u)...){}
#else
	Segment(){}
	Segment(Point_ const&a, Point_ const&b): data(a,b) {}
	//template<class A,class T1,class T2>
	  //Segment(A const&,T1 const&t1,T2 const&t2)
#endif
	Point_ source()const{return data.first;}
	Point_ target()const{return data.second;}
	Point_ operator[](int i)const{
		if((i%2)==0)
			return source();
		else
			return target();
	}
	Segment opposite()const{
		return Segment(target(),source());
	}
	//Vector_ vector()const{
	//	return Cv_()(data.first,data.second);
	//}
//	FT_ squared_length()const{
//		return Csd_()(data.first,data.second);
//	}
};

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

CGAL_KD_DEFAULT_TYPE(Segment_tag,(CGAL::Segment<K>),(Point_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Segment_tag>,(CartesianDKernelFunctors::Construct_segment<K>),(Segment_tag,Point_tag),(Construct_ttag<Point_tag>));
CGAL_KD_DEFAULT_FUNCTOR(Segment_extremity_tag,(CartesianDKernelFunctors::Segment_extremity<K>),(Segment_tag,Point_tag),());

} // namespace CGAL

#endif // CGAL_KERNELD_SEGMENTD_H
