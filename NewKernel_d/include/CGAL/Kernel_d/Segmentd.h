#ifndef CGAL_KERNELD_SEGMENTD_H
#define CGAL_KERNELD_SEGMENTD_H
#include <utility>
#include <CGAL/functor_tags.h>
#define Segmentd SegmentCd
namespace CGAL {
template <class R_> class Segmentd {
	typedef typename R_::FT FT_;
	typedef typename R_::template Type<Point_tag>::type Point_;
	typedef typename R_::template Type<Vector_tag>::type Vector_;
	typedef typename R_::template Functor<Construct_ttag<Vector_tag> >::type Cv_;
//	typedef typename R_::Compute_squared_distance Csd_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	typedef Segmentd<R_> Segment;
#ifdef CGAL_CXX0X
	template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Segmentd>>::value>::type>
	Segmentd(U&&...u):data(std::forward<U>(u)...){}
#else
	Segmentd(){}
	Segmentd(Point_ const&a, Point_ const&b): data(a,b) {}
#endif
	Point_ source()const{return data.first;}
	Point_ target()const{return data.second;}
	Point_ operator[](int i)const{
		if((i%2)==0)
			return data.first;
		else
			return data.second;
	}
	Segmentd opposite()const{
		return Segmentd(data.second,data.first);
	}
	Vector_ vector()const{
		return Cv_()(data.first,data.second);
	}
//	FT_ squared_length()const{
//		return Csd_()(data.first,data.second);
//	}
};

} // namespace CGAL

#undef Segmentd
#endif // CGAL_KERNELD_SEGMENTD_H
