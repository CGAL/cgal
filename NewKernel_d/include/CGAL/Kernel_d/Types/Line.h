#ifndef CGAL_KERNELD_TYPES_LINE_H
#define CGAL_KERNELD_TYPES_LINE_H
#include <utility>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
namespace CGAL {
template <class R_> class Line {
	typedef typename R_::FT FT_;
	typedef typename R_::Point Point_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	Line(){}
	Line(Point_ const&a, Point_ const&b): data(a,b) {}
	Point_ point(int i)const{
	  if(i==0) return data.first;
	  if(i==1) return data.second;
	  throw "not implemented";
	}
	Line opposite()const{
		return Line(data.second,data.first);
	}
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_line : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_line)
    typedef typename R_::Line result_type;
    typedef typename R_::Point Point;
    typedef typename R_::Vector Vector;
    typedef typename Get_functor<R_, Translated_point_tag>::type Tp_;
    //typedef typename Get_functor<R_, Difference_of_points_tag>::type Dp_;
    //typedef typename Get_functor<R_, Scaled_vector_tag>::type Sv_;
    result_type operator()(Point const&a, Point const&b)const{
      return result_type(a,b);
    }
    result_type operator()(Point const&a, typename First_if_different<Vector,Point>::Type const&b)const{
      Tp_ tp(this->kernel());
      return result_type(a,tp(a,b));
    }
  };
}

} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_LINE_H
