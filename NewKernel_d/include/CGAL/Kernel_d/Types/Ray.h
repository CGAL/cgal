#ifndef CGAL_KERNELD_TYPES_RAY_H
#define CGAL_KERNELD_TYPES_RAY_H
#include <utility>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
namespace CGAL {
template <class R_> class Ray {
	typedef typename R_::FT FT_;
	typedef typename R_::Point Point_;
	typedef typename R_::Vector Vector_;
	typedef std::pair<Point_,Vector_> Data_;
	Data_ data;
	public:
	Ray(){}
	Ray(Point_ const&a, Vector_ const&b): data(a,b) {}
	Point_ source()const{
	  return data.first;
	}
	// FIXME: return a R_::Direction?
	Vector_ direction()const{
	  return data.second;
	}
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_ray : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_ray)
    typedef typename R_::Ray result_type;
    typedef typename R_::Point Point;
    typedef typename R_::Vector Vector;
    typedef typename Get_functor<R_, Difference_of_points_tag>::type Dp_;
    //typedef typename Get_functor<R_, Translated_point_tag>::type Tp_;
    //typedef typename Get_functor<R_, Scaled_vector_tag>::type Sv_;
    result_type operator()(Point const&a, Vector const&b)const{
      return result_type(a,b);
    }
    result_type operator()(Point const&a, typename First_if_different<Point,Vector>::Type const&b)const{
      Dp_ dp(this->kernel());
      return result_type(a,dp(b,a));
    }
  };
}

} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_RAY_H
