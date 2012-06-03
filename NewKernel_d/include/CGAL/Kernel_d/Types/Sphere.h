#ifndef CGAL_KD_TYPE_SPHERE_H
#define CGAL_KD_TYPE_SPHERE_H
#include <CGAL/store_kernel.h>
namespace CGAL {
template <class R_> class Sphere {
	typedef typename R_::FT FT_;
	typedef typename R_::Point Point_;
	Point_ c_;
	FT_ r2_;

	public:
	Sphere(Point_ const&p, FT_ const&r2): c_(p), r2_(r2) {}
	// TODO: Add a piecewise constructor?

	Point_ center()const{return c_;}
	FT_ squared_radius()const{return r2_;}
};
namespace CartesianDKernelFunctors {
template <class R_> struct Construct_sphere : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_sphere)
  typedef typename R_::Sphere result_type;
  typedef typename R_::Point Point;
  typedef typename R_::FT FT;
  typedef typename R_::LA LA;
  result_type operator()(Point const&a, FT const&b)const{
    return result_type(a,b);
  }
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    throw "not implemented yet!";
  }
};
template <class R_> struct Center_of_sphere {
  CGAL_FUNCTOR_INIT_IGNORE(Center_of_sphere)
  typedef typename R_::Sphere Sphere;
  typedef typename R_::Point result_type;
  result_type operator()(Sphere const&s)const{
    return s.center();
  }
};
template <class R_> struct Squared_radius {
  CGAL_FUNCTOR_INIT_IGNORE(Squared_radius)
  typedef typename R_::Sphere Sphere;
  typedef typename R_::FT result_type;
  // TODO: Is_exact?
  result_type operator()(Sphere const&s)const{
    return s.squared_radius();
  }
};
}
} // namespace CGAL
#endif
