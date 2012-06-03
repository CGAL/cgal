#ifndef CGAL_KD_TYPE_SPHERE_H
#define CGAL_KD_TYPE_SPHERE_H
#include <CGAL/store_kernel.h>
namespace CGAL {
template <class R_> class Hyperplane {
	typedef typename R_::FT FT_;
	typedef typename R_::Point Point_;
	typedef typename R_::Vector Vector_;
	Vector_ v_;
	FT_ s_;

	public:
	Hyperplane(Vector_ const&v, FT_ const&s): v_(v), s_(s) {}
	// TODO: Add a piecewise constructor?

	Vector_ orthogonal_vector()const{return v_;}
	FT_ translation()const{return s_;}
};
namespace CartesianDKernelFunctors {
template <class R_> struct Construct_hyperplane : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_hyperplane)
  typedef typename R_::Hyperplane result_type;
  typedef typename R_::Point Point;
  typedef typename R_::Vector Vector;
  typedef typename R_::FT FT;
  typedef typename R_::LA LA;
  result_type operator()(Vector const&a, FT const&b)const{
    return result_type(a,b);
  }
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    throw "not implemented yet!";
  }
};
template <class R_> struct Orthogonal_vector {
  CGAL_FUNCTOR_INIT_IGNORE(Orthogonal_vector)
  typedef typename R_::Hyperplane Hyperplane;
  typedef typename R_::Point result_type;
  result_type operator()(Hyperplane const&s)const{
    return s.orthogonal_vector();
  }
};
template <class R_> struct Hyperplane_translation {
  CGAL_FUNCTOR_INIT_IGNORE(Hyperplane_translation)
  typedef typename R_::Hyperplane Hyperplane;
  typedef typename R_::FT result_type;
  // TODO: Is_exact?
  result_type operator()(Hyperplane const&s)const{
    return s.translation();
  }
};
}
} // namespace CGAL
#endif
