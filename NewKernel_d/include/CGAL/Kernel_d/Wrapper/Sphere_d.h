#ifndef CGAL_WRAPPER_SPHERE_D_H
#define CGAL_WRAPPER_SPHERE_D_H

#include <CGAL/representation_tags.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#ifndef CGAL_CXX0X
#include <boost/preprocessor/repetition.hpp>
#endif
#include <boost/utility/result_of.hpp>

namespace CGAL {

template <class R_>
class Sphere_d : public R_::Kernel_base::Sphere
{
  typedef typename R_::FT                    FT_;
  typedef typename R_::Kernel_base           Kbase;
  typedef typename R_::Point                 Point_;
  typedef typename Kbase::template Functor<Construct_ttag<Sphere_tag> >::type CSBase;
  typedef typename Kbase::template Functor<Center_of_sphere_tag>::type COSBase;
  typedef typename Kbase::template Functor<Squared_radius_tag>::type SRBase;

  typedef Sphere_d                            Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename R_::Sphere>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef typename Increment_dimension<Ambient_dimension,-1>::type Feature_dimension;

  typedef typename Kbase::Sphere      Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

#ifdef CGAL_CXX0X
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Sphere_d> >::value>::type> explicit Sphere_d(U&&...u)
	  : Rep(CSBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//	  : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Sphere_d(Eval_functor&&,F&&f,U&&...u)
	  : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Sphere_d(Rep const& v) : Rep(v) {}
  Sphere_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Sphere_d(Rep&& v) : Rep(std::move(v)) {}

#else

  Sphere_d() : Rep(CSBase()()) {}

  Sphere_d(Rep const& v) : Rep(v) {} // try not to use it

#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Sphere_d(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(CSBase()( \
	BOOST_PP_ENUM_PARAMS(N,t))) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Sphere_d(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(f(BOOST_PP_ENUM_PARAMS(N,t))) {}
  /*
  template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  Point_d(Eval_functor,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(Eval_functor(), BOOST_PP_ENUM_PARAMS(N,t)) {}
  */

  BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
 
#endif

    //TODO: if COSBase returns a reference to a base point, cast it to a
    //reference to a wrapper point. Ugly but should be safe.
    Point_ center()const{
      return Point_(Eval_functor(),COSBase(),rep());
    }
  FT_ squared_radius()const{
    return SRBase()(rep());
  }

};

} //namespace CGAL

#endif // CGAL_WRAPPER_SPHERE_D_H
