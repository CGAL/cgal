#ifndef CGAL_WRAPPER_POINT_RC_D_H
#define CGAL_WRAPPER_POINT_RC_D_H

#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#ifndef CGAL_CXX0X
#include <boost/preprocessor/repetition.hpp>
#endif
#include <boost/utility/result_of.hpp>

// no need for a fancy interface here, people can use the Point_d wrapper on
// top.

namespace CGAL {

template <class R_>
class Point_rc_d
{
  typedef typename R_::Kernel_base           Kbase;
  typedef typename Kbase::template Functor<Construct_ttag<Point_tag> >::type CPBase;
  typedef typename Kbase::template Functor<Compute_cartesian_coordinate_tag>::type CCBase;

  typedef Point_rc_d                            Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename R_::template Type<Point_tag>::type>::value));

public:
  typedef R_ R;

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Point_cartesian_const_iterator Cartesian_const_iterator;
  typedef typename Kbase::template Type<Point_tag>::type      Rep;
  typedef Handle_for<Rep> Data;

private:
  Data data;
public:

  const Rep& rep() const
  {
    return CGAL::get(data);
  }

#ifdef CGAL_CXX0X
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Point_rc_d> >::value>::type> explicit Point_rc_d(U&&...u)
	  : data(Eval_functor(),CPBase(),std::forward<U>(u)...){}

  template<class F,class...U> explicit Point_rc_d(Eval_functor&&,F&&f,U&&...u)
	  : data(Eval_functor(),std::forward<F>(f),std::forward<U>(u)...){}

  // try not to use these
  Point_rc_d(Rep const& v) : data(v) {}
  Point_rc_d(Rep& v) : data(static_cast<Rep const&>(v)) {}
  Point_rc_d(Rep&& v) : data(std::move(v)) {}

  // this one should be implicit
  Point_rc_d(Origin const& v)
    : data(Eval_functor(),CPBase(),v) {}
  Point_rc_d(Origin& v)
    : data(Eval_functor(),CPBase(),v) {}
  Point_rc_d(Origin&& v)
    : data(Eval_functor(),CPBase(),std::move(v)) {}

#else

  Point_rc_d() : data(Eval_functor(),CPBase()) {}

  Point_rc_d(Rep const& v) : data(v) {} // try not to use it

#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Point_rc_d(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : data(Eval_functor(),CPBase(),BOOST_PP_ENUM_PARAMS(N,t)) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Point_rc_d(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : data(Eval_functor(),f,BOOST_PP_ENUM_PARAMS(N,t)) {}

  BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
 
  // this one should be implicit
  Point_rc_d(Origin const& o)
    : data(Eval_functor(),CPBase(),o) {}

#endif

};

} //namespace CGAL

#endif // CGAL_WRAPPER_POINT_RC_D_H
