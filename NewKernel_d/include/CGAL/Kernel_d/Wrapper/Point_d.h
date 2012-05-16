#ifndef CGAL_WRAPPER_POINT_D_H
#define CGAL_WRAPPER_POINT_D_H

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

namespace CGAL {

template <class R_>
class Point_d : public R_::Kernel_base::Point
		// Deriving won't work if the point is just a __m256d.
{
  typedef typename R_::RT                    RT_;
  typedef typename R_::FT                    FT_;
  typedef typename R_::Kernel_base           Kbase;
  typedef typename R_::Vector                Vector_;
  typedef typename Kbase::template Functor<Construct_ttag<Point_tag> >::type CPBase;
  typedef typename Kbase::template Functor<Compute_point_cartesian_coordinate_tag>::type CCBase;

  typedef Point_d                            Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename R_::Point>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Point_cartesian_const_iterator Cartesian_const_iterator;
  typedef typename Kbase::Point      Rep;

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
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Point_d> >::value>::type> explicit Point_d(U&&...u)
	  : Rep(CPBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//	  : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Point_d(Eval_functor&&,F&&f,U&&...u)
	  : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Point_d(Rep const& v) : Rep(v) {}
  Point_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Point_d(Rep&& v) : Rep(std::move(v)) {}

  // this one should be implicit
  Point_d(Origin const& v)
    : Rep(CPBase()(v)) {}
  Point_d(Origin& v)
    : Rep(CPBase()(v)) {}
  Point_d(Origin&& v)
    : Rep(CPBase()(std::move(v))) {}

#else

  Point_d() : Rep(CPBase()()) {}

  Point_d(Rep const& v) : Rep(v) {} // try not to use it

#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Point_d(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(CPBase()( \
	BOOST_PP_ENUM_PARAMS(N,t))) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Point_d(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(f(BOOST_PP_ENUM_PARAMS(N,t))) {}
  /*
  template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  Point_d(Eval_functor,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(Eval_functor(), BOOST_PP_ENUM_PARAMS(N,t)) {}
  */

  BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
 
  // this one should be implicit
  Point_d(Origin const& o)
    : Rep(CPBase()(o)) {}

#endif

  typename boost::result_of<CCBase(Rep,int)>::type cartesian(int i)const{
	  return CCBase()(rep(),i);
  }


  /*
  Direction_d direction() const
  {
    return R().construct_direction_d_object()(*this);
  }

  Vector_d transform(const Aff_transformation_d &t) const
  {
    return t.transform(*this);
  }

  Vector_d operator/(const RT& c) const
  {
   return R().construct_divided_vector_d_object()(*this,c);
  }

  Vector_d operator/(const typename First_if_different<FT_,RT>::Type & c) const
  {
   return R().construct_divided_vector_d_object()(*this,c);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  x() const
  {
    return R().compute_x_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_y_3, Vector_3>::type
  y() const
  {
    return R().compute_y_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_z_3, Vector_3>::type
  z() const
  {
    return R().compute_z_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hx_3, Vector_3>::type
  hx() const
  {
    return R().compute_hx_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hy_3, Vector_3>::type
  hy() const
  {
    return R().compute_hy_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hz_3, Vector_3>::type
  hz() const
  {
    return R().compute_hz_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
  hw() const
  {
    return R().compute_hw_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) || (i == 2) );
    if (i==0) return x();
    if (i==1) return y();
    return z();
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 3) );
    if (i==0) return hx();
    if (i==1) return hy();
    if (i==2) return hz();
    return hw();
  }

  int dimension() const // bad idea?
  {
      return rep.dimension();
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(*this);
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(*this,3);
  }

  typename Qualified_result_of<typename R::Compute_squared_length_3, Vector_3>::type
  squared_length() const
  {
    return R().compute_squared_length_3_object()(*this);
  }
*/
};
#if 0
template <class R_> Point_d<R_>::Point_d(Point_d &)=default;
#endif

//TODO: IO

//template <class R_>
//Vector_d<R_> operator+(const Vector_d<R_>& v,const Vector_d<R_>& w) const
//{
//	return typename R::template Construct<Sum_of_vectors_tag>::type()(v,w);
//}
//
//template <class R_>
//Vector_d<R_> operator-(const Vector_d<R_>& v,const Vector_d<R_>& w) const
//{
//	return typename R::template Construct<Difference_of_vectors_tag>::type()(v,w);
//}

} //namespace CGAL

#endif // CGAL_WRAPPER_POINT_D_H
