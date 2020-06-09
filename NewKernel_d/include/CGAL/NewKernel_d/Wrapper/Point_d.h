// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_WRAPPER_POINT_D_H
#define CGAL_WRAPPER_POINT_D_H
#include <ostream>
#include <istream>
#include <CGAL/IO/io.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Wrap {

template <class R_>
class Point_d : public Get_type<typename R_::Kernel_base, Point_tag>::type
                // Deriving won't work if the point is just a __m256d.
                // Test boost/std::is_class for instance
{
  typedef typename Get_type<R_, RT_tag>::type                RT_;
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename R_::Kernel_base                Kbase;
  typedef typename Get_type<R_, Vector_tag>::type        Vector_;
  typedef typename Get_functor<Kbase, Construct_ttag<Point_tag> >::type CPBase;
  typedef typename Get_functor<Kbase, Compute_point_cartesian_coordinate_tag>::type CCBase;
  typedef typename Get_functor<Kbase, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CPI;


  typedef Point_d                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename Get_type<R_, Point_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename Get_type<Kbase, Point_tag>::type      Rep;
  //typedef typename CGAL::decay<typename boost::result_of<CPI(Rep,Begin_tag)>::type>::type Cartesian_const_iterator;

  const Rep& rep() const noexcept
  {
    return *this;
  }

  Rep& rep() noexcept
  {
    return *this;
  }

  typedef          R_                       R;

#if defined(BOOST_MSVC) && (BOOST_MSVC == 1900)
#  pragma warning(push)
#  pragma warning(disable: 4309)
#endif

  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Point_d> >::value>::type> explicit Point_d(U&&...u)
          : Rep(CPBase()(std::forward<U>(u)...)){}

#if defined(BOOST_MSVC) && (BOOST_MSVC == 1900)
#  pragma warning(pop)
#endif

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//          : Rep(Eval_functor(), std::forward<U>(u)...){}
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

  friend void swap(Self& a, Self& b)
#ifdef __cpp_lib_is_swappable
    noexcept(std::is_nothrow_swappable_v<Rep>)
#endif
    {
      using std::swap;
      swap(a.rep(), b.rep());
    }

  decltype(auto) cartesian(int i)const{
          return CCBase()(rep(),i);
  }
  decltype(auto) operator[](int i)const{
          return CCBase()(rep(),i);
  }

  decltype(auto) cartesian_begin()const{
          return CPI()(rep(),Begin_tag());
  }

  decltype(auto) cartesian_end()const{
          return CPI()(rep(),End_tag());
  }

  int dimension() const {
    typedef typename Get_functor<Kbase, Point_dimension_tag>::type PDBase;
    return PDBase()(rep());
  }

  friend auto operator==(Point_d const&p, Point_d const&q) {
    typedef typename Get_functor<Kbase, Equal_points_tag>::type EPBase;
    return EPBase()(p.rep(), q.rep());
  }

  friend auto operator!=(Point_d const&p, Point_d const&q) { return !(p==q); }

  friend std::ostream& operator <<(std::ostream& os, const Point_d& p)
  {
    auto b = p.cartesian_begin();
    auto e = p.cartesian_end();
    if(is_ascii(os))
    {
      os << p.dimension();
      for(; b != e; ++b){
        os << " " << *b;
      }
    }
    else
    {
      write(os, p.dimension());
      for(; b != e; ++b){
        write(os, *b);
      }
    }
    return os;
  }

  // TODO: test if the stream is binary or text?
  friend std::istream& operator>>(std::istream &is, Point_d & p)
  {
    int dim;
    if( is_ascii(is) )
      is >> dim;
    else
    {
      read(is, dim);
    }

    if(!is) return is;
    std::vector<FT_> coords(dim);
    if(is_ascii(is))
    {
      for(int i=0;i<dim;++i)
        is >> iformat(coords[i]);
    }
    else
    {
      for(int i=0;i<dim;++i)
        read(is, coords[i]);
    }

    // FIXME: with Epeck_d, currently, this stores pointers to coords which will soon be dead.
    if(is)
      p = Point_d(coords.begin(), coords.end());
    return is;
  }
};
#if 0
template <class R_> Point_d<R_>::Point_d(Point_d &)=default;
#endif


//template <class R_>
//Vector_d<R_> operator+(const Vector_d<R_>& v,const Vector_d<R_>& w) const
//{
//        return typename R::template Construct<Sum_of_vectors_tag>::type()(v,w);
//}
//
//template <class R_>
//Vector_d<R_> operator-(const Vector_d<R_>& v,const Vector_d<R_>& w) const
//{
//        return typename R::template Construct<Difference_of_vectors_tag>::type()(v,w);
//}

} //namespace Wrap
} //namespace CGAL
#endif // CGAL_WRAPPER_POINT_D_H
