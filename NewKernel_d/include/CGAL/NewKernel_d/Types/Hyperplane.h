// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_TYPE_HYPERPLANE_H
#define CGAL_KD_TYPE_HYPERPLANE_H
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/boost/iterator/counting_iterator.hpp>
namespace CGAL {
template <class R_> class Hyperplane {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Vector_tag>::type	Vector_;
	Vector_ v_;
	FT_ s_;

	public:
	Hyperplane(Vector_ const&v, FT_ const&s): v_(v), s_(s) {}
	// TODO: Add a piecewise constructor?

	Vector_ const& orthogonal_vector()const{return v_;}
	FT_ translation()const{return s_;}
};
namespace CartesianDKernelFunctors {
template <class R_> struct Construct_hyperplane : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_hyperplane)
  typedef typename Get_type<R_, Hyperplane_tag>::type	result_type;
  typedef typename Get_type<R_, Point_tag>::type	Point;
  typedef typename Get_type<R_, Vector_tag>::type	Vector;
  typedef typename Get_type<R_, FT_tag>::type FT;
  private:
  struct One {
    typedef int result_type;
    template<class T>int const& operator()(T const&)const{
      static const int one = 1;
      return one;
    }
  };
  public:

  result_type operator()(Vector const&a, FT const&b)const{
    return result_type(a,b);
  }
  // Not really needed
  result_type operator()()const{
    typename Get_functor<R_, Construct_ttag<Vector_tag> >::type cv(this->kernel());
    return result_type(cv(),0);
  }

  template <class Iter>
  result_type through(Iter f, Iter e)const{
    typedef typename R_::LA LA;
    typedef typename R_::Default_ambient_dimension D1;
    typedef typename R_::Max_ambient_dimension D2;
    typedef typename Increment_dimension<D1>::type D1i;
    typedef typename Increment_dimension<D2>::type D2i;

    typedef Eigen::Matrix<FT, Eigen_dimension<D1>::value, Eigen_dimension<D1i>::value,
	      Eigen::ColMajor|Eigen::AutoAlign, Eigen_dimension<D2>::value, Eigen_dimension<D2i>::value> Matrix;
    typedef Eigen::Matrix<FT, Eigen_dimension<D1i>::value, 1,
	      Eigen::ColMajor|Eigen::AutoAlign, Eigen_dimension<D2i>::value, 1> Vec;
    typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
    typename Get_functor<R_, Construct_ttag<Vector_tag> >::type cv(this->kernel());
    typename Get_functor<R_, Point_dimension_tag>::type pd(this->kernel());

    Point const& p0=*f;
    int d = pd(p0);
    Matrix m(d,d+1);
    for(int j=0;j<d;++j)
      m(0,j)=c(p0,j);
    // Write the point coordinates in lines.
    int i;
    for (i=1; ++f!=e; ++i) {
      Point const& p=*f;
      for(int j=0;j<d;++j)
	m(i,j)=c(p,j);
    }
    CGAL_assertion (i == d);
    for(i=0;i<d;++i)
      m(i,d)=-1;
    Eigen::FullPivLU<Matrix> lu(m);
    Vec res = lu.kernel().col(0);
    return this->operator()(cv(d,LA::vector_begin(res),LA::vector_end(res)-1),res(d));
  }
  template <class Iter>
  result_type operator()(Iter f, Iter e, Point const&p, CGAL::Oriented_side s=ON_ORIENTED_BOUNDARY)const{
    result_type ret = through(f, e);
    // I don't really like using ON_ORIENTED_BOUNDARY to mean that we don't care, we might as well not pass 'p' at all.
    if (s == ON_ORIENTED_BOUNDARY)
      return ret;
    typename Get_functor<R_, Oriented_side_tag>::type os(this->kernel());
    CGAL::Oriented_side o = os(ret, p);
    if (o == ON_ORIENTED_BOUNDARY || o == s)
      return ret;
    typename Get_functor<R_, Opposite_vector_tag>::type ov(this->kernel());
    typename Get_functor<R_, Construct_ttag<Vector_tag> >::type cv(this->kernel());
    return this->operator()(ov(ret.orthogonal_vector()), -ret.translation());
  }
};
template <class R_> struct Orthogonal_vector {
  CGAL_FUNCTOR_INIT_IGNORE(Orthogonal_vector)
  typedef typename Get_type<R_, Hyperplane_tag>::type		Hyperplane;
  typedef typename Get_type<R_, Vector_tag>::type const&	result_type;
  result_type operator()(Hyperplane const&s)const{
    return s.orthogonal_vector();
  }
};
template <class R_> struct Hyperplane_translation {
  CGAL_FUNCTOR_INIT_IGNORE(Hyperplane_translation)
  typedef typename Get_type<R_, Hyperplane_tag>::type	Hyperplane;
  typedef typename Get_type<R_, FT_tag>::type result_type;
  // TODO: Is_exact?
  result_type operator()(Hyperplane const&s)const{
    return s.translation();
  }
};
template <class R_> struct Value_at : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Value_at)
  typedef typename Get_type<R_, Hyperplane_tag>::type	Hyperplane;
  typedef typename Get_type<R_, Vector_tag>::type	Vector;
  typedef typename Get_type<R_, Point_tag>::type	Point;
  typedef typename Get_type<R_, FT_tag>::type		FT;
  typedef FT result_type;
  typedef typename Get_functor<R_, Scalar_product_tag>::type	Dot;
  typedef typename Get_functor<R_, Point_to_vector_tag>::type	P2V;
  result_type operator()(Hyperplane const&h, Point const&p)const{
    Dot dot(this->kernel());
    P2V p2v(this->kernel());
    return dot(h.orthogonal_vector(),p2v(p));
    // Use Orthogonal_vector to make it generic?
    // Copy the code from Scalar_product to avoid p2v?
  }
};
}
//TODO: Add a condition that the hyperplane type is the one from this file.
CGAL_KD_DEFAULT_TYPE(Hyperplane_tag,(CGAL::Hyperplane<K>),(Vector_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Hyperplane_tag>,(CartesianDKernelFunctors::Construct_hyperplane<K>),(Vector_tag,Hyperplane_tag),(Opposite_vector_tag,Oriented_side_tag));
CGAL_KD_DEFAULT_FUNCTOR(Orthogonal_vector_tag,(CartesianDKernelFunctors::Orthogonal_vector<K>),(Vector_tag,Hyperplane_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Hyperplane_translation_tag,(CartesianDKernelFunctors::Hyperplane_translation<K>),(Hyperplane_tag),());
CGAL_KD_DEFAULT_FUNCTOR(Value_at_tag,(CartesianDKernelFunctors::Value_at<K>),(Point_tag,Vector_tag,Hyperplane_tag),(Scalar_product_tag,Point_to_vector_tag));
} // namespace CGAL
#endif
