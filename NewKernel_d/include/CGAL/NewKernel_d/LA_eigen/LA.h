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
//
// Author(s)     : Marc Glisse

#ifndef CGAL_LA_EIGEN_H
#define CGAL_LA_EIGEN_H
#include <CGAL/config.h>
#ifndef CGAL_EIGEN3_ENABLED
#error Requires Eigen
#endif
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <Eigen/Dense>
#include <CGAL/NewKernel_d/LA_eigen/constructors.h>
#include <CGAL/iterator_from_indices.h>

namespace CGAL {

//FIXME: where could we use Matrix_base instead of Matrix?
// Dim_ real dimension
// Max_dim_ upper bound on the dimension
template<class NT_,class Dim_,class Max_dim_=Dim_> struct LA_eigen {
	typedef NT_ NT;
	typedef Dim_ Dimension;
	typedef Max_dim_ Max_dimension;
	enum { dimension = Eigen_dimension<Dimension>::value };
	enum { max_dimension = Eigen_dimension<Max_dimension>::value };
	template< class D2, class D3=D2 >
	struct Rebind_dimension {
	  typedef LA_eigen< NT, D2, D3 > Other;
	};
	template<class,class=void> struct Property : boost::false_type {};
	template<class D> struct Property<Has_vector_plus_minus_tag,D> : boost::true_type {};
	template<class D> struct Property<Has_vector_scalar_ops_tag,D> : boost::true_type {};
	template<class D> struct Property<Has_dot_product_tag,D> : boost::true_type {};

	typedef Eigen::Matrix<NT,Eigen_dimension<Dim_>::value,1,Eigen::ColMajor|Eigen::AutoAlign,Eigen_dimension<Max_dim_>::value,1> Vector;
	typedef Eigen::Matrix<NT,Eigen::Dynamic,1> Dynamic_vector;
	typedef Construct_eigen<Vector> Construct_vector;

#if (EIGEN_WORLD_VERSION>=3)
	typedef NT const* Vector_const_iterator;
#else
	typedef Iterator_from_indices<const type,const NT
#ifndef CGAL_CXX11
	  ,NT
#endif
	  > Vector_const_iterator;
#endif

	template<class Vec_>static Vector_const_iterator vector_begin(Vec_ const&a){
#if (EIGEN_WORLD_VERSION>=3)
	  return &a[0];
#else
	  return Vector_const_iterator(a,0);
#endif
	}

	template<class Vec_>static Vector_const_iterator vector_end(Vec_ const&a){
#if (EIGEN_WORLD_VERSION>=3)
	  // FIXME: Isn't that dangerous if a is an expression and not a concrete vector?
	  return &a[0]+a.size();
#else
	  return Vector_const_iterator(a,a.size());
#endif
	}

	typedef Eigen::Matrix<NT,dimension,dimension,Eigen::ColMajor|Eigen::AutoAlign,max_dimension,max_dimension> Square_matrix;
	typedef Eigen::Matrix<NT,dimension,Eigen::Dynamic,Eigen::ColMajor|Eigen::AutoAlign,max_dimension,Eigen::Dynamic> Dynamic_matrix;
		//TODO: don't pass on the values of Max_* for an expensive NT
                // typedef ... Constructor
                // typedef ... Accessor
#if 0
	private:
	template <class T> class Canonicalize_vector {
		typedef typename Dimension_eigen<T::SizeAtCompileTime>::type S1;
		typedef typename Dimension_eigen<T::MaxSizeAtCompileTime>::type S2;
		public:
		typedef typename Vector<S1,S2>::type type;
	};
	public:
#endif

	template<class Vec_>static int size_of_vector(Vec_ const&v){
		return (int)v.size();
	}

	template<class Vec_>static NT dot_product(Vec_ const&a,Vec_ const&b){
		return a.dot(b);
	}

	template<class Vec_> static int rows(Vec_ const&v) {
		return (int)v.rows();
	}
	template<class Vec_> static int columns(Vec_ const&v) {
		return (int)v.cols();
	}

	template<class Mat_> static NT determinant(Mat_ const&m,bool=false){
		return m.determinant();
	}

	template<class Mat_> static typename
	Same_uncertainty_nt<CGAL::Sign, NT>::type
	sign_of_determinant(Mat_ const&m,bool=false)
	{
		return CGAL::sign(m.determinant());
	}

	template<class Mat_> static int rank(Mat_ const&m){
		// return m.rank();
		// This one uses sqrt so cannot be used with Gmpq
		// TODO: use different algo for different NT?
		// Eigen::ColPivHouseholderQR<Mat_> decomp(m);
		Eigen::FullPivLU<Mat_> decomp(m);
		// decomp.setThreshold(0);
		return static_cast<int>(decomp.rank());
	}

	// m*a==b
	template<class DV, class DM, class V>
	static bool solve(DV&a, DM const&m, V const& b){
		//a = m.colPivHouseholderQr().solve(b);
		a = m.fullPivLu().solve(b);
		return b.isApprox(m*a);
	}

	static Dynamic_matrix basis(Dynamic_matrix const&m){
		return m.fullPivLu().image(m);
	}

	template<class Vec1,class Vec2> static Vector homogeneous_add(Vec1 const&a,Vec2 const&b){
		//TODO: use compile-time size when available
		int d=a.size();
		Vector v(d);
		v << b[d-1]*a.topRows(d-1)+a[d-1]*b.topRows(d-1), a[d-1]*b[d-1];
		return v;
	}

	template<class Vec1,class Vec2> static Vector homogeneous_sub(Vec1 const&a,Vec2 const&b){
		int d=a.size();
		Vector v(d);
		v << b[d-1]*a.topRows(d-1)-a[d-1]*b.topRows(d-1), a[d-1]*b[d-1];
		return v;
	}

	template<class Vec1,class Vec2> static std::pair<NT,NT> homogeneous_dot_product(Vec1 const&a,Vec2 const&b){
		int d=a.size();
		return make_pair(a.topRows(d-1).dot(b.topRows(d-1)), a[d-1]*b[d-1]);
	}

};
}
#endif
