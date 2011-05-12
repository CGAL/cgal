#ifndef CGAL_LA_EIGEN_H
#define CGAL_LA_EIGEN_H

#ifndef CGAL_USE_EIGEN
#error Requires Eigen
#endif
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <Eigen/Dense>
#include <CGAL/LA_eigen/constructors.h>
#include <CGAL/iterator_from_indices.h>

namespace CGAL {

//FIXME: where could we use Matrix_base instead of Matrix?
template<class NT_> struct LA_eigen {
	typedef NT_ NT;
	// Dim_ real dimension
	// Max_dim_ upper bound on the dimension
	template<class Dim_,class Max_dim_=Dim_> struct Vector {
		typedef Eigen::Matrix<NT,Eigen_dimension<Dim_>::value,1,Eigen::ColMajor|Eigen::AutoAlign,Eigen_dimension<Max_dim_>::value,1> type;
		typedef Construct_eigen<type> Constructor;
#if (EIGEN_WORLD_VERSION>=3)
		typedef NT const* const_iterator;
#else
		typedef Iterator_from_indices<const type,const NT
#ifndef CGAL_CXX0X
			,NT
#endif
			> const_iterator;
#endif

		template<class Vec_>static const_iterator vector_begin(Vec_ const&a){
#if (EIGEN_WORLD_VERSION>=3)
			return &a[0];
#else
			return const_iterator(a,0);
#endif
		}

		template<class Vec_>static const_iterator vector_end(Vec_ const&a){
#if (EIGEN_WORLD_VERSION>=3)
			return &a[0]+a.size();
#else
			return const_iterator(a,a.size());
#endif
		}

	};
        template<class Rows_,class Cols_=Rows_,class Max_rows_=Rows_,class Max_cols_=Cols_> struct Matrix {
		//TODO: don't pass on the values of Max_* for an expensive NT
                typedef Eigen::Matrix<NT,Eigen_dimension<Rows_>::value,Eigen_dimension<Cols_>::value,Eigen::ColMajor|Eigen::AutoAlign,Eigen_dimension<Max_rows_>::value,Eigen_dimension<Max_cols_>::value> type;
                // typedef ... Constructor
                // typedef ... Accessor
        };
	private:
	template <class T> class Canonicalize_vector {
		typedef typename Dimension_eigen<T::SizeAtCompileTime>::type S1;
		typedef typename Dimension_eigen<T::MaxSizeAtCompileTime>::type S2;
		public:
		typedef typename Vector<S1,S2>::type type;
	};
	public:

	template<class Vec_>static NT dot_product(Vec_ const&a,Vec_ const&b){
		return a.dot(b);
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

	template<class Vec1,class Vec2> static typename Canonicalize_vector<Vec1>::type homogeneous_add(Vec1 const&a,Vec2 const&b){
		//TODO: use compile-time size when available
		int d=a.size();
		typename Canonicalize_vector<Vec1>::type v(d);
		v << b[d-1]*a.topRows(d-1)+a[d-1]*b.topRows(d-1), a[d-1]*b[d-1];
		return v;
	}

	template<class Vec1,class Vec2> static typename Canonicalize_vector<Vec1>::type homogeneous_sub(Vec1 const&a,Vec2 const&b){
		int d=a.size();
		typename Canonicalize_vector<Vec1>::type v(d);
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
