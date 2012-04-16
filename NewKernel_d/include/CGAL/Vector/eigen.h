#ifndef CGAL_VECTOR_EIGEN_H
#define CGAL_VECTOR_EIGEN_H

#ifndef CGAL_USE_EIGEN
#error Requires Eigen
#endif
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <Eigen/Dense>
#include <CGAL/iterator_from_indices.h>
#include <CGAL/marcutils.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>


namespace CGAL {

// Dim_ real dimension
// Max_dim_ upper bound on the dimension
template<class NT_, class Dim_, class Max_dim_=Dim_> struct Eigen_vector {
	typedef NT_ NT;
	typedef Eigen::Matrix<NT,Eigen_dimension<Dim_>::value,1,Eigen::ColMajor|Eigen::AutoAlign,Eigen_dimension<Max_dim_>::value,1> type;

	struct Construct_vector {
		private:
			static void check_dim(int d){
				int m = Eigen_dimension<Max_dim_>::value;
				CGAL_assertion((m == Eigen::Dynamic) || (d <= m));
			}
		public:

		struct Dimension {
			// Initialize with NaN if possible?
			type operator()(int d) const {
				check_dim(d);
				return type(d);
			}
		};

		struct Iterator {
			template<typename Iter>
				type operator()(int d,Iter const& f,Iter const& e) const {
					check_dim(d);
					CGAL_assertion(d==std::distance(f,e));
					type a(d);
					// TODO: check the right way to do this
					std::copy(f,e,&a[0]);
					return a;
				}
		};

#if 0
		struct Iterator_add_one {
			template<typename Iter>
				type operator()(int d,Iter const& f,Iter const& e) const {
					check_dim(d);
					CGAL_assertion(d==std::distance(f,e)+1);
					type a(d);
					std::copy(f,e,&a[0]);
					a[d-1]=1;
					return a;
				}
		};
#endif

		struct Iterator_and_last {
			template<typename Iter,typename T>
				type operator()(int d,Iter const& f,Iter const& e,CGAL_FORWARDABLE(T) t) const {
					check_dim(d);
					CGAL_assertion(d==std::distance(f,e)+1);
					type a(d);
					std::copy(f,e,&a[0]);
					a[d-1]=CGAL_FORWARD(T,t);
					return a;
				}
		};

#ifdef CGAL_CXX0X
		struct Initializer_list {
			// Fix T==NT?
			template<class T>
				type operator()(std::initializer_list<T> l) const {
					return Iterator()(l.size(),l.begin(),l.end());
				}
		};
#endif

		struct Values {
#ifdef CGAL_CXX0X
			// TODO avoid going through Initializer_list which may cause extra copies. Possibly use forward_as_tuple.
			template<class...U>
				type operator()(U&&...u) const {
					check_dim(sizeof...(U)); // use static_assert
					return Initializer_list()({forward_safe<NT,U>(u)...});
				}
#else

#define CODE(Z,N,_) type operator()(BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
	check_dim(N); \
	type a(N); \
	a << BOOST_PP_ENUM_PARAMS(N,t); \
	return a; \
}
BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE

#endif
		};

		struct Values_divide {
#ifdef CGAL_CXX0X
			template<class H,class...U>
				type operator()(H const&h,U&&...u) const {
					check_dim(sizeof...(U)); // use static_assert
					return Initializer_list()({Rational_traits<NT>().make_rational(std::forward<U>(u),h)...});
				}
#else

#define VAR(Z,N,_) ( Rational_traits<NT>().make_rational( t##N ,h) )
#define CODE(Z,N,_) template <class H> type \
			operator()(H const&h, BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
				CGAL_assertion(N<=std::min(type::SizeAtCompileTime,type::MaxSizeAtCompileTime)); \
				type a(N); \
				a << BOOST_PP_ENUM(N,VAR,); \
				return a; \
			}
			BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE
#undef VAR

#endif
		};
	};



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

// Really needed?
	template<class Vec_>static int size_of_vector(Vec_ const&v){
		return v.size();
	}

// This complicates matter for little benefice
#if 0
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
#endif
};
}

#endif
