#ifndef CGAL_VECTOR_ARRAY_H
#define CGAL_VECTOR_ARRAY_H
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/marcutils.h>
#include <CGAL/array.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
namespace CGAL {

// May not be safe to use with dim!=max_dim.
// In that case, we should store the real dim next to the array.
template<class NT_,class Dim_,class Max_dim_=Dim_> struct Array_vector {
        typedef NT_ NT;
	typedef Dim_ Dimension;
	typedef Max_dim_ Max_dimension;
	template< class D2, class D3=D2 >
	struct Rebind_dimension {
	  typedef Array_vector< NT, D2, D3 > Other;
	};
	template<class> struct Property : boost::false_type {};

	static const unsigned d_=Max_dim_::value;
	typedef cpp0x::array<NT,d_> Vector;
	struct Construct_vector {
		struct Dimension {
			// Initialize with NaN if possible?
			Vector operator()(unsigned d) const {
				CGAL_assertion(d<=d_);
				return Vector();
			}
		};

		struct Iterator {
			template<typename Iter>
				Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e));
					CGAL_assertion(d<=d_);
					//TODO: optimize for forward iterators
					Vector a;
					std::copy(f,e,a.begin());
					return a;
				}
		};

#if 0
		struct Iterator_add_one {
			template<typename Iter>
				Vector operator()(unsigned d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					CGAL_assertion(d<=d_);
					//TODO: optimize
					Vector a;
					std::copy(f,e,a.begin());
					a.back()=1;
					return a;
				}
		};
#endif

		struct Iterator_and_last {
			template<typename Iter,typename T>
				Vector operator()(unsigned d,Iter const& f,Iter const& e,CGAL_FORWARDABLE(T) t) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					CGAL_assertion(d<=d_);
					//TODO: optimize for forward iterators
					Vector a;
					std::copy(f,e,a.begin());
					a.back()=CGAL_FORWARD(T,t);
					return a;
				}
		};

		struct Values {
#ifdef CGAL_CXX0X
			template<class...U>
				Vector operator()(U&&...u) const {
					static_assert(sizeof...(U)<=d_,"too many arguments");
					Vector a={{forward_safe<NT,U>(u)...}};
					return a;
				}
#else

#define CODE(Z,N,_) Vector operator()(BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
	CGAL_assertion(N<=d_); \
	Vector a={{BOOST_PP_ENUM_PARAMS(N,t)}}; \
	return a; \
}
BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE

#endif
		};

		struct Values_divide {
#ifdef CGAL_CXX0X
			template<class H,class...U>
				Vector operator()(H const& h,U&&...u) const {
					static_assert(sizeof...(U)<=d_,"too many arguments");
					Vector a={{Rational_traits<NT>().make_rational(std::forward<U>(u),h)...}};
					return a;
				}
#else

#define VAR(Z,N,_) Rational_traits<NT>().make_rational( t##N , h)
#define CODE(Z,N,_) template <class H> Vector \
			operator()(H const&h, BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
				CGAL_assertion(N<=d_); \
				Vector a={{BOOST_PP_ENUM(N,VAR,_)}}; \
				return a; \
			}
			BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE
#undef VAR

#endif
		};
	};

	typedef NT const* Vector_const_iterator;
	static Vector_const_iterator vector_begin(Vector const&a){
		return &a[0];
	}
	static Vector_const_iterator vector_end(Vector const&a){
		return &a[0]+d_; // Don't know the real size
	}
	static unsigned size_of_vector(Vector const&a){
		return d_; // Don't know the real size
	}

};

}
#endif
