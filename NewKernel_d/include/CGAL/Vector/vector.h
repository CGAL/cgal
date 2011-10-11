#ifndef CGAL_VECTOR_VECTOR_H
#define CGAL_VECTOR_VECTOR_H
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/marcutils.h>
#include <vector>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
namespace CGAL {

template<class NT_,class Dim_,class Max_dim_=Dim_> struct Vector_vector {
        typedef NT_ NT;
	typedef std::vector<NT> type;
	struct Constructor {
		struct Dimension {
			type operator()(int d) const {
				return type(d);
			}
		};

		struct Iterator {
			template<typename Iter>
				type operator()(int d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e));
					return type(f,e);
				}
		};

		// unneeded thanks to Iterator_and_last?
#if 0
		struct Iterator_add_one {
			template<typename Iter>
				type operator()(int d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					type a;
					a.reserve(d+1);
					a.insert(a.end(),f,e);
					a.push_back(1);
					return a;
				}
		};
#endif

		struct Iterator_and_last {
			template<typename Iter,typename T>
				type operator()(int d,Iter const& f,Iter const& e,CGAL_FORWARDABLE(T) t) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					type a;
					a.reserve(d+1);
					a.insert(a.end(),f,e);
					a.push_back(CGAL_FORWARD(T,t));
					return a;
				}
		};

		// useless, use a transform_iterator?
#if 0
		struct Iterator_and_last_divide {
			template<typename Iter,typename T>
				type operator()(int d,Iter f,Iter const& e,T const&t) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					type a;
					a.reserve(d+1);
					for(;f!=e;++f){
						a.push_back(*f/t);
					}
					return a;
				}
		};
#endif

		struct Values {
#ifdef CGAL_CXX0X
			template<class...U>
				type operator()(U&&...u) const {
					//TODO: check the right number of {}, g++ accepts one and two
					type a={forward_safe<NT,U>(u)...};
					return a;
				}
#else

#define VAR(Z,N,_) a.push_back(t##N);
#define CODE(Z,N,_) type operator()(BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
	type a; \
	a.reserve(N); \
	BOOST_PP_REPEAT(N,VAR,) \
	return a; \
}
BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE
#undef VAR

#endif
		};

		struct Values_divide {
#ifdef CGAL_CXX0X
			template<class H,class...U>
				type operator()(H const&h,U&&...u) const {
					//TODO: do we want to cast at some point?
					//e.g. to avoid 1/2 in integers
					// ==> use Rational_traits<NT>().make_rational(x,y) ?
					type a={Rational_traits<NT>().make_rational(std::forward<U>(u),h)...};
					return a;
				}
#else

#define VAR(Z,N,_) a.push_back(Rational_traits<NT>().make_rational( t##N ,h));
#define CODE(Z,N,_) template<class H> type \
			operator()(H const&h, BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
				type a; \
				a.reserve(N); \
				BOOST_PP_REPEAT(N,VAR,) \
				return a; \
			}
			BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE
#undef VAR

#endif
		};
	};
	typedef typename type::const_iterator const_iterator;
	static const_iterator vector_begin(type const&a){
		return a.begin();
	}
	static const_iterator vector_end(type const&a){
		return a.end();
	}
	static int size_of_vector(type const&a){
		return a.size();
	}
};


}
#endif

