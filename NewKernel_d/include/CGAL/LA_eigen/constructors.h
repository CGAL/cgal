#ifndef ecs_h
#define ecs_h
#include <CGAL/marcutils.h>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

namespace CGAL {
template <class Vector_> struct Construct_eigen {
	typedef Vector_ result_type;
	typedef typename Vector_::Scalar NT;

	struct Dimension {
	result_type operator()(int d) const {
		return result_type(d);
	}
	};

	struct Iterator {
	template<typename Iter>
	result_type operator()(int d,Iter const& f,Iter const& e) const {
		CGAL_assertion(d==std::distance(f,e));
		CGAL_assertion(d<=std::min(result_type::SizeAtCompileTime,result_type::MaxSizeAtCompileTime));
		result_type a(d);
		// TODO: check the right way to do this
		std::copy(f,e,&a[0]);
		return a;
	}
	};

#if 0
	struct Iterator_add_one {
	template<typename Iter>
	result_type operator()(int d,Iter const& f,Iter const& e) const {
		CGAL_assertion(d==std::distance(f,e)+1);
		CGAL_assertion(d<=std::min(result_type::SizeAtCompileTime,result_type::MaxSizeAtCompileTime));
		result_type a;
		std::copy(f,e,&a[0]);
		a[d-1]=1;
		return a;
	}
	};
#endif

	struct Iterator_and_last {
	template<typename Iter,typename T>
	result_type operator()(int d,Iter const& f,Iter const& e,CGAL_FORWARDABLE(T) t) const {
		CGAL_assertion(d==std::distance(f,e)+1);
		CGAL_assertion(d<=std::min(result_type::SizeAtCompileTime,result_type::MaxSizeAtCompileTime));
		result_type a;
		std::copy(f,e,&a[0]);
		a[d-1]=CGAL_FORWARD(T,t);
		return a;
	}
	};

#ifdef CGAL_CXX0X
	struct Initializer_list {
	template<class...U>
	result_type operator()(std::initializer_list<NT> l) const {
		return Iterator()(l.size(),l.begin(),l.end());
	}
	};
#endif

	struct Values {
#ifdef CGAL_CXX0X
	template<class...U>
	result_type operator()(U&&...u) const {
		return Initializer_list()({forward_safe<NT,U>(u)...});
	}
#else

#define CODE(Z,N,_) result_type operator()(BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
		CGAL_assertion(N<=std::min(result_type::SizeAtCompileTime,result_type::MaxSizeAtCompileTime)); \
	        result_type a(N); \
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
	result_type operator()(H const&h,U&&...u) const {
		return Initializer_list()({Rational_traits<NT>().make_rational(std::forward<U>(u),h)...});
	}
#else

#define VAR(Z,N,_) ( t##N / h )
#define CODE(Z,N,_) template <class H> result_type \
	operator()(H const&h, BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
		CGAL_assertion(N<=std::min(result_type::SizeAtCompileTime,result_type::MaxSizeAtCompileTime)); \
	        result_type a(N); \
	        a << BOOST_PP_ENUM(N,VAR,); \
	        return a; \
}
BOOST_PP_REPEAT_FROM_TO(1, 11, CODE, _ )
#undef CODE
#undef VAR

#endif
	};
};
}
#endif
