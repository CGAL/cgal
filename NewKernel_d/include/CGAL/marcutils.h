#ifndef marcutils
#define marcutils

#ifdef CGAL_CXX0X
#include <type_traits>
#include <utility>
#define CGAL_FORWARDABLE(T) T&&
#define CGAL_FORWARD(T,t) std::forward<T>(t)
#define CGAL_MOVE(t) std::move(t)
#define CGAL_CONSTEXPR constexpr
#else
#define CGAL_FORWARDABLE(T) T const&
#define CGAL_FORWARD(T,t) t
#define CGAL_MOVE(t) t
#define CGAL_CONSTEXPR
#endif
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Rational_traits.h>


namespace CGAL {
	// tell a function f(a,b,c) that its real argument is a(b,c)
	struct Eval_functor {};

	// forget the first argument. Useful to make something dependant
	// (and thus usable in SFINAE), although that's not a great design.
	template<class A,class B> struct Second_arg {
		typedef B type; 
	};      

	// like std::forward, except for basic types where it does a cast, to
	// avoid issues with narrowing conversions
#ifdef CGAL_CXX0X
	template<class T,class U,class V>
		typename std::conditional<std::is_arithmetic<T>::value&&std::is_arithmetic<typename std::remove_reference<U>::type>::value,T,U&&>::type
	       	forward_safe(V&& u) { return std::forward<U>(u); }
#else
	template<class T,class U> U const& forward_safe(U const& u) {
		return u;
	}
#endif

#ifdef CGAL_CXX0X
	template<class...> struct Constructible_from_each;
	template<class To,class From1,class...From> struct Constructible_from_each<To,From1,From...>{
		enum { value=std::is_convertible<From1,To>::value&&Constructible_from_each<To,From...>::value };
	};
	template<class To> struct Constructible_from_each<To>{
		enum { value=true };
	};
#else
// currently only used in C++0X code
#endif

	template<class T> struct Scale {
		T const& scale;
		Scale(T const& t):scale(t){}
		template<class FT>
#ifdef CGAL_CXX0X
		auto operator()(FT&& x)const->decltype(scale*std::forward<FT>(x))
#else
		FT operator()(FT const& x)const
#endif
		{
			return scale*CGAL_FORWARD(FT,x);
		}
	};
	template<class NT,class T> struct Divide {
#if !defined(CGAL_CXX0X) || !defined(BOOST_RESULT_OF_USE_DECLTYPE)
		// requires boost > 1.44
		// shouldn't be needed with C++0X
		//template<class> struct result;
		//template<class FT> struct result<Divide(FT)> {
		//	typedef FT type;
		//};
		typedef NT result_type;
#endif
		T const& scale;
		Divide(T const& t):scale(t){}
		template<class FT>
#ifdef CGAL_CXX0X
		//FIXME: gcc complains for Gmpq
		//auto operator()(FT&& x)const->decltype(Rational_traits<NT>().make_rational(std::forward<FT>(x),scale))
		NT operator()(FT&& x)const
#else
		NT operator()(FT const& x)const
#endif
		{
			return Rational_traits<NT>().
				make_rational(CGAL_FORWARD(FT,x),scale);
		}
	};

	template <class NT> struct has_cheap_constructor : boost::is_arithmetic<NT>{};
	template <bool p> struct has_cheap_constructor<Interval_nt<p> > {
		        enum { value=true };
	};

	// like std::multiplies but allows mixing types
	// in C++0x in doesn't need to be a template
	template < class Ret >
	struct multiplies {
		template<class A,class B>
#ifdef CGAL_CXX0X
		auto operator()(A&&a,B&&b)const->decltype(std::forward<A>(a)*std::forward<B>(b))
#else
		Ret operator()(A const& a, B const& b)const
#endif
		{
			return CGAL_FORWARD(A,a)*CGAL_FORWARD(B,b);
		}
	};
	template < class Ret >
	struct division {
		template<class A,class B>
#ifdef CGAL_CXX0X
		auto operator()(A&&a,B&&b)const->decltype(std::forward<A>(a)/std::forward<B>(b))
#else
		Ret operator()(A const& a, B const& b)const
#endif
		{
			return CGAL_FORWARD(A,a)/CGAL_FORWARD(B,b);
		}
	};

#ifdef CGAL_CXX0X
	using std::decay;
#else
	template<class T> struct decay : boost::remove_cv<typename boost::decay<T>::type> {};
#endif
}

#endif
