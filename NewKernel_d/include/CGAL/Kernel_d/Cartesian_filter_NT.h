#ifndef CGAL_KERNEL_D_CARTESIAN_FILTER_NT_H
#define CGAL_KERNEL_D_CARTESIAN_FILTER_NT_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Cartesian_change_FT.h>
#include <CGAL/internal/Exact_type_selector.h>

namespace CGAL {

template < typename Base_ >
struct Cartesian_filter_NT : public Base_
{
    CGAL_CONSTEXPR Cartesian_filter_NT(){}
    CGAL_CONSTEXPR Cartesian_filter_NT(int d):Base_(d){}
    typedef Base_ Kernel_base;
    typedef Cartesian_change_FT<Kernel_base,Interval_nt_advanced> K1;
    typedef typename internal::Exact_type_selector<typename Kernel_base::FT>::Type  Exact_nt;
    typedef Cartesian_change_FT<Kernel_base,Exact_nt> K2;

    template<class T,class D=void,class=typename map_functor_type<T>::type> struct Functor :
       Kernel_base::template Functor<T,D> {};
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
	    struct type {
		    //TODO: use compression (derive from a compressed_pair?)
		    typedef typename K1::template Functor<T>::type P1; P1 p1;
		    typedef typename K2::template Functor<T>::type P2; P2 p2;
		    typedef typename P2::result_type result_type;
		    type(){}
		    type(Cartesian_filter_NT const&k):p1(reinterpret_cast<K1 const&>(k)),p2(reinterpret_cast<K2 const&>(k)){}
		    //FIXME: if predicate's constructor takes a kernel as argument, how do we translate that? reinterpret_cast is really ugly and possibly unsafe.

#ifdef CGAL_CXX0X
		    template<class...U> result_type operator()(U&&...u)const{
			    {
			          Protect_FPU_rounding<true> p;
				  try {
					  typename P1::result_type res=p1(u...); // don't forward as u may be reused
					  if(is_certain(res)) return get_certain(res);
				  } catch (Uncertain_conversion_exception) {}
			    }
			    return p2(std::forward<U>(u)...);
		    }
#else
		    result_type operator()()const{ // does it make sense to have 0 argument?
			    {
			          Protect_FPU_rounding<true> p;
				  try {
					  typename P1::result_type res=p1();
					  if(is_certain(res)) return get_certain(res);
				  } catch (Uncertain_conversion_exception) {}
			    }
			    return p2();
		    }
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> result_type operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t))const{ \
			    { \
			          Protect_FPU_rounding<true> p; \
				  try { \
					  typename P1::result_type res=p1(BOOST_PP_ENUM_PARAMS(N,t)); \
					  if(is_certain(res)) return get_certain(res); \
				  } catch (Uncertain_conversion_exception) {} \
			    } \
			    return p2(BOOST_PP_ENUM_PARAMS(N,t)); \
		    }
		    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE

#endif
	    };
    };
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_FILTER_NT_H
