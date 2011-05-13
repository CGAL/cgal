#ifndef CGAL_KERNEL_D_CARTESIAN_WRAP_H
#define CGAL_KERNEL_D_CARTESIAN_WRAP_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Wrapper/Point_d.h>
#include <CGAL/Kernel_d/Wrapper/Vector_d.h>
#include <CGAL/Kernel_d/Wrapper/Segment_d.h>

namespace CGAL {

template < typename Base_ >
struct Cartesian_wrap : public Base_
{
    CGAL_CONSTEXPR Cartesian_wrap(){}
    CGAL_CONSTEXPR Cartesian_wrap(int d):Base_(d){}
    typedef Base_ Kernel_base;
    typedef Cartesian_wrap Self;

    template <class T,bool=false> struct map_type;
#define CGAL_Kernel_obj(X) typedef X##_d<Cartesian_wrap> X; \
    template<bool b> struct map_type<X##_tag,b> { typedef X type; };
#include <CGAL/Kernel_d/interface_macros.h>

    //TODO: adapt all functors
    //TODO: safely apply .rep() to the arguments (and transforming_iterator)
    template<class T,int i=0> struct Construct {
	    typedef typename Kernel_base::template Construct<T>::type B;
	    struct type {
		    typedef typename map_result_tag<T>::type result_tag;
		    typedef typename map_kernel_obj<Self,result_tag>::type result_type;
#ifdef CGAL_CXX0X
		    template<class...U> result_type operator()(U&&...u)const{
			    return result_type(Eval_functor(),B(),std::forward<U>(u)...);
		    }
#else
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> result_type \
		    operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const&u))const{ \
			    return result_type(Eval_functor(),B(),BOOST_PP_ENUM_PARAMS(N,u)); \
		    }
		    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#endif
	    };
    };
    template<int i> struct Construct<Construct_point_cartesian_const_iterator_tag,i> {
	    typedef typename Kernel_base::template Construct<Construct_point_cartesian_const_iterator_tag>::type type;
    };
    template<int i> struct Construct<Construct_vector_cartesian_const_iterator_tag,i> {
	    typedef typename Kernel_base::template Construct<Construct_vector_cartesian_const_iterator_tag>::type type;
    };
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_WRAP_H
