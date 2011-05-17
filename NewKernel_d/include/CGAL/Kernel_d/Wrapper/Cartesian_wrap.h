#ifndef CGAL_KERNEL_D_CARTESIAN_WRAP_H
#define CGAL_KERNEL_D_CARTESIAN_WRAP_H

#include <CGAL/basic.h>
#include <CGAL/is_iterator.h>
#include <CGAL/Kernel_d/Wrapper/Point_d.h>
#include <CGAL/Kernel_d/Wrapper/Vector_d.h>
#include <CGAL/Kernel_d/Wrapper/Segment_d.h>
#include <boost/mpl/or.hpp>

namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Is_wrapper)
template<class T,bool=has_Is_wrapper<T>::value> struct Is_wrapper {
	enum { value=false };
	typedef Tag_false type;
};
template<class T> struct Is_wrapper<T,true> {
	typedef typename T::Is_wrapper type;
	enum { value=type::value };
};

template<class T,bool=is_iterator<T>::value> struct Is_wrapper_iterator {
	enum { value=false };
	typedef Tag_false type;
};
template<class T> struct Is_wrapper_iterator<T,true> :
	Is_wrapper<typename std::iterator_traits<T>::value_type>
{ };

struct Forward_rep {
//TODO: make a good C++0X version with perfect forwarding
//#ifdef CGAL_CXX0X
//template <class T,class=typename std::enable_if<!Is_wrapper<typename std::decay<T>::type>::value&&!Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//T&& operator()(typename std::remove_reference<T>::type&& t) const {return static_cast<T&&>(t);};
//template <class T,class=typename std::enable_if<!Is_wrapper<typename std::decay<T>::type>::value&&!Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//T&& operator()(typename std::remove_reference<T>::type& t) const {return static_cast<T&&>(t);};
//
//template <class T,class=typename std::enable_if<Is_wrapper<typename std::decay<T>::type>::value>::type>
//typename Type_copy_cvref<T,typename std::decay<T>::type::Rep>::type&&
//operator()(T&& t) const {
//	return static_cast<typename Type_copy_cvref<T,typename std::decay<T>::type::Rep>::type&&>(t.rep());
//};
//
//template <class T,class=typename std::enable_if<Is_wrapper_iterator<typename std::decay<T>::type>::value>::type>
//transforming_iterator<Forward_rep,typename std::decay<T>::type>
//operator()(T&& t) const {
//	return make_transforming_iterator(std::forward<T>(t),Forward_rep());
//};
//#else
template <class T,bool=Is_wrapper<T>::value,bool=Is_wrapper_iterator<T>::value> struct result_;
template <class T> struct result_<T,false,false>{typedef T type;};
template <class T> struct result_<T,true,false>{typedef typename decay<T>::type::Rep type;};
template <class T> struct result_<T,false,true>{typedef transforming_iterator<Forward_rep,typename decay<T>::type> type;};
template<class> struct result;
template<class T> struct result<Forward_rep(T)> : result_<T> {};

template <class T> typename boost::disable_if<boost::mpl::or_<Is_wrapper<T>,Is_wrapper_iterator<T> >,T>::type const& operator()(T const& t) const {return t;}

template <class T> typename boost::enable_if<Is_wrapper<T>,T>::type::Rep const& operator()(T const& t) const {return t.rep();}

template <class T> transforming_iterator<Forward_rep,typename boost::enable_if<Is_wrapper_iterator<T>,T>::type> operator()(T const& t) const {return make_transforming_iterator(t,Forward_rep());}
//#endif
};
}

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

    //Translate the arguments
    template<class T,class D=void,class=typename map_functor_type<T>::type,bool=boost::is_same<typename Kernel_base::template Functor<T>::type,Null_functor>::value> struct Functor {
	    typedef typename Kernel_base::template Functor<T>::type B;
	    struct type {
		    typedef typename B::result_type result_type;
#ifdef CGAL_CXX0X
		    template<class...U> result_type operator()(U&&...u)const{
			    return B()(internal::Forward_rep()(u)...);
		    }
#else
#define VAR(Z,N,_) internal::Forward_rep()(u##N)
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> result_type \
		    operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const&u))const{ \
			    return B()(BOOST_PP_ENUM(N,VAR,)); \
		    }
		    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
#undef VAR
#endif
	    };
    };

    //Translate both the arguments and the result
    template<class T,class D,class C> struct Functor<T,D,C,true> {
	    typedef Null_functor type;
    };

    template<class T,class D> struct Functor<T,D,Construct_tag,false> {
	    typedef typename Kernel_base::template Functor<T>::type B;
	    struct type {
		    typedef typename map_result_tag<T>::type result_tag;
		    typedef typename map_kernel_obj<Self,result_tag>::type result_type;
#ifdef CGAL_CXX0X
		    template<class...U> result_type operator()(U&&...u)const{
			    return result_type(Eval_functor(),B(),internal::Forward_rep()(u)...);
		    }
#else
#define VAR(Z,N,_) internal::Forward_rep()(u##N)
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> result_type \
		    operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const&u))const{ \
			    return result_type(Eval_functor(),B(),BOOST_PP_ENUM(N,VAR,)); \
		    }
		    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
#undef VAR
#endif
	    };
    };

};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_WRAP_H
