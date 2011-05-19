#ifndef CGAL_STORE_KERNEL_H
#define CGAL_STORE_KERNEL_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_empty.hpp>

namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Do_not_store_kernel)
template<class T,bool=has_Do_not_store_kernel<T>::value> struct Do_not_store_kernel {
	enum { value=false };
	typedef Tag_false type;
};
template<class T> struct Do_not_store_kernel<T,true> {
	typedef typename T::Do_not_store_kernel type;
	enum { value=type::value };
};
}

template<class R_,bool=boost::is_empty<R_>::value||internal::Do_not_store_kernel<R_>::value>
struct Store_kernel {
	Store_kernel(){}
	Store_kernel(R_ const&){}
	enum { kernel_is_stored = false };
	R_ kernel()const{return R_();}
	typedef R_ reference_type;
};
template<class R_>
struct Store_kernel<R_,false> {
	Store_kernel(){ CGAL_warning_msg(true,"I should know my kernel"); }
	Store_kernel(R_ const& r):rp(&r){}
	enum { kernel_is_stored = true };
	R_ const& kernel()const{return *rp;}
	typedef R_ const& reference_type;
	private:
	R_ const* rp;
};

//For a second kernel. TODO: find something more elegant
template<class R_,bool=boost::is_empty<R_>::value||internal::Do_not_store_kernel<R_>::value>
struct Store_kernel2 {
	Store_kernel2(){}
	Store_kernel2(R_ const&){}
	enum { kernel2_is_stored = false };
	R_ kernel2()const{return R_();}
	typedef R_ reference2_type;
};
template<class R_>
struct Store_kernel2<R_,false> {
	Store_kernel2(){ CGAL_warning_msg(true,"I should know my kernel"); }
	Store_kernel2(R_ const& r):rp(&r){}
	enum { kernel2_is_stored = true };
	R_ const& kernel2()const{return *rp;}
	typedef R_ const& reference2_type;
	private:
	R_ const* rp;
};
}
#define CGAL_BASE_INIT(X,Y) \
	X():Y(){} \
	X(R_ const&r):Y(r){}
#define CGAL_FUNCTOR_INIT_STORE(X) CGAL_BASE_INIT(X,Store_kernel<R_>)
#define CGAL_FUNCTOR_INIT_IGNORE(X) \
	X(){} \
	X(R_ const&){}

#endif // CGAL_STORE_KERNEL_H
