#ifndef CGAL_KD_TYPE_AFF_TRANSFORMATION_H
#define CGAL_KD_TYPE_AFF_TRANSFORMATION_H
#include <CGAL/store_kernel.h>
#include <boost/preprocessor/repetition.hpp>

// Dummy, that's all the Kernel_d concept requires, so a useful class will wait.

namespace CGAL {
template<class R_>
struct Aff_transformation {
  typedef R_ R;
};
namespace CartesianDKernelFunctors {
template<class R_> struct Construct_aff_transformation {
  CGAL_FUNCTOR_INIT_IGNORE(Construct_aff_transformation)
  typedef R_ R;
  typedef typename Get_type<R, Aff_transformation_tag>::type result_type;
#ifdef CGAL_CXX0X
  template<class...T>
  result_type operator()(T&&...)const{return result_type();}
#else
  result_type operator()()const{
    return result_type();
  }
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> \
  result_type operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const& BOOST_PP_INTERCEPT))const{ \
    return result_type(); \
  }
  BOOST_PP_REPEAT_FROM_TO(1, 9, CODE, _ )
#undef CODE

#endif
};
}
CGAL_KD_DEFAULT_TYPE(Aff_transformation_tag,(CGAL::Aff_transformation<K>),(),());
CGAL_KD_DEFAULT_FUNCTOR(Construct_ttag<Aff_transformation_tag>,(CartesianDKernelFunctors::Construct_aff_transformation<K>),(Aff_transformation_tag),());

}
#endif
