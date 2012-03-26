#ifndef CGAL_DEFINE_KERNEL_TYPES_H
#define CGAL_DEFINE_KERNEL_TYPES_H
#include <CGAL/functor_tags.h>
#include <CGAL/typeset.h>
#ifdef CGAL_CXX0X
#include <type_traits>
#else
#include <boost/type_traits.hpp>
#endif

namespace CGAL {
  template<class K, class Base=K, class List=typename K::Object_list> struct Define_kernel_types;
  template<class K, class Base>
    struct Define_kernel_types <K, Base, typeset<> > : Base {};
  template<class K>
    struct Define_kernel_types <K, void, typeset<> > {};
  template<class K, class Base, class List>
    struct Define_kernel_types :
      Typedef_tag_type<typename List::head,
        typename K::template Type<typename List::head>::type,
	Define_kernel_types<K, Base, typename List::tail>
      > {};
}
#endif
