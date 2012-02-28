#ifndef CGAL_TYPESET_H
#define CGAL_TYPESET_H
#ifdef CGAL_CXX0X
#include <type_traits>
#else
#include <boost/type_traits/conditional.hpp>
#endif

// Sometimes using tuple just to list types is overkill (takes forever to
// instantiate).

namespace CGAL {
#ifdef CGAL_CXX0X
  template<class...> struct typeset;
  template<class T,class...U> struct typeset<T,U...> {
    typedef T head;
    typedef typeset<U...> tail;
    template<class X> struct contains {
      enum { value = std::conditional<std::is_same<T,X>::value,std::true_type,tail::contains<X>>::type::value };
    };
    template<class X> struct add {
      typedef std::conditional<contains<X>::value,typeset<T,U...>,typeset<T,U...,X>>::type type;
    };
  };
  template<> struct typeset<> {
    template<class X> struct contains : boost::false_type {};
    template<class X> struct add {
      typedef typeset<X> type;
    };
  };
#else
    // ???
#endif
}
#endif
