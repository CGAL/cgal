#ifndef CGAL_TYPESET_H
#define CGAL_TYPESET_H
#ifdef CGAL_CXX0X
#include <type_traits>
#else
#include <boost/type_traits.hpp>
#endif

// Sometimes using tuple just to list types is overkill (takes forever to
// instantiate).

namespace CGAL {
#ifdef CGAL_CXX0X
  template<class...> struct typeset;
  template<class H,class...U> struct typeset<H,U...> {
    typedef H head;
    typedef typeset<U...> tail;
    typedef typeset type;
    template<class X> using contains = typename
      std::conditional<
        std::is_same<H,X>::value,
        std::true_type,
        typename tail::template contains<X>
      >::type;
    template<class X> using add = typename
      std::conditional<
        contains<X>::value,
        typeset<H,U...>,
	typeset<H,U...,X>
      >::type;
  };
  template<> struct typeset<> {
    typedef typeset type;
    template<class X> using contains = std::false_type;
    template<class X> using add = typeset<X>;
  };
#else
  template<class,class> struct typeset;
  template<class H=void, class T=typename
    boost::mpl::if_<boost::is_same<H,void>, void, typeset<void, void> >::type >
  struct typeset {
    typedef typeset type;
    typedef H head;
    typedef T tail;
    template<class X> struct contains :
      boost::mpl::if_<boost::is_same<H,X>,boost::true_type,typename tail::template contains<X> >::type
    {};
    template<class X,class=void> struct add :
      //boost::mpl::if_<boost::is_same<H,X>,typeset,typeset<X,typeset> >::type
      typeset<H,typename tail::template add<X>::type>
    {};
    template<class V> struct add<H,V> : typeset {};
  };
  template<> struct typeset<> {
    typedef typeset type;
    template<class X> struct contains : boost::false_type {};
    template<class X> struct add : typeset<X> {};
  };
#endif
}
#endif
