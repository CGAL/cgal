#ifndef _COMPOSITE_VISITOR_H_
#define _COMPOSITE_VISITOR_H_

#include <CGAL/compiler_config.h>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#include <tuple>

namespace CGAL {

template<typename T, std::size_t i, class U >
struct overload_helper : public overload_helper<T, i - 1, U> {
  typedef typename std::tuple_element<i - 1, T>::type inner_type;
 
  typename inner_type::result_type operator()(typename inner_type::argument_type x) {
    return std::get<i - 1>(static_cast<U*>(this)->t)(x); }
 
  using overload_helper<T, i - 1, U>::operator();
};
 
template<typename T, class U>
struct overload_helper<T, 1, U> {
  typedef typename std::tuple_element<0 ,T>::type inner_type;
 
  typename inner_type::result_type operator()(typename inner_type::argument_type x) {
    return std::get<0>(static_cast<U*>(this)->t)(x); }
};
 
template<typename T>
struct overload : public overload_helper<T, std::tuple_size<T>::value, overload<T> >
{
  // 
  typedef typename std::tuple_element<0, T>::type::result_type result_type;

  T t;
  overload(T&& t) : t(t) {}
  overload(const T& t) : t(t) {}
};

template<typename T>
overload<T> make_overload(T&& t) { return overload<T>(t); }

} // namespace CGAL
#endif // !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#endif /* _COMPOSITE_VISITOR_H_ */

