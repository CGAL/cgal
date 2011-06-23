#ifndef CGAL_OVERLOAD_H
#define CGAL_OVERLOAD_H

#include <CGAL/compiler_config.h>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#include <tuple>

namespace CGAL {

template<typename T, std::size_t i, class U >
struct Overload_helper : public Overload_helper<T, i - 1, U> {
  typedef typename std::tuple_element<i - 1, T>::type inner_type;
 
  typename inner_type::result_type operator()(typename inner_type::argument_type x) {
    return std::get<i - 1>(static_cast<U*>(this)->t)(x); }
 
  using Overload_helper<T, i - 1, U>::operator();
};
 
template<typename T, class U>
struct Overload_helper<T, 1, U> {
  typedef typename std::tuple_element<0 ,T>::type inner_type;
 
  typename inner_type::result_type operator()(typename inner_type::argument_type x) {
    return std::get<0>(static_cast<U*>(this)->t)(x); }
};
 
template<typename T>
struct Overload : public Overload_helper<T, std::tuple_size<T>::value, Overload<T> >
{
  // 
  typedef typename std::tuple_element<0, T>::type::result_type result_type;

  T t;
  Overload(T&& t) : t(t) {}
  Overload(const T& t) : t(t) {}
};

template<typename T>
Overload<T> make_overload(T&& t) { return Overload<T>(t); }

} // namespace CGAL
#endif // !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#endif /* CGAL_OVERLOAD_H */

