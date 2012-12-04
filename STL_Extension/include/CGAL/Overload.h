#ifndef CGAL_OVERLOAD_H
#define CGAL_OVERLOAD_H

#include <CGAL/compiler_config.h>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)

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
 
template<typename... Args>
struct Overload : public Overload_helper<std::tuple<Args...>, std::tuple_size<std::tuple<Args...> >::value, Overload<Args...> >
{
  typedef std::tuple<Args...> Tuple_type;
  typedef typename std::tuple_element<0, Tuple_type>::type::result_type result_type;


  Tuple_type t;

  Overload(Args&&... args) : t(std::forward<Args>(args)...) {}
  explicit Overload(const Tuple_type& t) : t(t) {}

  Overload(Tuple_type&& t) : t(t) {}
};

template<typename... Args>
Overload<Args...> make_overload(Args&&... args) { return Overload<Args...>{ std::forward_as_tuple(args...) }; }

} // namespace CGAL
#endif // !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)

#endif /* CGAL_OVERLOAD_H */

