// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef POINTER_TYPE_ID_DWA2002222_HPP
# define POINTER_TYPE_ID_DWA2002222_HPP

# include <boost/python/type_id.hpp>
# include <boost/type_traits/composite_traits.hpp>

namespace boost { namespace python { namespace converter { 

namespace detail
{
  template <bool is_ref = false>
  struct pointer_typeid_select
  {
      template <class T>
      static inline type_info execute(T*(*)() = 0)
      {
          return type_id<T>();
      }
  };

  template <>
  struct pointer_typeid_select<true>
  {
      template <class T>
      static inline type_info execute(T* const volatile&(*)() = 0)
      {
          return type_id<T>();
      }
    
      template <class T>
      static inline type_info execute(T*volatile&(*)() = 0)
      {
          return type_id<T>();
      }
    
      template <class T>
      static inline type_info execute(T*const&(*)() = 0)
      {
          return type_id<T>();
      }

      template <class T>
      static inline type_info execute(T*&(*)() = 0)
      {
          return type_id<T>();
      }
  };
}

// Usage: pointer_type_id<T>()
//
// Returns a type_info associated with the type pointed
// to by T, which may be a pointer or a reference to a pointer.
template <class T>
type_info pointer_type_id(T(*)() = 0)
{
    return detail::pointer_typeid_select<
          is_reference<T>::value
        >::execute((T(*)())0);
}

}}} // namespace boost::python::converter

#endif // POINTER_TYPE_ID_DWA2002222_HPP
