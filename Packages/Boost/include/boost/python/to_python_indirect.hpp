// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TO_PYTHON_INDIRECT_DWA200221_HPP
# define TO_PYTHON_INDIRECT_DWA200221_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/type_traits/object_traits.hpp>
# include <boost/python/object/pointer_holder.hpp>
# include <boost/python/object/instance.hpp>
# include <boost/python/converter/registered.hpp>
# include <boost/python/detail/unwind_type.hpp>
# include <boost/python/detail/none.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/python/object/make_ptr_instance.hpp>
# include <memory>

namespace boost { namespace python {

template <class T, class MakeHolder>
struct to_python_indirect
{
    PyObject* operator()(T ptr) const;
 private:
    static PyTypeObject* type();
};

//
// implementations
//
namespace detail
{
  struct make_owning_holder
  {
      typedef PyObject* result_type;
      template <class T>
      static result_type execute(T* p)
      {
          // can't use auto_ptr with Intel 5 and VC6 Dinkum library
          // for some reason. We get link errors against the auto_ptr
          // copy constructor.
# if defined(__ICL) && __ICL < 600 
          typedef boost::shared_ptr<T> smart_pointer;
# else 
          typedef std::auto_ptr<T> smart_pointer;
# endif
          typedef objects::pointer_holder<smart_pointer, T> holder_t;

          smart_pointer ptr(p);
          return objects::make_ptr_instance<T, holder_t>::execute(ptr);
      }
  };

  struct make_reference_holder
  {
      typedef PyObject* result_type;
      template <class T>
      static result_type execute(T* p)
      {
          typedef objects::pointer_holder<T*, T> holder_t;
          return objects::make_ptr_instance<T, holder_t>::execute(p);
      }
  };

  struct get_pointer_class
  {
      typedef PyTypeObject* result_type;
      template <class T>
      static result_type execute(T* p)
      {
          BOOST_STATIC_ASSERT(is_class<T>::value);
          return converter::registered<T>::converters.class_object;
      }
  };

  // null_pointer_to_none -- return none() for null pointers and 0 for all other types/values
  // 
  // Uses simulated partial ordering
  template <class T>
  inline PyObject* null_pointer_to_none(T&, int)
  {
      return 0;
  }

  // overload for pointers
  template <class T>
  inline PyObject* null_pointer_to_none(T* x, long)
  {
      return x == 0 ? python::detail::none() : 0;
  }
}

template <class T, class MakeHolder>
inline PyObject* to_python_indirect<T,MakeHolder>::operator()(T x) const
{
    BOOST_STATIC_ASSERT(is_pointer<T>::value || is_reference<T>::value);
        
    PyObject* const null_result = detail::null_pointer_to_none(x, 1L);
    if (null_result != 0)
        return null_result;

    return detail::unwind_type<MakeHolder>(x);
}

template <class T, class MakeHolder>
inline PyTypeObject* to_python_indirect<T,MakeHolder>::type()
{
    return detail::unwind_type<detail::get_pointer_class,T>();
}

}} // namespace boost::python

#endif // TO_PYTHON_INDIRECT_DWA200221_HPP
