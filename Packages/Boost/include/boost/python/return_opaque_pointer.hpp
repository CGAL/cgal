// Copyright Gottfried Ganﬂauge 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
/*
 * Generic Return value converter generator for opaque C++-pointers
 */
# ifndef RETURN_OPAQUE_POINTER_HPP_
# define RETURN_OPAQUE_POINTER_HPP_

# include <boost/python/detail/prefix.hpp>
# include <boost/python/opaque_pointer_converter.hpp>
# include <boost/python/detail/indirect_traits.hpp>
# include <boost/mpl/if.hpp>

namespace boost { namespace python {

namespace detail
{
  template <class Pointer>
  struct opaque_conversion_holder
  {
      inline PyObject *operator()(Pointer p) const
      {
          static opaque_pointer_converter<Pointer> converter (
              typeid (Pointer).name());

          return converter.convert(p);
      }
  };

  template <class R>
  struct return_opaque_pointer_requires_a_pointer_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;
}
    
struct return_opaque_pointer
{
    template <class R>
    struct apply
    {
        BOOST_STATIC_CONSTANT(
            bool, ok = is_pointer<R>::value);
        
        typedef typename mpl::if_c<
            ok
          , detail::opaque_conversion_holder<R>
          , detail::return_opaque_pointer_requires_a_pointer_type<R>
        >::type type;
    };
};
}} // namespace boost::python
# endif // RETURN_OPAQUE_POINTER_HPP_
