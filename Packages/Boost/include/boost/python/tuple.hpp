// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TUPLE_20020706_HPP
#define TUPLE_20020706_HPP

# include <boost/python/detail/prefix.hpp>

#include <boost/python/object.hpp>
#include <boost/python/converter/pytype_object_mgr_traits.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

namespace boost { namespace python {

namespace detail
{
  struct BOOST_PYTHON_DECL tuple_base : object
  {
   protected:
      tuple_base();
      tuple_base(object_cref sequence);
      
      BOOST_PYTHON_FORWARD_OBJECT_CONSTRUCTORS(tuple_base, object)

   private:
      static detail::new_reference call(object const&);
  };
}

class tuple : public detail::tuple_base
{
    typedef detail::tuple_base base;
 public:
    tuple() {}

    template <class T>
    explicit tuple(T const& sequence)
        : base(object(sequence))
    {
    }

 public: // implementation detail -- for internal use only
    BOOST_PYTHON_FORWARD_OBJECT_CONSTRUCTORS(tuple, base)
};

//
// Converter Specializations    // $$$ JDG $$$ moved here to prevent
//                              // G++ bug complaining specialization
                                // provided after instantiation
namespace converter
{
  template <>
  struct object_manager_traits<tuple>
      : pytype_object_manager_traits<&PyTuple_Type,tuple>
  {
  };
}

// for completeness
inline tuple make_tuple() { return tuple(); }

# define BOOST_PP_ITERATION_PARAMS_1 (3, (1, BOOST_PYTHON_MAX_ARITY, <boost/python/detail/make_tuple.hpp>))
# include BOOST_PP_ITERATE()

}}  // namespace boost::python

#endif

