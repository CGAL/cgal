// Copyright Gottfried Ganﬂauge 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
/*
 * Generic Conversion of opaque C++-pointers to a Python-Wrapper.
 */
# ifndef OPAQUE_POINTER_CONVERTER_HPP_
# define OPAQUE_POINTER_CONVERTER_HPP_

# include <boost/python/detail/prefix.hpp>
# include <boost/python/lvalue_from_pytype.hpp>
# include <boost/python/to_python_converter.hpp>
# include <boost/python/detail/dealloc.hpp>
# include <boost/python/detail/none.hpp>
# include <boost/type_traits/remove_pointer.hpp>
# include <boost/type_traits/is_pointer.hpp>

// opaque_pointer_converter --
//
// usage: opaque_pointer_converter<Pointer>("name")
//
// registers to- and from- python conversions for a type Pointer,
// and a corresponding Python type called "name".
//
// Note:
// In addition you need to define specializations for type_id
// on the type pointed to by Pointer using
// BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(Pointee)
//
// For an example see libs/python/test/opaque.cpp
//
namespace boost { namespace python {
    namespace detail {
        template <class R>
        struct opaque_pointer_converter_requires_a_pointer_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
        {}
# endif
        ;
    }

template <class Pointer>
struct opaque_pointer_converter
    : to_python_converter<
          Pointer, opaque_pointer_converter<Pointer> >
{
    BOOST_STATIC_CONSTANT(
        bool, ok = is_pointer<Pointer>::value);
        
    typedef typename mpl::if_c<
        ok
        , Pointer
        , detail::opaque_pointer_converter_requires_a_pointer_type<Pointer>
    >::type ptr_type;

private:
    struct instance;

public:
    explicit opaque_pointer_converter(char const* name)
    {
        type_object.tp_name = const_cast<char *> (name);

        lvalue_from_pytype<
            opaque_pointer_converter<ptr_type>,
            &opaque_pointer_converter<ptr_type>::type_object
        >();
    }

    static PyObject* convert(ptr_type x)
    {
        PyObject *result = 0;
        
        if (x != 0) {
            instance *o = PyObject_New (instance, &type_object);

            o->x   = x;
            result = &o->base_;
        } else {
            result = detail::none();
        }
        
        return (result);
    }

    static typename ::boost::remove_pointer<ptr_type>::type&
    execute(instance &p_)
    {
        return *p_.x;
    }

private:
    static PyTypeObject type_object;

    // This is a POD so we can use PyObject_Del on it, for example.
    struct instance
    {
        PyObject base_;
        ptr_type x;
    };
};

template <class Pointer>
PyTypeObject opaque_pointer_converter<Pointer>::type_object =
{
    PyObject_HEAD_INIT(NULL)
    0,
    0,
    sizeof(typename opaque_pointer_converter<Pointer>::instance),
    0,
    ::boost::python::detail::dealloc
};
}} // namespace boost::python
#  ifdef BOOST_MSVC
// MSC works without this workaround, but needs another one ...
#  define BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(Pointee)      \
     BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(Pointee)
#  else
#   define BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(Pointee)                     \
    namespace boost { namespace python {                                        \
    template<>                                                                  \
    inline type_info type_id<Pointee>(BOOST_PYTHON_EXPLICIT_TT_DEF(Pointee))    \
    {                                                                           \
        return type_info (typeid (Pointee *));                                  \
    }                                                                           \
    template<>                                                                  \
    inline type_info type_id<const volatile Pointee&>(                          \
        BOOST_PYTHON_EXPLICIT_TT_DEF(const volatile Pointee&))                  \
    {                                                                           \
        return type_info (typeid (Pointee *));                                  \
    }                                                                           \
}}
#  endif
# endif    // OPAQUE_POINTER_CONVERTER_HPP_
