// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef BUILTIN_CONVERTERS_DWA2002124_HPP
# define BUILTIN_CONVERTERS_DWA2002124_HPP
# include <boost/python/detail/prefix.hpp>
# include <boost/python/detail/none.hpp>
# include <boost/python/handle.hpp>
# include <boost/implicit_cast.hpp>
# include <string>
# include <complex>
# include <boost/limits.hpp>

// Since all we can use to decide how to convert an object to_python
// is its C++ type, there can be only one such converter for each
// type. Therefore, for built-in conversions we can bypass registry
// lookups using explicit specializations of arg_to_python and
// result_to_python.

namespace boost { namespace python {

namespace converter
{
  template <class T> struct arg_to_python;
  BOOST_PYTHON_DECL PyObject* do_return_to_python(char);
  BOOST_PYTHON_DECL PyObject* do_return_to_python(char const*);
  BOOST_PYTHON_DECL PyObject* do_return_to_python(PyObject*);
  BOOST_PYTHON_DECL PyObject* do_arg_to_python(PyObject*);
}

// Provide specializations of to_python_value
template <class T> struct to_python_value;

namespace detail
{
  // Since there's no registry lookup, always report the existence of
  // a converter.
  struct builtin_to_python
  {
      // This information helps make_getter() decide whether to try to
      // return an internal reference or not. I don't like it much,
      // but it will have to serve for now.
      BOOST_STATIC_CONSTANT(bool, uses_registry = false);
  };
}

// Use expr to create the PyObject corresponding to x
# define BOOST_PYTHON_RETURN_TO_PYTHON_BY_VALUE(T, expr)        \
    template <> struct to_python_value<T&>                      \
        : detail::builtin_to_python                             \
    {                                                           \
        inline PyObject* operator()(T const& x) const           \
        {                                                       \
            return (expr);                                      \
        }                                                       \
    };                                                          \
    template <> struct to_python_value<T const&>                \
        : detail::builtin_to_python                             \
    {                                                           \
        inline PyObject* operator()(T const& x) const           \
        {                                                       \
            return (expr);                                      \
        }                                                       \
    };

# define BOOST_PYTHON_ARG_TO_PYTHON_BY_VALUE(T, expr)   \
    namespace converter                                 \
    {                                                   \
      template <> struct arg_to_python< T >             \
        : handle<>                                      \
      {                                                 \
          arg_to_python(T const& x)                     \
            : python::handle<>(expr) {}                 \
      };                                                \
    } 

// Specialize argument and return value converters for T using expr
# define BOOST_PYTHON_TO_PYTHON_BY_VALUE(T, expr)       \
        BOOST_PYTHON_RETURN_TO_PYTHON_BY_VALUE(T,expr)  \
        BOOST_PYTHON_ARG_TO_PYTHON_BY_VALUE(T,expr)

// Specialize converters for signed and unsigned T to Python Int
# define BOOST_PYTHON_TO_INT(T)                                         \
    BOOST_PYTHON_TO_PYTHON_BY_VALUE(signed T, ::PyInt_FromLong(x))      \
    BOOST_PYTHON_TO_PYTHON_BY_VALUE(                                    \
        unsigned T                                                      \
        , static_cast<unsigned long>(x) > static_cast<unsigned long>(   \
                std::numeric_limits<long>::max())                       \
        ? ::PyLong_FromUnsignedLong(x)                                  \
        : ::PyInt_FromLong(x))

// Bool is not signed.
BOOST_PYTHON_TO_PYTHON_BY_VALUE(bool, ::PyInt_FromLong(x))

// note: handles signed char and unsigned char, but not char (see below)
BOOST_PYTHON_TO_INT(char)

BOOST_PYTHON_TO_INT(short)
BOOST_PYTHON_TO_INT(int)
BOOST_PYTHON_TO_INT(long)

// using Python's macro instead of Boost's - we don't seem to get the
// config right all the time.
# ifdef HAVE_LONG_LONG 
BOOST_PYTHON_TO_PYTHON_BY_VALUE(signed BOOST_PYTHON_LONG_LONG, ::PyLong_FromLongLong(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(unsigned BOOST_PYTHON_LONG_LONG, ::PyLong_FromUnsignedLongLong(x))
# endif
    
# undef BOOST_TO_PYTHON_INT

BOOST_PYTHON_TO_PYTHON_BY_VALUE(char, converter::do_return_to_python(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(char const*, converter::do_return_to_python(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(std::string, ::PyString_FromStringAndSize(x.data(),implicit_cast<int>(x.size())))
# ifndef BOOST_NO_STD_WSTRING
BOOST_PYTHON_TO_PYTHON_BY_VALUE(std::wstring, ::PyUnicode_FromWideChar(x.data(),implicit_cast<int>(x.size())))
# endif 
BOOST_PYTHON_TO_PYTHON_BY_VALUE(float, ::PyFloat_FromDouble(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(double, ::PyFloat_FromDouble(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(long double, ::PyFloat_FromDouble(x))
BOOST_PYTHON_RETURN_TO_PYTHON_BY_VALUE(PyObject*, converter::do_return_to_python(x))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(std::complex<float>, ::PyComplex_FromDoubles(x.real(), x.imag()))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(std::complex<double>, ::PyComplex_FromDoubles(x.real(), x.imag()))
BOOST_PYTHON_TO_PYTHON_BY_VALUE(std::complex<long double>, ::PyComplex_FromDoubles(x.real(), x.imag()))

# undef BOOST_PYTHON_RETURN_TO_PYTHON_BY_VALUE
# undef BOOST_PYTHON_ARG_TO_PYTHON_BY_VALUE
# undef BOOST_PYTHON_TO_PYTHON_BY_VALUE
# undef BOOST_PYTHON_TO_INT
    
namespace converter
{ 

  void initialize_builtin_converters();

}

}} // namespace boost::python::converter

#endif // BUILTIN_CONVERTERS_DWA2002124_HPP
