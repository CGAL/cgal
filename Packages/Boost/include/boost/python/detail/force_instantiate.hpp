// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef FORCE_INSTANTIATE_DWA200265_HPP
# define FORCE_INSTANTIATE_DWA200265_HPP

namespace boost { namespace python { namespace detail { 

// Allows us to force the argument to be instantiated without
// incurring unused variable warnings

# if !defined(BOOST_MSVC) || BOOST_MSVC == 1200 || _MSC_FULL_VER > 13102196

template <class T>
inline void force_instantiate(T const&) {}

# else

#  pragma optimize("g", off)
inline void force_instantiate_impl(...) {}
#  pragma optimize("", on)
template <class T>
inline void force_instantiate(T const& x)
{
    detail::force_instantiate_impl(&x);
}
# endif

}}} // namespace boost::python::detail

#endif // FORCE_INSTANTIATE_DWA200265_HPP
