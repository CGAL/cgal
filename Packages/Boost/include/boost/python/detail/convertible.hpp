// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CONVERTIBLE_DWA2002614_HPP
# define CONVERTIBLE_DWA2002614_HPP

# if defined(__EDG_VERSION__) && __EDG_VERSION__ <= 241
#  include <boost/mpl/if.hpp>
#  include <boost/type_traits/conversion_traits.hpp>
# endif 

// Supplies a runtime is_convertible check which can be used with tag
// dispatching to work around the Metrowerks Pro7 limitation with boost::is_convertible
namespace boost { namespace python { namespace detail { 

typedef char* yes_convertible;
typedef int* no_convertible;

template <class Target>
struct convertible
{
# if !defined(__EDG_VERSION__) || __EDG_VERSION__ > 241 || __EDG_VERSION__ == 238
    static inline no_convertible check(...) { return 0; }
    static inline yes_convertible check(Target) { return 0; }
# else
    template <class X>
    static inline typename mpl::if_c<
        is_convertible<X,Target>::value
        , yes_convertible
        , no_convertible
        >::type check(X const&) { return 0; }
# endif 
};

}}} // namespace boost::python::detail

#endif // CONVERTIBLE_DWA2002614_HPP
