// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TYPE_LIST_DWA2002913_HPP
# define TYPE_LIST_DWA2002913_HPP

# include <boost/config.hpp>
# include <boost/python/detail/preprocessor.hpp>
# include <boost/preprocessor/arithmetic/inc.hpp>

# if BOOST_PYTHON_MAX_ARITY + 2 > BOOST_PYTHON_MAX_BASES
#  define BOOST_PYTHON_LIST_SIZE BOOST_PP_INC(BOOST_PP_INC(BOOST_PYTHON_MAX_ARITY))
# else
#  define BOOST_PYTHON_LIST_SIZE BOOST_PYTHON_MAX_BASES
# endif

// Compute the MPL vector header to use for lists up to BOOST_PYTHON_LIST_SIZE in length
# if BOOST_PYTHON_LIST_SIZE > 48
#  error Arities above 48 not supported by Boost.Python due to MPL internal limit
# elif BOOST_PYTHON_LIST_SIZE > 38
#  include <boost/mpl/vector/vector50.hpp>
# elif BOOST_PYTHON_LIST_SIZE > 28
#  include <boost/mpl/vector/vector40.hpp>
# elif BOOST_PYTHON_LIST_SIZE > 18
#  include <boost/mpl/vector/vector30.hpp>
# elif BOOST_PYTHON_LIST_SIZE > 8
#  include <boost/mpl/vector/vector20.hpp>
# else
#  include <boost/mpl/vector/vector10.hpp>
# endif

# ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#  include <boost/python/detail/type_list_impl.hpp>
# else
#  include <boost/python/detail/type_list_impl_no_pts.hpp>
# endif

#endif // TYPE_LIST_DWA2002913_HPP
