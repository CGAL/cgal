# /* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2002.  Permission to copy, use,        *
#  *     modify, sell, and distribute this software is granted provided       *
#  *     this copyright notice appears in all copies.  This software is       *
#  *     provided "as is" without express or implied warranty, and with       *
#  *     no claim at to its suitability for any purpose.                      *
#  *                                                                          *
#  ************************************************************************** */
#
# /* See http://www.boost.org for most recent version. */
#
# ifndef BOOST_PREPROCESSOR_CONFIG_CONFIG_HPP
# define BOOST_PREPROCESSOR_CONFIG_CONFIG_HPP
#
# /* BOOST_PP_CONFIG_FLAGS */
#
# define BOOST_PP_CONFIG_STRICT() 0x0001
# define BOOST_PP_CONFIG_IDEAL() 0x0002
#
# define BOOST_PP_CONFIG_MSVC() 0x0004
# define BOOST_PP_CONFIG_MWCC() 0x0008
# define BOOST_PP_CONFIG_BCC() 0x0010
# define BOOST_PP_CONFIG_EDG() 0x0020
#
# ifndef BOOST_PP_CONFIG_FLAGS
#    if defined(__SPIRIT_PP__) || defined(__MWERKS__) && __MWERKS__ >= 0x3200
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_STRICT())
#    elif defined(__EDG__) || defined(__EDG_VERSION__)
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_EDG() | BOOST_PP_CONFIG_STRICT())
#    elif defined(__MWERKS__)
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_MWCC())
#    elif defined(__BORLANDC__) || defined(__IBMC__) || defined(__IBMCPP__) || defined(__SUNPRO_CC)
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_BCC())
#    elif defined(_MSC_VER)
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_MSVC())
#    else
#        define BOOST_PP_CONFIG_FLAGS() (BOOST_PP_CONFIG_STRICT())
#    endif
# endif
#
# /* BOOST_PP_CONFIG_EXTENDED_LINE_INFO */
#
# ifndef BOOST_PP_CONFIG_EXTENDED_LINE_INFO
#    define BOOST_PP_CONFIG_EXTENDED_LINE_INFO 0
# endif
#
# /* BOOST_PP_CONFIG_ERRORS */
#
# ifndef BOOST_PP_CONFIG_ERRORS
#    ifdef NDEBUG
#        define BOOST_PP_CONFIG_ERRORS 0
#    else
#        define BOOST_PP_CONFIG_ERRORS 1
#    endif
# endif
#
# endif
