// (C) Copyright Thomas Witt    2002.
// Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

// no include guard multiple inclusion intended

//
// This is a temporary workaround until the bulk of this is
// available in boost config.
// 23/02/03 thw
//

#undef BOOST_NO_IS_CONVERTIBLE
#undef BOOST_NO_IS_CONVERTIBLE_TEMPLATE
#undef BOOST_NO_STRICT_ITERATOR_INTEROPERABILITY
#undef BOOST_ARG_DEPENDENT_TYPENAME
#undef BOOST_NO_LVALUE_RETURN_DETECTION
#undef BOOST_NO_ONE_WAY_ITERATOR_INTEROP

#ifdef BOOST_ITERATOR_CONFIG_DEF
# undef BOOST_ITERATOR_CONFIG_DEF
#else
# error missing or nested #include config_def
#endif 
