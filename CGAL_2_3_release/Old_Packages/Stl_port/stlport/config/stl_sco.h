// SCO UDK 7 compiler (UnixWare 7x, OSR 5, UnixWare 2x)


//#  define __SGI_STL_NO_ARROW_OPERATOR 1
//#  define __STL_NO_NEW_STYLE_CASTS 1

// #  define __STL_NEED_MUTABLE 1

.. #  define __STL_NO_BAD_ALLOC 1

// allocator::rebind used not to work properly
// #  define __STL_NO_MEMBER_TEMPLATE_CLASSES 1
// #  define __STL_NO_MEMBER_TEMPLATE_KEYWORD 1

#  define __STL_NO_FRIEND_TEMPLATES 1
#  define __STL_NO_QUALIFIED_FRIENDS 1


// #  define __STL_NO_DEFAULT_NON_TYPE_PARAM 1

//#  define __STL_HAS_NO_NEW_IOSTREAMS 1
//#  define __STL_HAS_NO_NEW_C_HEADERS 1

// ???
//#  define __STL_STATIC_CONST_INIT_BUG 1

// ???
//#  define __STL_LINK_TIME_INSTANTIATION 1

// ???
#  define __STL_NO_TEMPLATE_CONVERSIONS 1

#     define __STL_LONG_LONG 1

#     if defined(_REENTRANT)
#           define _UITHREADS     /* if      UnixWare < 7.0.1 */
#           define __STL_UITHREADS
#     endif /* _REENTRANT */
