//----------------------------------------------------------------------//
//             STLport fix for MSVC
//----------------------------------------------------------------------//
#ifndef CGAL_STL_CONFIGVC7

#ifdef _MSC_VER
#if _MSC_VER < 1300
#error "this is not MSVC 7.0 !"
#else  // VC7

// <NBulteau@jouve.fr> : suppressed "truncated debug info" warning
#   pragma warning(disable:4786)

// decorated name length exceeded warning only once.
#   pragma warning(once:4503)

# include <iterator> // including iterator-traits fixes from vc7/ dir
#define __STL_BEGIN_NAMESPACE namespace std {
#define __STL_END_NAMESPACE }
#endif //VC7
#else  // not _MSC_VER
#error "this is not MSVC !"
#endif // _MSC_VER

#endif // CGAL_STL_CONFIGVC7
