/*
 * Copyright (c) 1999 
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/*
 * Purpose of this file :
 *
 * To hold user-definable portion of STLport settings.
 * If you use STLport in several projects, you may want
 * to place copy of this file modified for each project
 * in project include directory, so that you have different
 * settings for different projects.
 *
 */


//==========================================================
// User-settable macros that control compilation:
//              Features selection
//==========================================================


/* __STL_USE_NEW_IOSTREAMS: if defined, then the STL will use new,
 *   standard-conforming iostreams (e.g. the <iosfwd> header).  If not
 *   defined, the STL will use old cfront-style iostreams (e.g. the
 *   <iostream.h> header). Normally it is being set automatically.
 */
// #define   __STL_USE_NEW_IOSTREAMS	1

/* 
 * Uncomment to suppress using new-style streams even if they are
 * available.
 * Beware - __STL_USE_OWN_NAMESPACE depends on this macro, too.
 * Do that only if you are absolutely sure backwards-compatible 
 * <iostream.h> is not actually a wrapper with <iostream>
 * Hint : In VC++ 5.x, they are not.
 */

// #define   __STL_NO_NEW_IOSTREAMS	1

/*
 * __STL_STD_REBUILD - define this if you are actually rebuilding
 * vendor's C++ support library with STLport installed.
 */
// # define __STL_STD_REBUILD 1

/* 
 * __STL_USE_OWN_NAMESPACE/__STL_NO_OWN_NAMESPACE
 * If defined, STLport uses STLport:: namespace, else std::
 * Defining this helps A LOT in resolving problems with 
 * vendor C++ standard library interaction. 
 * The reason you have to use separate namespace is that new-style IO
 * compiled library may have its own idea about STL stuff (string, vector, etc.),
 * so redefining them in the same namespace would break ODR and may cause
 * undefined behaviour. Rule of thumb is - if new-style iostreams are
 * available, there could be a conflict. Otherwise you should be OK.
 * This flag is going to be defined in stl_config.h if __STL_USE_NEW_IOSTREAMS is on.
 * But you may wish to force it anyway.
 * Alternatively, you may want to disable it setting __STL_NO_OWN_NAMESPACE on.
 */
// #  define __STL_USE_OWN_NAMESPACE  1
// #  define __STL_NO_OWN_NAMESPACE  1


/* 
 * __STLPORT_NAMESPACE : This is the namespace STLport uses. 
 * Do NOT try to change it to "std".
 * In case you defined __STL_USE_OWN_NAMESPACE, STLport reside there.
 * If you put STLport in std:: (__STL_NO_OWN_NAMESPACE), stlport::
 * namespace is still available and is equivalent to std::
 * STLport also defines user-level macro STLPORT (=__STLPORT_NAMESPACE)
 * which always denotes STLport namespace and is intended to be used in 
 * application's code for portability.
 */

#  define __STLPORT_NAMESPACE stlport

/*
 * If __STL_USE_OWN_NAMESPACE is in effect, STLport will try to rename std:: for the user
 * to stlport::. If you don't want this feature, or if it does not quite work for your
 * compiler, please define the following switch :
 */
// # define __STL_DONT_REDEFINE_STD 1


/*
 * __STL_WHOLE_NATIVE_STD : only meaningful in __STL_USE_OWN_NAMESPACE mode.
 * Normally, STLport only imports necessary components from native std:: namespace -
 * those not yet provided by STLport (<iostream>, <complex>, etc.) 
 * and their dependencies (<string>, <stdexcept>). 
 * You might want everything from std:: being available in std:: namespace when you
 * include corresponding STLport header (like STLport <map> provides std::map as well, etc.),
 * if you are going to use both stlport:: and std:: components in your code.
 * Otherwise this option is not recommended as it increases the size of your object files
 * and slows down compilation.
 */
// # define __STL_WHOLE_NATIVE_STD

/* 
 * __STL_USE_SGI_STRING : Forces use of SGI string even if
 * native <string> is available. Unless you defined __STL_USE_OWN_NAMESPACE,
 * STLport uses native <string> if new iostreams are being used, 
 * as part of compiled runtime library depends on it.
 * You may force use of SGI string uncommenting this macro.
 * IMPORTANT:
 * DO NOT use SGI <string> with native <iostream> unless you recompile 
 * standard C++ runtime library with STLport installed, or
 * (better) defined __STL_USE_OWN_NAMESPACE
 */

// #define  __STL_USE_SGI_STRING  1


/* 
 * Edit relative path below (or put full path) to get native 
 * compiler vendor's headers included. Default is "../include"
 * Hint : never install STLport in the directory that ends with "include"
 */
// #  undef __STL_NATIVE_INCLUDE_PATH
// #  define __STL_NATIVE_INCLUDE_PATH ../include
// same for C library headers like <cstring>
// #  undef __STL_NATIVE_CPP_C_INCLUDE_PATH
// #  define __STL_NATIVE_CPP_C_INCLUDE_PATH ../include
// same for C headers like <string.h>
// #  undef __STL_NATIVE_C_INCLUDE_PATH
// #  define __STL_NATIVE_C_INCLUDE_PATH ../include

/* 
 * Set __STL_DEBUG to turn the "Debug Mode" on.
 * That gets you checked iterators/ranges in the manner
 * of "Safe STL". Very useful for debugging. Thread-safe.
 */
// #define   __STL_DEBUG 1


/* __STL_ASSERTIONS: if defined, then enable runtime checking through the
 * __stl_assert macro.
 */
//#define __STL_ASSERTIONS 1

/*
 * Uncomment this to force all debug diagnostic to be directed through a
 * user-defined global function:
 *	void __stl_debug_message(const char * format_str, ...)
 * instead of predefined STLport routine. 
 * This allows you to take control of debug message output.
 * Default routine calls fprintf(stderr,...)
 * Note : If you set this macro, you must supply __stl_debug_message 
 * function definition somewhere.
 */
//#define __STL_DEBUG_MESSAGE 1

/*
 * Uncomment this to force all failed assertions to be executed through
 * user-defined global function:
 *	void __stl_debug_terminate(void). This allows
 * you to take control of assertion behaviour for debugging purposes.
 * Default routine throws unique exception if __STL_USE_EXCEPTIONS is set,
 * calls abort() otherwise.
 * Note : If you set this macro, you must supply __stl_debug_terminate 
 * function definition somewhere.
 */
//#define __STL_DEBUG_TERMINATE 1

/*
 * Comment this out to enable throwing exceptions from default __stl_debug_terminate()
 * instead of calling abort().
 */
#define __STL_NO_DEBUG_EXCEPTIONS 1

/* 
 * Uncomment that to disable exception handling code 
 */
#define   __STL_NO_EXCEPTIONS 1

/*
 * __STL_NO_NAMESPACES: if defined, don't put the library in namespace
 * stlport:: or std::, even if the compiler supports namespaces
 */

// #define   __STL_NO_NAMESPACES 1

/* 
 * __STL_NO_RELOPS_NAMESPACE: if defined, don't put the relational
 * operator templates (>, <=. >=, !=) in namespace std::rel_ops, even
 * if the compiler supports namespaces and partial ordering of
 * function templates.
 */

// #define __STL_NO_RELOPS_NAMESPACE 1

/* _REENTRANT: define this if your project is multithreaded.
 * STLport uses MT-safe allocator support then. 
*/ 
// #define _REENTRANT

/* 
 * _NOTHREADS: if defined, STLport don't use any 
 * multithreading support.
 */
// #define _NOTHREADS

/* _PTHREADS: if defined, use Posix threads for multithreading support. */
// #define _PTHREADS

/*
 * __STL_NO_NEW_C_HEADERS:  if defined, STLport does not 
 * use new-style headers : <cstdlib>, etc. 
 */
// #define   __STL_NO_NEW_C_HEADERS 1

/* 
 * Uncomment __STL_USE_NEWALLOC to force allocator<T> to use plain "new"
 * instead of SGI optimized node allocator engine.
 */
// #define   __STL_USE_NEWALLOC   1

/* 
 * Uncomment __STL_USE_MALLOC to force allocator<T> to use plain "malloc" 
 * instead of SGI optimized node allocator engine.
 */
// #define   __STL_USE_MALLOC 1

/*
 * Set __STL_DEBUG_ALLOC to use allocators that perform memory debugging,
 * such as padding/checking for memory consistency 
 */
// #define   __STL_DEBUG_ALLOC 1


/*
 * Use this option to catch uninitialized members in your classes.
 * When it is set, construct() and destroy() fill the class storage
 * with __STL_SHRED_BYTE (see below). 
 * Note : __STL_DEBUG and __STL_DEBUG_ALLOC don't set this option automatically.
 */

// # define __STL_DEBUG_UNINITIALIZED 1

/*
 * Uncomment and provide a definition for the byte with which raw memory
 * will be filled if __STL_DEBUG_ALLOC or __STL_DEBUG_UNINITIALIZED is defined. 
 * Choose a value which is likely to cause a noticeable problem if dereferenced 
 * or otherwise abused. A good value may already be defined for your platform; see
 * stl_config.h
 */
// #define __STL_SHRED_BYTE 0xFF

/*
 *  This macro prevents instantiation of at() member function
 *  for containers (vector and deque).
 *  We do not instantiate at() that does not throw range errors -
 *  if this macro is defined, at() method is not defined.
 *
 */
// # define __STL_DONT_THROW_RANGE_ERRORS

//==========================================================
// Compatibility section
//==========================================================

/*
 * __STL_USE_SGI_ALLOCATORS is a hook so that users can disable use of
 * allocator<T> as default parameter for containers, and use SGI
 * raw allocators as default ones, without having to edit library headers.
 * Use of this macro is discouraged.
 */
// #define   __STL_USE_SGI_ALLOCATORS 1

/* 
 * This definition precludes SGI reverse_iterator to be compatible with
 * other parts of MSVC library. (With partial specialization, it just
 * has no effect).
 * Use it _ONLY_ if you use SGI-style reverse_iterator<> template explicitly
 */
// #    define __STL_NO_MSVC50_COMPATIBILITY 1


/* 
 * You should define this macro if compiling with MFC - STLport <stl_config.h>
 * then include <afx.h> instead of <windows.h> to get synchronisation primitives 
 *
 */

// # define __STL_USE_MFC 1

/* 
 * Use abbreviated class names for linker benefit (don't affect interface).
 * This option is obsolete, but should work in this release.
 *
 */
// # define __STL_USE_ABBREVS

/*
 * Use minimum set of default arguments on template classes that have more
 * than one - for example map<>, set<>.
 * This has effect only if __STL_LIMITED_DEFAULT_TEMPLATES is on.
 * If __STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS is set, you'll be able to compile
 * set<T> with those compilers, but you'll have to use __set__<T, less<T>>
 *
 * Affects : map<>, multimap<>, set<>, multiset<>, hash_*<>, 
 * queue<>, priority_queue<>, stack<>, istream_iterator<>
 */

// # define __STL_MINIMUM_DEFAULT_TEMPLATE_PARAMS 1

/*
 * By default, STLport uses proxy technique to enable operator -> for
 * iterators even for those compilers that check the return type of
 * unused instantiations. If this causes problems for you, turn the following
 * switch on to disable proxy ->() operators. This actually should be done
 * in compiler-dependant section, not here. 
 * auto_ptr implements proxy operator even if they are disabled in general,
 * as it is very unlikely that you instantiate auto_ptr<> on pointers and other builtins.
 * However, if this is the case, uncomment second line.
 */

// # define __STL_NO_PROXY_ARROW_OPERATOR 1
// # define __STL_NO_AUTO_PTR_PROXY_ARROW_OPERATOR 1



/*
 * Turn __STL_USE_DECLSPEC on if your project includes multiple DLLs and you want to
 * configure one of them to instantiate STLport exports. 
 * Note : you should definitely do that if you use STLport default node allocator
 * and pass STL objects across DLL boundaries.
 *
 * To do so : you should define __STL_USE_DECLSPEC in all compilation units;
 *                       define __STL_DESIGNATED_DLL when compiling DLL which is designated
 *                       to instantiate STLport exports (other components will import symbols
 *                       from there). This designated DLL should include at least <string> header.
 *  
 * Note : that far, it was only tested with Microsoft Visual C++ compiler.
 */
// # define __STL_USE_DECLSPEC   1
// # define __STL_DESIGNATED_DLL 1


/*
 * Experimental switch for embedded systems where no iostreams are available
 * at all
 */

// # define __STL_NO_IOSTREAMS 1

//==========================================================



// Local Variables:
// mode:C++
// End:
