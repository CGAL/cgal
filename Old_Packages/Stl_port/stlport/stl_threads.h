/*
 * Copyright (c) 1997-1999
 * Silicon Graphics Computer Systems, Inc.
 *
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

// WARNING: This is an internal header file, included by other C++
// standard library headers.  You should not attempt to use this header
// file directly.
// Stl_config.h should be included before this file.


#ifndef __SGI_STL_INTERNAL_THREADS_H
#define __SGI_STL_INTERNAL_THREADS_H

// Supported threading models are native SGI, pthreads, uithreads
// (similar to pthreads, but based on an earlier draft of the Posix
// threads standard), and Win32 threads.  Uithread support by Jochen
// Schlick, 1999, and Solaris threads generalized to them.

#ifndef __STL_CONFIG_H
#include <stl_config.h>
#endif

# if ! defined (__STLPORT_CSTDDEF)
#  include <cstddef>
# endif

#if defined(__STL_SGI_THREADS)

#include <mutex.h>
#include <time.h>

#elif defined(__STL_PTHREADS)

#include <pthread.h>

#elif defined(__STL_WIN32THREADS)

# if !defined (__STLPORT_WINDOWS_H_INCLUDED) && ! defined (_WINDOWS_)

#  ifndef __STL_MSVC
#   define NOMINMAX
#   ifdef __STL_USE_MFC
#    include <afx.h>
#   else
#    include <windows.h>
#   endif

#  else 

// This section serves as a replacement for windows.h header for Visual C++

extern "C" {

#   if (defined(_M_MRX000) || defined(_M_ALPHA) \
       || (defined(_M_PPC) && (_MSC_VER >= 1000))) && !defined(RC_INVOKED)
#    define InterlockedIncrement       _InterlockedIncrement
#    define InterlockedDecrement       _InterlockedDecrement
#    define InterlockedExchange        _InterlockedExchange
#    define __STL_STDCALL
#   else
#    ifdef _MAC
#     define __STL_STDCALL _cdecl
#    else
#     define __STL_STDCALL __stdcall
#    endif
#   endif

__STL_IMPORT_DECLSPEC long __STL_STDCALL InterlockedIncrement(long*);
__STL_IMPORT_DECLSPEC long __STL_STDCALL InterlockedDecrement(long*);
__STL_IMPORT_DECLSPEC long __STL_STDCALL InterlockedExchange(long*, long);

__STL_IMPORT_DECLSPEC void __STL_STDCALL Sleep(unsigned long);

#   if defined (InterlockedIncrement)
#    pragma intrinsic(_InterlockedIncrement)
#    pragma intrinsic(_InterlockedDecrement)
#    pragma intrinsic(_InterlockedExchange)
#   endif
} /* extern "C" */

#  endif /* STL_MSVC */

#  undef min
#  undef max
#  define __STLPORT_WINDOWS_H_INCLUDED
# endif /* _WINDOWS_ */

#elif defined (__STL_UITHREADS)
// this inclusion is potential hazard to bring up all sorts
// of old-style headers. Let's assume vendor already know how
// to deal with that.
#include <ctime>
#include <cstdio>
#include <cwchar>
#include <synch.h>
#endif

__STL_BEGIN_NAMESPACE


// Class _Refcount_Base provides a type, _RC_t, a data member,
// _M_ref_count, and member functions _M_incr and _M_decr, which perform
// atomic preincrement/predecrement.  The constructor initializes 
// _M_ref_count.

// Hack for SGI o32 compilers.
#if defined(__STL_SGI_THREADS) && !defined(__add_and_fetch) && \
    (__mips < 3 || !(defined (_ABIN32) || defined(_ABI64)))
#  define __add_and_fetch(__l,__v) add_then_test((unsigned long*)__l,__v)  
#  define __test_and_set(__l,__v)  test_and_set(__l,__v)
#endif /* o32 */

struct __STL_CLASS_DECLSPEC _Refcount_Base
{
  // The type _RC_t
# ifdef __STL_WIN32THREADS
  typedef long _RC_t;
# else
  typedef size_t _RC_t;
#endif
  
  // The data member _M_ref_count
   volatile _RC_t _M_ref_count;

  // Constructor
# ifdef __STL_PTHREADS
  pthread_mutex_t _M_ref_count_lock;
  _Refcount_Base(_RC_t __n) : _M_ref_count(__n)
    { pthread_mutex_init(&_M_ref_count_lock, 0); }
# elif defined (__STL_UITHREADS)
  mutex_t _M_ref_count_lock;
  _Refcount_Base(_RC_t __n) : _M_ref_count(__n)
  { mutex_init(&_M_ref_count_lock,0,NULL); }
# else
  _Refcount_Base(_RC_t __n) : _M_ref_count(__n) {}
# endif

  // _M_incr and _M_decr
# ifdef __STL_SGI_THREADS
  void _M_incr() {  __add_and_fetch(&_M_ref_count, 1); }
  _RC_t _M_decr() { return __add_and_fetch(&_M_ref_count, (size_t) -1); }
# elif defined (__STL_WIN32THREADS)
   void _M_incr() { InterlockedIncrement((_RC_t*)&_M_ref_count); }
  _RC_t _M_decr() { return InterlockedDecrement((_RC_t*)&_M_ref_count); }
# elif defined(__STL_PTHREADS)
  void _M_incr() {
    pthread_mutex_lock(&_M_ref_count_lock);
    ++_M_ref_count;
    pthread_mutex_unlock(&_M_ref_count_lock);
  }
  _RC_t _M_decr() {
    pthread_mutex_lock(&_M_ref_count_lock);
    volatile _RC_t __tmp = --_M_ref_count;
    pthread_mutex_unlock(&_M_ref_count_lock);
    return __tmp;
  }
# elif defined (__STL_UITHREADS)
  void _M_incr() {
    mutex_lock(&_M_ref_count_lock);
    ++_M_ref_count;
    mutex_unlock(&_M_ref_count_lock);
  }
  _RC_t _M_decr() {
    mutex_lock(&_M_ref_count_lock);
    volatile _RC_t __tmp = --_M_ref_count;
    mutex_unlock(&_M_ref_count_lock);
    return __tmp;
  }
# else  /* No threads */
  void _M_incr() { ++_M_ref_count; }
  _RC_t _M_decr() { return --_M_ref_count; }
# endif
};

// Atomic swap on unsigned long
// This is guaranteed to behave as though it were atomic only if all
// possibly concurrent updates use _Atomic_swap.
// In some cases the operation is emulated with a lock.
# ifdef __STL_SGI_THREADS
inline unsigned long _Atomic_swap(unsigned long * __p, unsigned long __q) {
#   	if __mips < 3 || !(defined (_ABIN32) || defined(_ABI64))
  return test_and_set(__p, __q);
#       else
  return __test_and_set(__p, (unsigned long)__q);
#       endif
}
# elif defined(__STL_WIN32THREADS)
inline unsigned long _Atomic_swap(unsigned long * __p, unsigned long __q) {
  return (unsigned long) InterlockedExchange((long*)__p, (long)__q);
}
# elif defined(__STL_PTHREADS)
// We use a template here only to get a unique initialized instance.
template<int __dummy>
struct _Swap_lock_struct {
  static pthread_mutex_t _S_swap_lock;
};


// This should be portable, but performance is expected
// to be quite awful.  This really needs platform specific
// code.
inline unsigned long _Atomic_swap(unsigned long * __p, unsigned long __q) {
  pthread_mutex_lock(&_Swap_lock_struct<0>::_S_swap_lock);
  unsigned long __result = *__p;
  *__p = __q;
  pthread_mutex_unlock(&_Swap_lock_struct<0>::_S_swap_lock);
  return __result;
}
# elif defined (__STL_UITHREADS)
// any better solutions ?
// We use a template here only to get a unique initialized instance.
template<int __dummy>
struct _Swap_lock_struct {
  static mutex_t _S_swap_lock;
};

// This should be portable, but performance is expected
// to be quite awful.  This really needs platform specific
// code.
inline unsigned long _Atomic_swap(unsigned long * __p, unsigned long __q) {
  mutex_lock(&_Swap_lock_struct<0>::_S_swap_lock);
  unsigned long __result = *__p;
  *__p = __q;
  mutex_unlock(&_Swap_lock_struct<0>::_S_swap_lock);
  return __result;
}
# else
static inline unsigned long _Atomic_swap(unsigned long * __p, unsigned long __q) {
  unsigned long __result = *__p;
  *__p = __q;
  return __result;
}
# endif

// Locking class.  The constructor initializes the lock, the destructor
// destroys it.  There are explicit member functions for acquiring and
// releasing the lock.

// VC++ does not like {} initializers for structs with constructors,
// and, as {} functionality is crucial for default allocator, we
// make this base class just for it. It exports _Init and _Destroy functions to
// _STL_mutex_lock. For non-static cases, clients should use  _STL_mutex_lock.

struct __STL_CLASS_DECLSPEC _STL_mutex_base
{
#if defined(__STL_SGI_THREADS) || defined(__STL_WIN32THREADS)
  // It should be relatively easy to get this to work on any modern Unix.
  volatile unsigned long _M_lock;
  inline void _M_initialize() { _M_lock=0; }
  inline void _M_destroy() {}

  static inline void _S_nsec_sleep(int __log_nsec) {
#     ifdef __STL_SGI_THREADS
          struct timespec __ts;
          /* Max sleep is 2**27nsec ~ 60msec      */
          __ts.tv_sec = 0;
          __ts.tv_nsec = 1 << __log_nsec;
          nanosleep(&__ts, 0);
#     elif defined(__STL_WIN32THREADS)
	  if (__log_nsec <= 20) {
	      Sleep(0);
	  } else {
	      Sleep(1 << (__log_nsec - 20));
	  }
#     else
#	error unimplemented
#     endif
  }

  void _M_acquire_lock() {
    // spins if we suspect uniprocessor
# define __low_spin_max 30
    // spins for multiprocessor
# define __high_spin_max 1000

    static unsigned __spin_max = __low_spin_max;
    static unsigned __last_spins = 0;

    volatile unsigned long* __lock = &this->_M_lock;

    if (_Atomic_swap((unsigned long*)__lock, 1)) {
      unsigned __my_spin_max;
      unsigned __my_last_spins;
      volatile unsigned __junk = 17; 	// Value doesn't matter.
      unsigned  __i;

      __my_spin_max = __spin_max;
      __my_last_spins = __last_spins;
      __junk = 17;
      
      for (__i = 0; __i < __my_spin_max; ++__i) {
	if (__i < __my_last_spins/2 || *__lock) {
	  __junk *= __junk; __junk *= __junk;
	  __junk *= __junk; __junk *= __junk;
	} else {
	  if (!_Atomic_swap((unsigned long*)__lock, 1)) {
	    // got it!
	    // Spinning worked.  Thus we're probably not being scheduled
	    // against the other process with which we were contending.
	    // Thus it makes sense to spin longer the next time.
	    __last_spins = __i;
	    __spin_max = __high_spin_max;
	    return;
	  }
	}
      }

      // We are probably being scheduled against the other process.  Sleep.
      __spin_max = __low_spin_max;

      for (__i = 0 ;; ++__i) {
	int __log_nsec = __i + 6;
	
	if (__log_nsec > 27) __log_nsec = 27;
	if (!_Atomic_swap((unsigned long *)__lock, 1)) {
	  break;
	}
	_S_nsec_sleep(__log_nsec);
      }

    } /* first _Atomic_swap */

  }

  inline void _M_release_lock() {
    volatile unsigned long* __lock = &_M_lock;
#   if defined(__STL_SGI_THREADS) && defined(__GNUC__) && __mips >= 3
        asm("sync");
        *__lock = 0;
#   elif defined(__STL_SGI_THREADS) && __mips >= 3 \
	 && (defined (_ABIN32) || defined(_ABI64))
        __lock_release(__lock);
#   else 
        *__lock = 0;
        // This is not sufficient on many multiprocessors, since
        // writes to protected variables and the lock may be reordered.
#   endif
  }

// We no longer use win32 critical sections.
// They appear to be slower in the contention-free case,
// and they appear difficult to initialize without introducing a race.

#elif defined(__STL_PTHREADS)
  pthread_mutex_t _M_lock;
  inline void _M_initialize() {
	pthread_mutex_init(&_M_lock, NULL);	
  }
  inline void _M_destroy() {
        pthread_mutex_destroy(&_M_lock);
  }
  inline void _M_acquire_lock() { pthread_mutex_lock(&_M_lock); }
  inline void _M_release_lock() { pthread_mutex_unlock(&_M_lock); }

#elif defined (__STL_UITHREADS)
  mutex_t _M_lock;
  inline void _M_initialize() {
    mutex_init(&_M_lock,0,NULL);	
  }
  inline void _M_destroy() {
    mutex_destroy(&_M_lock);
  }
  inline void _M_acquire_lock() { mutex_lock(&_M_lock); }
  inline void _M_release_lock() { mutex_unlock(&_M_lock); }

#else /* No threads */
  inline void _M_initialize() {}
  inline void _M_destroy() {}
  inline void _M_acquire_lock() {}
  inline void _M_release_lock() {}
#endif
};

// Well - behaving class, does not need static initializer
struct __STL_CLASS_DECLSPEC _STL_mutex_lock : public _STL_mutex_base {
  inline _STL_mutex_lock() {
    _M_initialize();
  }

  inline ~_STL_mutex_lock() {
    _M_destroy();
  }
private:
  _STL_mutex_lock(const _STL_mutex_lock&);
  void operator=(const _STL_mutex_lock&);
};

#ifdef __STL_PTHREADS
// Pthreads locks must be statically initialized to something other than
// the default value of zero.
#   define __STL_MUTEX_INITIALIZER = { PTHREAD_MUTEX_INITIALIZER }
#elif defined(__STL_SGI_THREADS) || defined(__STL_WIN32THREADS)
#   define __STL_MUTEX_INITIALIZER = { 0 }
#elif defined (__STL_UITHREADS)
# define   __STL_MUTEX_INITIALIZER = { DEFAULTMUTEX } 
#else
#   define __STL_MUTEX_INITIALIZER
#endif


// A locking class that uses _STL_mutex_lock.  The constructor takes
// a reference to an _STL_mutex_lock, and acquires a lock.  The destructor
// releases the lock.
// It's not clear that this is exactly the right functionality.
// It will probably change in the future.

struct __STL_CLASS_DECLSPEC _STL_auto_lock
{
  _STL_mutex_base& _M_lock;
  
  _STL_auto_lock(_STL_mutex_base& __lock) : _M_lock(__lock)
    { _M_lock._M_acquire_lock(); }
  ~_STL_auto_lock() { _M_lock._M_release_lock(); }

private:
  void operator=(const _STL_auto_lock&);
  _STL_auto_lock(const _STL_auto_lock&);
};

__STL_END_NAMESPACE

# if !defined (__STL_LINK_TIME_INSTANTIATION)
#  include <stl_threads.c>
# endif

#endif /* __SGI_STL_INTERNAL_THREADS_H */

// Local Variables:
// mode:C++
// End:
