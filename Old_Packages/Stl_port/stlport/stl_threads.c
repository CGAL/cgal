/*
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
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
#ifndef __STL_THREADS_C
#define __STL_THREADS_C

__STL_BEGIN_NAMESPACE

# if defined (__STL_PTHREADS)
# if ( __STL_STATIC_TEMPLATE_DATA > 0 )
    template<int __dummy>
    pthread_mutex_t
    _Swap_lock_struct<__dummy>::_S_swap_lock = PTHREAD_MUTEX_INITIALIZER;
#  else
    __DECLARE_INSTANCE(pthread_mutex_t, _Swap_lock_struct<0>::_S_swap_lock, 
		       =PTHREAD_MUTEX_INITIALIZER);
# endif /* ( __STL_STATIC_TEMPLATE_DATA > 0 ) */
# elif defined (__STL_UITHREADS)
# if ( __STL_STATIC_TEMPLATE_DATA > 0 )
    template<int __dummy>
    mutex_t
    _Swap_lock_struct<__dummy>::_S_swap_lock = DEFAULTMUTEX;
#  else
    __DECLARE_INSTANCE(mutex_t, _Swap_lock_struct<0>::_S_swap_lock, =DEFAULTMUTEX);
# endif /* ( __STL_STATIC_TEMPLATE_DATA > 0 ) */
#endif /* __STL_PTHREADS */

__STL_END_NAMESPACE

#endif /*  __STL_THREADS_C */

// Local Variables:
// mode:C++
// End:
