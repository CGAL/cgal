/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         regex_synch.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Thread synchronisation for regex code.
  *                Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */

#ifndef BOOST_REGEX_SYNCH_HPP
#define BOOST_REGEX_SYNCH_HPP

#ifndef BOOST_REGEX_CONFIG_HPP
#include <boost/regex/config.hpp>
#endif

#if defined(BOOST_HAS_THREADS)
#  if defined(BOOST_HAS_WINTHREADS)
#     include <windows.h>
#  elif defined(BOOST_HAS_BETHREADS)
#     include <OS.h>
#     include <cassert>
#  elif defined(BOOST_HAS_PTHREADS)
#     include <pthread.h>
#  else
#     error "Unknown threading API"
#  endif
#endif


namespace boost{
   namespace re_detail{

#ifdef __BORLANDC__
   #pragma option push -a8 -b -Vx -Ve -pc
#endif

void BOOST_REGEX_CALL re_init_threads();
void BOOST_REGEX_CALL re_free_threads();

#ifdef BOOST_HAS_THREADS

#  ifdef BOOST_HAS_BETHREADS

typedef sem_id CRITICAL_SECTION;

inline void BOOST_REGEX_CALL InitializeCriticalSection(CRITICAL_SECTION* ps)
{
    *ps = create_sem(1, "regex++");
    assert(*ps > 0);
}

inline void BOOST_REGEX_CALL DeleteCriticalSection(CRITICAL_SECTION* ps)
{
    int t = delete_sem(*ps);
    assert(t == B_NO_ERROR);
}

inline void BOOST_REGEX_CALL EnterCriticalSection(CRITICAL_SECTION* ps)
{
   status_t t = acquire_sem(*ps);
   assert(t == B_NO_ERROR);
}

inline void BOOST_REGEX_CALL LeaveCriticalSection(CRITICAL_SECTION* ps)
{
    status_t t = release_sem(*ps);
    assert(t == B_NO_ERROR);
}

#  elif defined(BOOST_HAS_PTHREADS)

typedef pthread_mutex_t CRITICAL_SECTION;

inline void BOOST_REGEX_CALL InitializeCriticalSection(CRITICAL_SECTION* ps)
{
   pthread_mutex_init(ps, 0);
}

inline void BOOST_REGEX_CALL DeleteCriticalSection(CRITICAL_SECTION* ps)
{
   pthread_mutex_destroy(ps);
}

inline void BOOST_REGEX_CALL EnterCriticalSection(CRITICAL_SECTION* ps)
{
   pthread_mutex_lock(ps);
}

inline void BOOST_REGEX_CALL LeaveCriticalSection(CRITICAL_SECTION* ps)
{
   pthread_mutex_unlock(ps);
}

#  elif !defined(BOOST_HAS_WINTHREADS)
#    error "Unknown threading API"
#  endif

template <class Lock>
class lock_guard
{
   typedef Lock lock_type;
public:
   lock_guard(lock_type& m, bool aq = true)
      : mut(m), owned(false){ acquire(aq); }

   ~lock_guard()
   { acquire(false); }

   void BOOST_REGEX_CALL acquire(bool aq = true)
   {
      if(aq && !owned)
      {
         mut.acquire(true);
         owned = true;
      }
      else if(!aq && owned)
      {
         mut.acquire(false);
         owned = false;
      }
   }
private:
   lock_type& mut;
   bool owned;
   // VC6 warning suppression:
   lock_guard& operator=(const lock_guard&);
};


class critical_section
{
public:
   critical_section()
   { InitializeCriticalSection(&hmutex);}

   critical_section(const critical_section&)
   { InitializeCriticalSection(&hmutex);}

   const critical_section& BOOST_REGEX_CALL operator=(const critical_section&)
   {return *this;}

   ~critical_section()
   {DeleteCriticalSection(&hmutex);}

private:

   void BOOST_REGEX_CALL acquire(bool aq)
   { if(aq) EnterCriticalSection(&hmutex);
      else LeaveCriticalSection(&hmutex);
   }

   CRITICAL_SECTION hmutex;

public:
   typedef lock_guard<critical_section> ro_guard;
   typedef lock_guard<critical_section> rw_guard;

   friend class lock_guard<critical_section>;
};

inline bool BOOST_REGEX_CALL operator==(const critical_section&, const critical_section&)
{
   return false;
}

inline bool BOOST_REGEX_CALL operator<(const critical_section&, const critical_section&)
{
   return true;
}

typedef lock_guard<critical_section> cs_guard;

BOOST_REGEX_DECL extern critical_section* p_re_lock;
BOOST_REGEX_DECL extern unsigned int re_lock_count;

#define BOOST_REGEX_GUARD(inst) boost::re_detail::critical_section::rw_guard g(inst);

#else  // BOOST_HAS_THREADS

#define BOOST_REGEX_GUARD(inst)

#endif // BOOST_HAS_THREADS

#ifdef __BORLANDC__
  #pragma option pop
#endif

} // namespace re_detail
} // namespace boost

#endif // sentry







