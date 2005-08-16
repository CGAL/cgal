#ifndef BOOST_DETAIL_SP_COUNTED_BASE_GCC_IA64_HPP_INCLUDED
#define BOOST_DETAIL_SP_COUNTED_BASE_GCC_IA64_HPP_INCLUDED

//
//  detail/sp_counted_base_gcc_ia64.hpp - g++ on IA64
//
//  Copyright (c) 2001, 2002, 2003 Peter Dimov and Multi Media Ltd.
//  Copyright 2004-2005 Peter Dimov
//  Copyright 2005 Ben Hutchings
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//
//  Lock-free algorithm by Alexander Terekhov
//

#include <typeinfo>

namespace boost
{

namespace detail
{

inline void atomic_increment( long * pw )
{
    // ++*pw;

    long tmp;

    // No barrier is required here but fetchadd always has an acquire or
    // release barrier associated with it.  We choose release as it should be
    // cheaper.
    __asm__ ("fetchadd8.rel %0=[%2],1" :
         "=r"(tmp), "=m"(*pw) :
         "r"(pw));
}

inline long atomic_decrement( long * pw )
{
    // return --*pw;

    long rv;

    __asm__ ("     fetchadd8.rel %0=[%2],-1 ;; \n"
             "     cmp.eq        p7,p0=1,%0 ;; \n"
             "(p7) ld8.acq       %0=[%2]    " :
             "=&r"(rv), "=m"(*pw) :
             "r"(pw) :
             "p7");

    return rv;
}

inline long atomic_conditional_increment( long * pw )
{
    // if( *pw != 0 ) ++*pw;
    // return *pw;

    long rv, tmp, tmp2;

    __asm__ ("0:   ld8          %0=[%4]           ;; \n"
         "     cmp.eq       p7,p0=0,%0        ;; \n"
         "(p7) br.cond.spnt 1f                \n"
         "     mov          ar.ccv=%0         \n"
         "     add          %1=1,%0           ;; \n"
         "     cmpxchg8.acq %2=[%4],%1,ar.ccv ;; \n"
         "     cmp.ne       p7,p0=%0,%2       ;; \n"
         "(p7) br.cond.spnt 0b                \n"
         "     mov          %0=%1             ;; \n"
         "1:" : 
         "=&r"(rv), "=&r"(tmp), "=&r"(tmp2), "=m"(*pw) :
         "r"(pw) :
         "ar.ccv", "p7");

    return rv;
}

class sp_counted_base
{
private:

    sp_counted_base( sp_counted_base const & );
    sp_counted_base & operator= ( sp_counted_base const & );

    long use_count_;        // #shared
    long weak_count_;       // #weak + (#shared != 0)

public:

    sp_counted_base(): use_count_( 1 ), weak_count_( 1 )
    {
    }

    virtual ~sp_counted_base() // nothrow
    {
    }

    // dispose() is called when use_count_ drops to zero, to release
    // the resources managed by *this.

    virtual void dispose() = 0; // nothrow

    // destroy() is called when weak_count_ drops to zero.

    virtual void destroy() // nothrow
    {
        delete this;
    }

    virtual void * get_deleter( std::type_info const & ti ) = 0;

    void add_ref_copy()
    {
        atomic_increment( &use_count_ );
    }

    bool add_ref_lock() // true on success
    {
        return atomic_conditional_increment( &use_count_ ) != 0;
    }

    void release() // nothrow
    {
        if( atomic_decrement( &use_count_ ) == 0 )
        {
            dispose();
            weak_release();
        }
    }

    void weak_add_ref() // nothrow
    {
        atomic_increment( &weak_count_ );
    }

    void weak_release() // nothrow
    {
        if( atomic_decrement( &weak_count_ ) == 0 )
        {
            destroy();
        }
    }

    long use_count() const // nothrow
    {
        return static_cast<long const volatile &>( use_count_ );
    }
};

} // namespace detail

} // namespace boost

#endif  // #ifndef BOOST_DETAIL_SP_COUNTED_BASE_GCC_IA64_HPP_INCLUDED
