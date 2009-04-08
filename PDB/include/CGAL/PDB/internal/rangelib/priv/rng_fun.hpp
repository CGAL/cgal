// Boost. Iterable Range Library (rangelib)
//
// Copyright 2003-2004 John Torjo (john@torjo.com) and Matthew Wilson (matthew@synesis.com.au)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.


#ifndef CGAL_PDB_BOOST_RTL_RNGFUN_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_RNGFUN_HPP_INCLUDED


#include <boost/config.hpp>

namespace CGAL { namespace PDB { namespace internal { namespace rangelib { namespace rng {

/* 
    The following will fail to compile, if your compiler does not support partial specialization.

    void add_item_summary( item_summary & summary, sold_items_info & info);
    ....
    rng::for_each( sliced_byeach(v, ptr_fun(add_item_summary), slice_into_days),
        show_item_summary);


   Use this instead:
    rng::for_each( sliced_byeach(v, rng::fun(add_item_summary), slice_into_days),
                                    ^^^^^^^^                   
        show_item_summary);
*/

// note: we need the type of Arg1, without the reference
//
// Example: 'int&' -> 'int'
// (compilers that don't support partial specialization, can't do this)
template<class Arg1, class Arg2, class result>
struct fun_t2_base {
    typedef Arg1 first_argument_type;

    typedef result (*func)(Arg1&,Arg2);
    fun_t2_base( func f) : m_f(f) {}
    result operator()( Arg1 & arg1, Arg2 arg2) const {
        return m_f( arg1, arg2);
    }
    func m_f;
};

template<class Arg1, class Arg2, class Arg3, class result>
struct fun_t3_base {
    typedef Arg1 first_argument_type;

    typedef result (*func)(Arg1&,Arg2,Arg3);
    fun_t3_base( func f) : m_f(f) {}
    result operator()( Arg1 & arg1, Arg2 arg2, Arg3 arg3) const {
        return m_f( arg1, arg2, arg3);
    }
    func m_f;
};



template<class Arg1, class Arg2>
struct void_fun_t2_base {
    typedef Arg1 first_argument_type;

    typedef void (*func)(Arg1&,Arg2);
    void_fun_t2_base( func f) : m_f(f) {}
    void operator()( Arg1 & arg1, Arg2 arg2) const {
        m_f( arg1, arg2);
    }
    func m_f;
};

template<class Arg1, class Arg2, class Arg3>
struct void_fun_t3_base {
    typedef Arg1 first_argument_type;

    typedef void (*func)(Arg1&,Arg2,Arg3);
    void_fun_t3_base( func f) : m_f(f) {}
    void operator()( Arg1 & arg1, Arg2 arg2, Arg3 arg3) const {
        m_f( arg1, arg2, arg3);
    }
    func m_f;
};


namespace detail {
    template< class T> struct fun_finder_impl {
        template< class arg1, class arg2> struct fun2_type {
            typedef fun_t2_base<arg1,arg2,T> result;
        };

        template< class arg1, class arg2, class arg3> struct fun3_type {
            typedef fun_t3_base<arg1,arg2,arg3,T> result;
        };
    };

    template<> struct fun_finder_impl<void> {
        template< class arg1, class arg2> struct fun2_type {
            typedef void_fun_t2_base<arg1,arg2> result;
        };

        template< class arg1, class arg2, class arg3> struct fun3_type {
            typedef void_fun_t3_base<arg1,arg2,arg3> result;
        };
    };

    template<class result, class arg1, class arg2> struct fun_finder2 {
        typedef fun_finder_impl<result> finder_impl_type;
        typedef typename finder_impl_type::template fun2_type<arg1,arg2> finder_type;
        typedef typename finder_type::result type;
    };
    template<class result, class arg1, class arg2, class arg3> struct fun_finder3 {
        typedef fun_finder_impl<result> finder_impl_type;
        typedef typename finder_impl_type::template fun3_type<arg1,arg2,arg3> finder_type;
        typedef typename finder_type::result type;
    };
}

template<class Arg1, class Arg2, class result>
struct fun_t2 : public detail::fun_finder2<result,Arg1,Arg2>::type {
    typedef typename detail::fun_finder2<result,Arg1,Arg2>::type base;
    /*
    painfully discovering once again that VC6 ***** big time

    template< class func> fun_t2( func f) : base(f) {
    }
    */
    
    typedef result (*func)(Arg1&,Arg2);
    fun_t2( func f) : base(f) {
    }
};


template<class Arg1, class Arg2, class Arg3, class result>
struct fun_t3 : public detail::fun_finder3<result,Arg1,Arg2,Arg3>::type {
    typedef typename detail::fun_finder3<result,Arg1,Arg2,Arg3>::type base;
    
    typedef result (*func)(Arg1&,Arg2,Arg3);
    fun_t3( func f) : base(f) {
    }
};



/*
    VC6 users
    (may God have mercy on you;) )

    Because VC6 chokes (amongst many others) on this:

    template<class Arg1, class Arg2> inline void fun( void (pf)(Arg1&, Arg2) ) {}
    void f(int&, long) {}
    int main(int argc, char* argv[]) {
        fun(&f);
    }

  You should use func<reference_type> 
  (specify the reference type).

    The above should become
    int main(int argc, char* argv[]) {
        fun<int>(&f);
    }
*/
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template<class Arg1, class Arg2, class result> inline
fun_t2<Arg1,Arg2,result> fun( result (*pf)(Arg1&, Arg2) ) {
    return fun_t2<Arg1,Arg2,result>(pf);
}
#else
template<class Arg, class Arg1, class Arg2, class result> inline
fun_t2<Arg,Arg2,result> fun( result (*pf)(Arg1, Arg2) ) {
    return fun_t2<Arg,Arg2,result>(pf);
}

#endif

#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template<class Arg1, class Arg2, class Arg3, class result> inline
fun_t3<Arg1,Arg2,Arg3,result> fun( result (*pf)(Arg1&, Arg2, Arg3) ) {
    return fun_t3<Arg1,Arg2,Arg3,result>(pf);
}
#else
template<class Arg, class Arg1, class Arg2, class Arg3, class result> inline
fun_t3<Arg,Arg2,Arg3,result> fun( result (*pf)(Arg1, Arg2, Arg3) ) {
    return fun_t3<Arg,Arg2,Arg3,result>(pf);
}
#endif



}}}}}


#endif

