#ifndef BOOST_BIND_HPP_INCLUDED
#define BOOST_BIND_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  bind.hpp - binds function objects to arguments
//
//  Copyright (c) 2001-2004 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2001 David Abrahams
//  Copyright (c) 2005 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  See http://www.boost.org/libs/bind/bind.html for documentation.
//

#include <boost/config.hpp>
#include <boost/ref.hpp>
#include <boost/mem_fn.hpp>
#include <boost/type.hpp>
#include <boost/bind/arg.hpp>
#include <boost/detail/workaround.hpp>

// Borland-specific bug, visit_each() silently fails to produce code

#if defined(__BORLANDC__)
#  define BOOST_BIND_VISIT_EACH boost::visit_each
#else
#  define BOOST_BIND_VISIT_EACH visit_each
#endif

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4512) // assignment operator could not be generated
#endif

namespace boost
{

namespace _bi // implementation details
{

// result_traits

template<class R, class F> struct result_traits
{
    typedef R type;
};

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)

struct unspecified {};

template<class F> struct result_traits<unspecified, F>
{
    typedef typename F::result_type type;
};

template<class F> struct result_traits< unspecified, reference_wrapper<F> >
{
    typedef typename F::result_type type;
};

#endif

// ref_compare

template<class T> bool ref_compare(T const & a, T const & b, long)
{
    return a == b;
}

template<class T> bool ref_compare(reference_wrapper<T> const & a, reference_wrapper<T> const & b, int)
{
    return a.get_pointer() == b.get_pointer();
}

// bind_t forward declaration for listN

template<class R, class F, class L> class bind_t;

// value

template<class T> class value
{
public:

    value(T const & t): t_(t) {}

    T & get() { return t_; }
    T const & get() const { return t_; }

    bool operator==(value const & rhs) const
    {
        return t_ == rhs.t_;
    }

private:

    T t_;
};

// type

template<class T> class type {};

// unwrap

template<class F> inline F & unwrap(F * f, long)
{
    return *f;
}

template<class F> inline F & unwrap(reference_wrapper<F> * f, int)
{
    return f->get();
}

template<class F> inline F & unwrap(reference_wrapper<F> const * f, int)
{
    return f->get();
}

#if !( defined(__MWERKS__) && BOOST_WORKAROUND(__MWERKS__, <= 0x3004) )

template<class R, class T> inline _mfi::dm<R, T> unwrap(R T::* * pm, int)
{
    return _mfi::dm<R, T>(*pm);
}

#if !BOOST_WORKAROUND(__IBMCPP__, BOOST_TESTED_AT(600))
// IBM/VisualAge 6.0 is not able to handle this overload.
template<class R, class T> inline _mfi::dm<R, T> unwrap(R T::* const * pm, int)
{
    return _mfi::dm<R, T>(*pm);
}
#endif


#endif

// listN

class list0
{
public:

    list0() {}

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A &, long)
    {
        return unwrap(&f, 0)();
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A &, long) const
    {
        return unwrap(&f, 0)();
    }

    template<class F, class A> void operator()(type<void>, F & f, A &, int)
    {
        unwrap(&f, 0)();
    }

    template<class F, class A> void operator()(type<void>, F const & f, A &, int) const
    {
        unwrap(&f, 0)();
    }

    template<class V> void accept(V &) const
    {
    }

    bool operator==(list0 const &) const
    {
        return true;
    }
};

template<class A1> class list1
{
public:

    explicit list1(A1 a1): a1_(a1) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }

    template<class T> T & operator[] ( _bi::value<T> & v ) const { return v.get(); }

    template<class T> T const & operator[] ( _bi::value<T> const & v ) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
    }

    bool operator==(list1 const & rhs) const
    {
        return ref_compare(a1_, rhs.a1_, 0);
    }

private:

    A1 a1_;
};

template<class A1, class A2> class list2
{
public:

    list2(A1 a1, A2 a2): a1_(a1), a2_(a2) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
    }

    bool operator==(list2 const & rhs) const
    {
        return ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
};

template<class A1, class A2, class A3> class list3
{
public:

    list3(A1 a1, A2 a2, A3 a3): a1_(a1), a2_(a2), a3_(a3) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
    }

    bool operator==(list3 const & rhs) const
    {
        return ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
};

template<class A1, class A2, class A3, class A4> class list4
{
public:

    list4(A1 a1, A2 a2, A3 a3, A4 a4): a1_(a1), a2_(a2), a3_(a3), a4_(a4) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
    }

    bool operator==(list4 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
};

template<class A1, class A2, class A3, class A4, class A5> class list5
{
public:

    list5(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5): a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }
    A5 operator[] (boost::arg<5>) const { return a5_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }
    A5 operator[] (boost::arg<5> (*) ()) const { return a5_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
        BOOST_BIND_VISIT_EACH(v, a5_, 0);
    }

    bool operator==(list5 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0) && ref_compare(a5_, rhs.a5_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
    A5 a5_;
};

template<class A1, class A2, class A3, class A4, class A5, class A6> class list6
{
public:

    list6(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6): a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5), a6_(a6) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }
    A5 operator[] (boost::arg<5>) const { return a5_; }
    A6 operator[] (boost::arg<6>) const { return a6_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }
    A5 operator[] (boost::arg<5> (*) ()) const { return a5_; }
    A6 operator[] (boost::arg<6> (*) ()) const { return a6_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
        BOOST_BIND_VISIT_EACH(v, a5_, 0);
        BOOST_BIND_VISIT_EACH(v, a6_, 0);
    }

    bool operator==(list6 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0) && ref_compare(a5_, rhs.a5_, 0) && ref_compare(a6_, rhs.a6_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
    A5 a5_;
    A6 a6_;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7> class list7
{
public:

    list7(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7): a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5), a6_(a6), a7_(a7) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }
    A5 operator[] (boost::arg<5>) const { return a5_; }
    A6 operator[] (boost::arg<6>) const { return a6_; }
    A7 operator[] (boost::arg<7>) const { return a7_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }
    A5 operator[] (boost::arg<5> (*) ()) const { return a5_; }
    A6 operator[] (boost::arg<6> (*) ()) const { return a6_; }
    A7 operator[] (boost::arg<7> (*) ()) const { return a7_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
        BOOST_BIND_VISIT_EACH(v, a5_, 0);
        BOOST_BIND_VISIT_EACH(v, a6_, 0);
        BOOST_BIND_VISIT_EACH(v, a7_, 0);
    }

    bool operator==(list7 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0) && ref_compare(a5_, rhs.a5_, 0) && ref_compare(a6_, rhs.a6_, 0) &&
            ref_compare(a7_, rhs.a7_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
    A5 a5_;
    A6 a6_;
    A7 a7_;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8> class list8
{
public:

    list8(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8): a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5), a6_(a6), a7_(a7), a8_(a8) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }
    A5 operator[] (boost::arg<5>) const { return a5_; }
    A6 operator[] (boost::arg<6>) const { return a6_; }
    A7 operator[] (boost::arg<7>) const { return a7_; }
    A8 operator[] (boost::arg<8>) const { return a8_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }
    A5 operator[] (boost::arg<5> (*) ()) const { return a5_; }
    A6 operator[] (boost::arg<6> (*) ()) const { return a6_; }
    A7 operator[] (boost::arg<7> (*) ()) const { return a7_; }
    A8 operator[] (boost::arg<8> (*) ()) const { return a8_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
        BOOST_BIND_VISIT_EACH(v, a5_, 0);
        BOOST_BIND_VISIT_EACH(v, a6_, 0);
        BOOST_BIND_VISIT_EACH(v, a7_, 0);
        BOOST_BIND_VISIT_EACH(v, a8_, 0);
    }

    bool operator==(list8 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0) && ref_compare(a5_, rhs.a5_, 0) && ref_compare(a6_, rhs.a6_, 0) &&
            ref_compare(a7_, rhs.a7_, 0) && ref_compare(a8_, rhs.a8_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
    A5 a5_;
    A6 a6_;
    A7 a7_;
    A8 a8_;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9> class list9
{
public:

    list9(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8, A9 a9): a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5), a6_(a6), a7_(a7), a8_(a8), a9_(a9) {}

    A1 operator[] (boost::arg<1>) const { return a1_; }
    A2 operator[] (boost::arg<2>) const { return a2_; }
    A3 operator[] (boost::arg<3>) const { return a3_; }
    A4 operator[] (boost::arg<4>) const { return a4_; }
    A5 operator[] (boost::arg<5>) const { return a5_; }
    A6 operator[] (boost::arg<6>) const { return a6_; }
    A7 operator[] (boost::arg<7>) const { return a7_; }
    A8 operator[] (boost::arg<8>) const { return a8_; }
    A9 operator[] (boost::arg<9>) const { return a9_; }

    A1 operator[] (boost::arg<1> (*) ()) const { return a1_; }
    A2 operator[] (boost::arg<2> (*) ()) const { return a2_; }
    A3 operator[] (boost::arg<3> (*) ()) const { return a3_; }
    A4 operator[] (boost::arg<4> (*) ()) const { return a4_; }
    A5 operator[] (boost::arg<5> (*) ()) const { return a5_; }
    A6 operator[] (boost::arg<6> (*) ()) const { return a6_; }
    A7 operator[] (boost::arg<7> (*) ()) const { return a7_; }
    A8 operator[] (boost::arg<8> (*) ()) const { return a8_; }
    A9 operator[] (boost::arg<9> (*) ()) const { return a9_; }

    template<class T> T & operator[] (_bi::value<T> & v) const { return v.get(); }

    template<class T> T const & operator[] (_bi::value<T> const & v) const { return v.get(); }

    template<class T> T & operator[] (reference_wrapper<T> const & v) const { return v.get(); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> & b) const { return b.eval(*this); }

    template<class R, class F, class L> typename result_traits<R, F>::type operator[] (bind_t<R, F, L> const & b) const { return b.eval(*this); }

    template<class R, class F, class A> R operator()(type<R>, F & f, A & a, long)
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_], a[a9_]);
    }

    template<class R, class F, class A> R operator()(type<R>, F const & f, A & a, long) const
    {
        return unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_], a[a9_]);
    }

    template<class F, class A> void operator()(type<void>, F & f, A & a, int)
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_], a[a9_]);
    }

    template<class F, class A> void operator()(type<void>, F const & f, A & a, int) const
    {
        unwrap(&f, 0)(a[a1_], a[a2_], a[a3_], a[a4_], a[a5_], a[a6_], a[a7_], a[a8_], a[a9_]);
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, a1_, 0);
        BOOST_BIND_VISIT_EACH(v, a2_, 0);
        BOOST_BIND_VISIT_EACH(v, a3_, 0);
        BOOST_BIND_VISIT_EACH(v, a4_, 0);
        BOOST_BIND_VISIT_EACH(v, a5_, 0);
        BOOST_BIND_VISIT_EACH(v, a6_, 0);
        BOOST_BIND_VISIT_EACH(v, a7_, 0);
        BOOST_BIND_VISIT_EACH(v, a8_, 0);
        BOOST_BIND_VISIT_EACH(v, a9_, 0);
    }

    bool operator==(list9 const & rhs) const
    {
        return
            ref_compare(a1_, rhs.a1_, 0) && ref_compare(a2_, rhs.a2_, 0) && ref_compare(a3_, rhs.a3_, 0) &&
            ref_compare(a4_, rhs.a4_, 0) && ref_compare(a5_, rhs.a5_, 0) && ref_compare(a6_, rhs.a6_, 0) &&
            ref_compare(a7_, rhs.a7_, 0) && ref_compare(a8_, rhs.a8_, 0) && ref_compare(a9_, rhs.a9_, 0);
    }

private:

    A1 a1_;
    A2 a2_;
    A3 a3_;
    A4 a4_;
    A5 a5_;
    A6 a6_;
    A7 a7_;
    A8 a8_;
    A9 a9_;
};

// bind_t

#ifndef BOOST_NO_VOID_RETURNS

template<class R, class F, class L> class bind_t
{
public:

    typedef bind_t this_type;

    bind_t(F f, L const & l): f_(f), l_(l) {}

#define BOOST_BIND_RETURN return
#include <boost/bind/bind_template.hpp>
#undef BOOST_BIND_RETURN

};

#else

template<class R> struct bind_t_generator
{

template<class F, class L> class implementation
{
public:

    typedef implementation this_type;

    implementation(F f, L const & l): f_(f), l_(l) {}

#define BOOST_BIND_RETURN return
#include <boost/bind/bind_template.hpp>
#undef BOOST_BIND_RETURN

};

};

template<> struct bind_t_generator<void>
{

template<class F, class L> class implementation
{
private:

    typedef void R;

public:

    typedef implementation this_type;

    implementation(F f, L const & l): f_(f), l_(l) {}

#define BOOST_BIND_RETURN
#include <boost/bind/bind_template.hpp>
#undef BOOST_BIND_RETURN

};

};

template<class R2, class F, class L> class bind_t: public bind_t_generator<R2>::BOOST_NESTED_TEMPLATE implementation<F, L>
{
public:

    bind_t(F f, L const & l): bind_t_generator<R2>::BOOST_NESTED_TEMPLATE implementation<F, L>(f, l) {}

};

#endif

// function_equal

#ifndef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP

// put overloads in _bi, rely on ADL

# ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING

template<class R, class F, class L> bool function_equal( bind_t<R, F, L> const & a, bind_t<R, F, L> const & b )
{
    return a.compare(b);
}

# else

template<class R, class F, class L> bool function_equal_impl( bind_t<R, F, L> const & a, bind_t<R, F, L> const & b, int )
{
    return a.compare(b);
}

# endif // #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#else // BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP

// put overloads in boost

} // namespace _bi

# ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING

template<class R, class F, class L> bool function_equal( _bi::bind_t<R, F, L> const & a, _bi::bind_t<R, F, L> const & b )
{
    return a.compare(b);
}

# else

template<class R, class F, class L> bool function_equal_impl( _bi::bind_t<R, F, L> const & a, _bi::bind_t<R, F, L> const & b, int )
{
    return a.compare(b);
}

# endif // #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING

namespace _bi
{

#endif // BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP

// add_value

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) || (__SUNPRO_CC >= 0x530)

template<class T> struct add_value
{
    typedef _bi::value<T> type;
};

template<class T> struct add_value< value<T> >
{
    typedef _bi::value<T> type;
};

template<class T> struct add_value< reference_wrapper<T> >
{
    typedef reference_wrapper<T> type;
};

template<int I> struct add_value< arg<I> >
{
    typedef boost::arg<I> type;
};

template<int I> struct add_value< arg<I> (*) () >
{
    typedef boost::arg<I> (*type) ();
};

template<class R, class F, class L> struct add_value< bind_t<R, F, L> >
{
    typedef bind_t<R, F, L> type;
};

#else

template<int I> struct _avt_0;

template<> struct _avt_0<1>
{
    template<class T> struct inner
    {
        typedef T type;
    };
};

template<> struct _avt_0<2>
{
    template<class T> struct inner
    {
        typedef value<T> type;
    };
};

typedef char (&_avt_r1) [1];
typedef char (&_avt_r2) [2];

template<class T> _avt_r1 _avt_f(value<T>);
template<class T> _avt_r1 _avt_f(reference_wrapper<T>);
template<int I> _avt_r1 _avt_f(arg<I>);
template<int I> _avt_r1 _avt_f(arg<I> (*) ());
template<class R, class F, class L> _avt_r1 _avt_f(bind_t<R, F, L>);

_avt_r2 _avt_f(...);

template<class T> struct add_value
{
    static T t();
    typedef typename _avt_0<sizeof(_avt_f(t()))>::template inner<T>::type type;
};

#endif

// list_av_N

template<class A1> struct list_av_1
{
    typedef typename add_value<A1>::type B1;
    typedef list1<B1> type;
};

template<class A1, class A2> struct list_av_2
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef list2<B1, B2> type;
};

template<class A1, class A2, class A3> struct list_av_3
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef list3<B1, B2, B3> type;
};

template<class A1, class A2, class A3, class A4> struct list_av_4
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef list4<B1, B2, B3, B4> type;
};

template<class A1, class A2, class A3, class A4, class A5> struct list_av_5
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef typename add_value<A5>::type B5;
    typedef list5<B1, B2, B3, B4, B5> type;
};

template<class A1, class A2, class A3, class A4, class A5, class A6> struct list_av_6
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef typename add_value<A5>::type B5;
    typedef typename add_value<A6>::type B6;
    typedef list6<B1, B2, B3, B4, B5, B6> type;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7> struct list_av_7
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef typename add_value<A5>::type B5;
    typedef typename add_value<A6>::type B6;
    typedef typename add_value<A7>::type B7;
    typedef list7<B1, B2, B3, B4, B5, B6, B7> type;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8> struct list_av_8
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef typename add_value<A5>::type B5;
    typedef typename add_value<A6>::type B6;
    typedef typename add_value<A7>::type B7;
    typedef typename add_value<A8>::type B8;
    typedef list8<B1, B2, B3, B4, B5, B6, B7, B8> type;
};

template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9> struct list_av_9
{
    typedef typename add_value<A1>::type B1;
    typedef typename add_value<A2>::type B2;
    typedef typename add_value<A3>::type B3;
    typedef typename add_value<A4>::type B4;
    typedef typename add_value<A5>::type B5;
    typedef typename add_value<A6>::type B6;
    typedef typename add_value<A7>::type B7;
    typedef typename add_value<A8>::type B8;
    typedef typename add_value<A9>::type B9;
    typedef list9<B1, B2, B3, B4, B5, B6, B7, B8, B9> type;
};

// operator!

struct logical_not
{
    template<class V> bool operator()(V const & v) const { return !v; }
};

template<class R, class F, class L>
    bind_t< bool, logical_not, list1< bind_t<R, F, L> > >
    operator! (bind_t<R, F, L> const & f)
{
    typedef list1< bind_t<R, F, L> > list_type;
    return bind_t<bool, logical_not, list_type> ( logical_not(), list_type(f) );
}

// relational operators

#define BOOST_BIND_OPERATOR( op, name ) \
\
struct name \
{ \
    template<class V, class W> bool operator()(V const & v, W const & w) const { return v op w; } \
}; \
 \
template<class R, class F, class L, class A2> \
    bind_t< bool, name, list2< bind_t<R, F, L>, typename add_value<A2>::type > > \
    operator op (bind_t<R, F, L> const & f, A2 a2) \
{ \
    typedef typename add_value<A2>::type B2; \
    typedef list2< bind_t<R, F, L>, B2> list_type; \
    return bind_t<bool, name, list_type> ( name(), list_type(f, a2) ); \
}

BOOST_BIND_OPERATOR( ==, equal )
BOOST_BIND_OPERATOR( !=, not_equal )

BOOST_BIND_OPERATOR( <, less )
BOOST_BIND_OPERATOR( <=, less_equal )

BOOST_BIND_OPERATOR( >, greater )
BOOST_BIND_OPERATOR( >=, greater_equal )

#undef BOOST_BIND_OPERATOR

#if defined(__GNUC__) && BOOST_WORKAROUND(__GNUC__, < 3)

// resolve ambiguity with rel_ops

#define BOOST_BIND_OPERATOR( op, name ) \
\
template<class R, class F, class L> \
    bind_t< bool, name, list2< bind_t<R, F, L>, bind_t<R, F, L> > > \
    operator op (bind_t<R, F, L> const & f, bind_t<R, F, L> const & g) \
{ \
    typedef list2< bind_t<R, F, L>, bind_t<R, F, L> > list_type; \
    return bind_t<bool, name, list_type> ( name(), list_type(f, g) ); \
}

BOOST_BIND_OPERATOR( !=, not_equal )
BOOST_BIND_OPERATOR( <=, less_equal )
BOOST_BIND_OPERATOR( >, greater )
BOOST_BIND_OPERATOR( >=, greater_equal )

#endif

} // namespace _bi

// visit_each

template<class V, class T> void visit_each(V & v, _bi::value<T> const & t, int)
{
    BOOST_BIND_VISIT_EACH(v, t.get(), 0);
}

template<class V, class R, class F, class L> void visit_each(V & v, _bi::bind_t<R, F, L> const & t, int)
{
    t.accept(v);
}

// bind

#ifndef BOOST_BIND
#define BOOST_BIND bind
#endif

// generic function objects

template<class R, class F>
    _bi::bind_t<R, F, _bi::list0>
    BOOST_BIND(F f)
{
    typedef _bi::list0 list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type());
}

template<class R, class F, class A1>
    _bi::bind_t<R, F, typename _bi::list_av_1<A1>::type>
    BOOST_BIND(F f, A1 a1)
{
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type(a1));
}

template<class R, class F, class A1, class A2>
    _bi::bind_t<R, F, typename _bi::list_av_2<A1, A2>::type>
    BOOST_BIND(F f, A1 a1, A2 a2)
{
    typedef typename _bi::list_av_2<A1, A2>::type list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type(a1, a2));
}

template<class R, class F, class A1, class A2, class A3>
    _bi::bind_t<R, F, typename _bi::list_av_3<A1, A2, A3>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3)
{
    typedef typename _bi::list_av_3<A1, A2, A3>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3));
}

template<class R, class F, class A1, class A2, class A3, class A4>
    _bi::bind_t<R, F, typename _bi::list_av_4<A1, A2, A3, A4>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4)
{
    typedef typename _bi::list_av_4<A1, A2, A3, A4>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5>
    _bi::bind_t<R, F, typename _bi::list_av_5<A1, A2, A3, A4, A5>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5)
{
    typedef typename _bi::list_av_5<A1, A2, A3, A4, A5>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6>
    _bi::bind_t<R, F, typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6)
{
    typedef typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    _bi::bind_t<R, F, typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7)
{
    typedef typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
    _bi::bind_t<R, F, typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8)
{
    typedef typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
    _bi::bind_t<R, F, typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8, A9 a9)
{
    typedef typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8, a9));
}

// generic function objects, alternative syntax

template<class R, class F>
    _bi::bind_t<R, F, _bi::list0>
    BOOST_BIND(boost::type<R>, F f)
{
    typedef _bi::list0 list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type());
}

template<class R, class F, class A1>
    _bi::bind_t<R, F, typename _bi::list_av_1<A1>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1)
{
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type(a1));
}

template<class R, class F, class A1, class A2>
    _bi::bind_t<R, F, typename _bi::list_av_2<A1, A2>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2)
{
    typedef typename _bi::list_av_2<A1, A2>::type list_type;
    return _bi::bind_t<R, F, list_type> (f, list_type(a1, a2));
}

template<class R, class F, class A1, class A2, class A3>
    _bi::bind_t<R, F, typename _bi::list_av_3<A1, A2, A3>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3)
{
    typedef typename _bi::list_av_3<A1, A2, A3>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3));
}

template<class R, class F, class A1, class A2, class A3, class A4>
    _bi::bind_t<R, F, typename _bi::list_av_4<A1, A2, A3, A4>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4)
{
    typedef typename _bi::list_av_4<A1, A2, A3, A4>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5>
    _bi::bind_t<R, F, typename _bi::list_av_5<A1, A2, A3, A4, A5>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5)
{
    typedef typename _bi::list_av_5<A1, A2, A3, A4, A5>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6>
    _bi::bind_t<R, F, typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6)
{
    typedef typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    _bi::bind_t<R, F, typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7)
{
    typedef typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
    _bi::bind_t<R, F, typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8)
{
    typedef typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8));
}

template<class R, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
    _bi::bind_t<R, F, typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type>
    BOOST_BIND(boost::type<R>, F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8, A9 a9)
{
    typedef typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type list_type;
    return _bi::bind_t<R, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8, a9));
}

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)

// adaptable function objects

template<class F>
    _bi::bind_t<_bi::unspecified, F, _bi::list0>
    BOOST_BIND(F f)
{
    typedef _bi::list0 list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type> (f, list_type());
}

template<class F, class A1>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_1<A1>::type>
    BOOST_BIND(F f, A1 a1)
{
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type> (f, list_type(a1));
}

template<class F, class A1, class A2>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_2<A1, A2>::type>
    BOOST_BIND(F f, A1 a1, A2 a2)
{
    typedef typename _bi::list_av_2<A1, A2>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type> (f, list_type(a1, a2));
}

template<class F, class A1, class A2, class A3>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_3<A1, A2, A3>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3)
{
    typedef typename _bi::list_av_3<A1, A2, A3>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3));
}

template<class F, class A1, class A2, class A3, class A4>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_4<A1, A2, A3, A4>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4)
{
    typedef typename _bi::list_av_4<A1, A2, A3, A4>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4));
}

template<class F, class A1, class A2, class A3, class A4, class A5>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_5<A1, A2, A3, A4, A5>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5)
{
    typedef typename _bi::list_av_5<A1, A2, A3, A4, A5>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4, a5));
}

template<class F, class A1, class A2, class A3, class A4, class A5, class A6>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6)
{
    typedef typename _bi::list_av_6<A1, A2, A3, A4, A5, A6>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6));
}

template<class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7)
{
    typedef typename _bi::list_av_7<A1, A2, A3, A4, A5, A6, A7>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7));
}

template<class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8)
{
    typedef typename _bi::list_av_8<A1, A2, A3, A4, A5, A6, A7, A8>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8));
}

template<class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
    _bi::bind_t<_bi::unspecified, F, typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type>
    BOOST_BIND(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8, A9 a9)
{
    typedef typename _bi::list_av_9<A1, A2, A3, A4, A5, A6, A7, A8, A9>::type list_type;
    return _bi::bind_t<_bi::unspecified, F, list_type>(f, list_type(a1, a2, a3, a4, a5, a6, a7, a8, a9));
}

#endif // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)

// function pointers

#define BOOST_BIND_CC
#define BOOST_BIND_ST

#include <boost/bind/bind_cc.hpp>

#undef BOOST_BIND_CC
#undef BOOST_BIND_ST

#ifdef BOOST_BIND_ENABLE_STDCALL

#define BOOST_BIND_CC __stdcall
#define BOOST_BIND_ST

#include <boost/bind/bind_cc.hpp>

#undef BOOST_BIND_CC
#undef BOOST_BIND_ST

#endif

#ifdef BOOST_BIND_ENABLE_FASTCALL

#define BOOST_BIND_CC __fastcall
#define BOOST_BIND_ST

#include <boost/bind/bind_cc.hpp>

#undef BOOST_BIND_CC
#undef BOOST_BIND_ST

#endif

#ifdef BOOST_BIND_ENABLE_PASCAL

#define BOOST_BIND_ST pascal
#define BOOST_BIND_CC

#include <boost/bind/bind_cc.hpp>

#undef BOOST_BIND_ST
#undef BOOST_BIND_CC

#endif

// member function pointers

#define BOOST_BIND_MF_NAME(X) X
#define BOOST_BIND_MF_CC

#include <boost/bind/bind_mf_cc.hpp>

#undef BOOST_BIND_MF_NAME
#undef BOOST_BIND_MF_CC

#ifdef BOOST_MEM_FN_ENABLE_CDECL

#define BOOST_BIND_MF_NAME(X) X##_cdecl
#define BOOST_BIND_MF_CC __cdecl

#include <boost/bind/bind_mf_cc.hpp>

#undef BOOST_BIND_MF_NAME
#undef BOOST_BIND_MF_CC

#endif

#ifdef BOOST_MEM_FN_ENABLE_STDCALL

#define BOOST_BIND_MF_NAME(X) X##_stdcall
#define BOOST_BIND_MF_CC __stdcall

#include <boost/bind/bind_mf_cc.hpp>

#undef BOOST_BIND_MF_NAME
#undef BOOST_BIND_MF_CC

#endif

#ifdef BOOST_MEM_FN_ENABLE_FASTCALL

#define BOOST_BIND_MF_NAME(X) X##_fastcall
#define BOOST_BIND_MF_CC __fastcall

#include <boost/bind/bind_mf_cc.hpp>

#undef BOOST_BIND_MF_NAME
#undef BOOST_BIND_MF_CC

#endif

// data member pointers

/*

#if defined(__GNUC__) && (__GNUC__ == 2)

namespace _bi
{

template<class T> struct add_cref
{
    typedef T const & type;
};

template<class T> struct add_cref< T & >
{
    typedef T const & type;
};

template<> struct add_cref<void>
{
    typedef void type;
};

} // namespace _bi

template<class R, class T, class A1>
_bi::bind_t< typename _bi::add_cref<R>::type, _mfi::dm<R, T>, typename _bi::list_av_1<A1>::type >
    BOOST_BIND(R T::*f, A1 a1)
{
    typedef _mfi::dm<R, T> F;
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<typename _bi::add_cref<R>::type, F, list_type>(F(f), list_type(a1));
}

#else

template<class R, class T, class A1>
_bi::bind_t< R const &, _mfi::dm<R, T>, typename _bi::list_av_1<A1>::type >
    BOOST_BIND(R T::*f, A1 a1)
{
    typedef _mfi::dm<R, T> F;
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<R const &, F, list_type>(F(f), list_type(a1));
}

#endif

*/

template<class R, class T, class A1>
_bi::bind_t< R, _mfi::dm<R, T>, typename _bi::list_av_1<A1>::type >
    BOOST_BIND(R T::*f, A1 a1)
{
    typedef _mfi::dm<R, T> F;
    typedef typename _bi::list_av_1<A1>::type list_type;
    return _bi::bind_t<R, F, list_type>( F(f), list_type(a1) );
}

} // namespace boost

#ifndef BOOST_BIND_NO_PLACEHOLDERS

# include <boost/bind/placeholders.hpp>

#endif

#ifdef BOOST_MSVC
# pragma warning(default: 4512) // assignment operator could not be generated
# pragma warning(pop)
#endif

#endif // #ifndef BOOST_BIND_HPP_INCLUDED
