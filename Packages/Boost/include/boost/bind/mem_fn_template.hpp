//
//  bind/mem_fn_template.hpp
//
//  Do not include this header directly
//
//  Copyright (c) 2001 Peter Dimov and Multi Media Ltd.
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//
//  See http://www.boost.org/libs/bind/mem_fn.html for documentation.
//

// mf0

template<class R, class T BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf0)
{
public:

    typedef R result_type;
    typedef T * argument_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) ())
    F f_;

    template<class U> R call(U & u, T const *) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)();
    }

    template<class U> R call(U & u, void const *) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)();
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf0)(F f): f_(f) {}

    R operator()(T * p) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)();
    }

    template<class U> R operator()(U & u) const
    {
        BOOST_MEM_FN_RETURN call(u, &u);
    }

    R operator()(T & t) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)();
    }
};

// cmf0

template<class R, class T BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf0)
{
public:

    typedef R result_type;
    typedef T const * argument_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) () const)
    F f_;

    template<class U> R call(U & u, T const *) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)();
    }

    template<class U> R call(U & u, void const *) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)();
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf0)(F f): f_(f) {}

    template<class U> R operator()(U const & u) const
    {
        BOOST_MEM_FN_RETURN call(u, &u);
    }

    R operator()(T const & t) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)();
    }
};

// mf1

template<class R, class T, class A1 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf1)
{
public:

    typedef R result_type;
    typedef T * first_argument_type;
    typedef A1 second_argument_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1))
    F f_;

    template<class U, class B1> R call(U & u, T const *, B1 & b1) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1);
    }

    template<class U, class B1> R call(U & u, void const *, B1 & b1) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf1)(F f): f_(f) {}

    R operator()(T * p, A1 a1) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1);
    }

    template<class U> R operator()(U & u, A1 a1) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1);
    }

    R operator()(T & t, A1 a1) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1);
    }
};

// cmf1

template<class R, class T, class A1 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf1)
{
public:

    typedef R result_type;
    typedef T const * first_argument_type;
    typedef A1 second_argument_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1) const)
    F f_;

    template<class U, class B1> R call(U & u, T const *, B1 & b1) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1);
    }

    template<class U, class B1> R call(U & u, void const *, B1 & b1) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf1)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1);
    }

    R operator()(T const & t, A1 a1) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1);
    }
};

// mf2

template<class R, class T, class A1, class A2 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf2)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2))
    F f_;

    template<class U, class B1, class B2> R call(U & u, T const *, B1 & b1, B2 & b2) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2);
    }

    template<class U, class B1, class B2> R call(U & u, void const *, B1 & b1, B2 & b2) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf2)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2);
    }

    R operator()(T & t, A1 a1, A2 a2) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2);
    }
};

// cmf2

template<class R, class T, class A1, class A2 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf2)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2) const)
    F f_;

    template<class U, class B1, class B2> R call(U & u, T const *, B1 & b1, B2 & b2) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2);
    }

    template<class U, class B1, class B2> R call(U & u, void const *, B1 & b1, B2 & b2) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf2)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2);
    }

    R operator()(T const & t, A1 a1, A2 a2) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2);
    }
};

// mf3

template<class R, class T, class A1, class A2, class A3 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf3)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3))
    F f_;

    template<class U, class B1, class B2, class B3> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3);
    }

    template<class U, class B1, class B2, class B3> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf3)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3);
    }
};

// cmf3

template<class R, class T, class A1, class A2, class A3 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf3)
{
public:

    typedef R result_type;

private:

    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3) const)
    F f_;

    template<class U, class B1, class B2, class B3> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3);
    }

    template<class U, class B1, class B2, class B3> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3);
    }

public:

    explicit BOOST_MEM_FN_NAME(cmf3)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3);
    }
};

// mf4

template<class R, class T, class A1, class A2, class A3, class A4 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf4)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4))
    F f_;

    template<class U, class B1, class B2, class B3, class B4> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4);
    }

    template<class U, class B1, class B2, class B3, class B4> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf4)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3, A4 a4) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3, A4 a4) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3, A4 a4) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4);
    }
};

// cmf4

template<class R, class T, class A1, class A2, class A3, class A4 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf4)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4) const)
    F f_;

    template<class U, class B1, class B2, class B3, class B4> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4);
    }

    template<class U, class B1, class B2, class B3, class B4> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf4)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3, A4 a4) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3, A4 a4) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4);
    }
};

// mf5

template<class R, class T, class A1, class A2, class A3, class A4, class A5 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf5)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5))
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf5)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4, a5);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5);
    }
};

// cmf5

template<class R, class T, class A1, class A2, class A3, class A4, class A5 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf5)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5) const)
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf5)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5);
    }
};

// mf6

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf6)
{
public:

    typedef R result_type;

private:

    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6))
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6);
    }

public:

    explicit BOOST_MEM_FN_NAME(mf6)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4, a5, a6);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6);
    }
};

// cmf6

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf6)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6) const)
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf6)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6);
    }
};

// mf7

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf7)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6, A7))
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6, b7);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6, b7);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf7)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4, a5, a6, a7);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6, a7);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6, a7);
    }
};

// cmf7

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf7)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6, A7) const)
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6, b7);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6, b7);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf7)(F f): f_(f) {}

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6, a7);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6, a7);
    }
};

// mf8

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(mf8)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6, A7, A8))
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7, class B8> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7, B8 & b8) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6, b7, b8);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7, class B8> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7, B8 & b8) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6, b7, b8);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(mf8)(F f): f_(f) {}

    R operator()(T * p, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4, a5, a6, a7, a8);
    }

    template<class U> R operator()(U & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6, a7, a8);
    }

    R operator()(T & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6, a7, a8);
    }
};

// cmf8

template<class R, class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8 BOOST_MEM_FN_CLASS_F> class BOOST_MEM_FN_NAME(cmf8)
{
public:

    typedef R result_type;

private:
    
    BOOST_MEM_FN_TYPEDEF(R (BOOST_MEM_FN_CC T::*F) (A1, A2, A3, A4, A5, A6, A7, A8) const)
    F f_;

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7, class B8> R call(U & u, T const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7, B8 & b8) const
    {
        BOOST_MEM_FN_RETURN (u.*f_)(b1, b2, b3, b4, b5, b6, b7, b8);
    }

    template<class U, class B1, class B2, class B3, class B4, class B5, class B6, class B7, class B8> R call(U & u, void const *, B1 & b1, B2 & b2, B3 & b3, B4 & b4, B5 & b5, B6 & b6, B7 & b7, B8 & b8) const
    {
        BOOST_MEM_FN_RETURN (get_pointer(u)->*f_)(b1, b2, b3, b4, b5, b6, b7, b8);
    }

public:
    
    explicit BOOST_MEM_FN_NAME(cmf8)(F f): f_(f) {}

    R operator()(T const * p, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN (p->*f_)(a1, a2, a3, a4, a5, a6, a7, a8);
    }

    template<class U> R operator()(U const & u, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN call(u, &u, a1, a2, a3, a4, a5, a6, a7, a8);
    }

    R operator()(T const & t, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8) const
    {
        BOOST_MEM_FN_RETURN (t.*f_)(a1, a2, a3, a4, a5, a6, a7, a8);
    }
};

