//
//  bind/bind_template.hpp
//
//  Do not include this header directly.
//
//  Copyright (c) 2001 Peter Dimov and Multi Media Ltd.
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//
//  See http://www.boost.org/libs/bind/bind.html for documentation.
//

    typedef typename result_traits<R, F>::type result_type;

    result_type operator()()
    {
        list0 a;
        BOOST_BIND_EVALUATE;
    }

    result_type operator()() const
    {
        list0 a;
        BOOST_BIND_EVALUATE;
    }

    template<class A1> result_type operator()(A1 & a1)
    {
        list1<A1 &> a(a1);
        BOOST_BIND_EVALUATE;
    }

    template<class A1> result_type operator()(A1 & a1) const
    {
        list1<A1 &> a(a1);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2> result_type operator()(A1 & a1, A2 & a2)
    {
        list2<A1 &, A2 &> a(a1, a2);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2> result_type operator()(A1 & a1, A2 & a2) const
    {
        list2<A1 &, A2 &> a(a1, a2);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3> result_type operator()(A1 & a1, A2 & a2, A3 & a3)
    {
        list3<A1 &, A2 &, A3 &> a(a1, a2, a3);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3> result_type operator()(A1 & a1, A2 & a2, A3 & a3) const
    {
        list3<A1 &, A2 &, A3 &> a(a1, a2, a3);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4)
    {
        list4<A1 &, A2 &, A3 &, A4 &> a(a1, a2, a3, a4);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4) const
    {
        list4<A1 &, A2 &, A3 &, A4 &> a(a1, a2, a3, a4);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5)
    {
        list5<A1 &, A2 &, A3 &, A4 &, A5 &> a(a1, a2, a3, a4, a5);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5) const
    {
        list5<A1 &, A2 &, A3 &, A4 &, A5 &> a(a1, a2, a3, a4, a5);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6)
    {
        list6<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &> a(a1, a2, a3, a4, a5, a6);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6) const
    {
        list6<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &> a(a1, a2, a3, a4, a5, a6);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7)
    {
        list7<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &> a(a1, a2, a3, a4, a5, a6, a7);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7) const
    {
        list7<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &> a(a1, a2, a3, a4, a5, a6, a7);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7, A8 & a8)
    {
        list8<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &, A8 &> a(a1, a2, a3, a4, a5, a6, a7, a8);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7, A8 & a8) const
    {
        list8<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &, A8 &> a(a1, a2, a3, a4, a5, a6, a7, a8);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7, A8 & a8, A9 & a9)
    {
        list9<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &, A8 &, A9 &> a(a1, a2, a3, a4, a5, a6, a7, a8, a9);
        BOOST_BIND_EVALUATE;
    }

    template<class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9> result_type operator()(A1 & a1, A2 & a2, A3 & a3, A4 & a4, A5 & a5, A6 & a6, A7 & a7, A8 & a8, A9 & a9) const
    {
        list9<A1 &, A2 &, A3 &, A4 &, A5 &, A6 &, A7 &, A8 &, A9 &> a(a1, a2, a3, a4, a5, a6, a7, a8, a9);
        BOOST_BIND_EVALUATE;
    }

    template<class A> result_type eval(A & a)
    {
        BOOST_BIND_EVALUATE;
    }

    template<class A> result_type eval(A & a) const
    {
        BOOST_BIND_EVALUATE;
    }

    template<class V> void accept(V & v) const
    {
        BOOST_BIND_VISIT_EACH(v, f_, 0);
        l_.accept(v);
    }

private:

    F f_;
    L l_;
