// preprocessed version of 'boost/mpl/vector/vector30.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20
    >
struct vector21
{
    typedef aux::vector_tag<21> tag;
    typedef vector21 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    

    typedef void_ item21;
    typedef T20 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,21> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 20> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector21<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector20<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<21>
{
    template< typename V > struct result_
    {
        typedef typename V::item21 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 21> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,21 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 21> >
    : size_traits< aux::vector_tag< 21> >
{
};

template<>
struct clear_traits< aux::vector_tag< 21> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21
    >
struct vector22
{
    typedef aux::vector_tag<22> tag;
    typedef vector22 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    

    typedef void_ item22;
    typedef T21 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,22> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 21> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector22<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector21<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<22>
{
    template< typename V > struct result_
    {
        typedef typename V::item22 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 22> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,22 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 22> >
    : size_traits< aux::vector_tag< 22> >
{
};

template<>
struct clear_traits< aux::vector_tag< 22> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22
    >
struct vector23
{
    typedef aux::vector_tag<23> tag;
    typedef vector23 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    

    typedef void_ item23;
    typedef T22 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,23> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 22> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector23<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector22<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<23>
{
    template< typename V > struct result_
    {
        typedef typename V::item23 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 23> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,23 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 23> >
    : size_traits< aux::vector_tag< 23> >
{
};

template<>
struct clear_traits< aux::vector_tag< 23> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23
    >
struct vector24
{
    typedef aux::vector_tag<24> tag;
    typedef vector24 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    

    typedef void_ item24;
    typedef T23 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,24> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 23> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector24<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector23<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<24>
{
    template< typename V > struct result_
    {
        typedef typename V::item24 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 24> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,24 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 24> >
    : size_traits< aux::vector_tag< 24> >
{
};

template<>
struct clear_traits< aux::vector_tag< 24> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    >
struct vector25
{
    typedef aux::vector_tag<25> tag;
    typedef vector25 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    

    typedef void_ item25;
    typedef T24 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,25> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 24> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector25<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector24<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<25>
{
    template< typename V > struct result_
    {
        typedef typename V::item25 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 25> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,25 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 25> >
    : size_traits< aux::vector_tag< 25> >
{
};

template<>
struct clear_traits< aux::vector_tag< 25> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    , typename T25
    >
struct vector26
{
    typedef aux::vector_tag<26> tag;
    typedef vector26 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    typedef T25 item25;
    

    typedef void_ item26;
    typedef T25 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,26> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 25> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector26<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            , typename Vector::item24
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector25<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            , typename Vector::item25
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<26>
{
    template< typename V > struct result_
    {
        typedef typename V::item26 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 26> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,26 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 26> >
    : size_traits< aux::vector_tag< 26> >
{
};

template<>
struct clear_traits< aux::vector_tag< 26> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    , typename T25, typename T26
    >
struct vector27
{
    typedef aux::vector_tag<27> tag;
    typedef vector27 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    typedef T25 item25;
    typedef T26 item26;
    

    typedef void_ item27;
    typedef T26 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,27> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 26> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector27<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            , typename Vector::item24, typename Vector::item25
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector26<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            , typename Vector::item25, typename Vector::item26
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<27>
{
    template< typename V > struct result_
    {
        typedef typename V::item27 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 27> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,27 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 27> >
    : size_traits< aux::vector_tag< 27> >
{
};

template<>
struct clear_traits< aux::vector_tag< 27> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    , typename T25, typename T26, typename T27
    >
struct vector28
{
    typedef aux::vector_tag<28> tag;
    typedef vector28 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    typedef T25 item25;
    typedef T26 item26;
    typedef T27 item27;
    

    typedef void_ item28;
    typedef T27 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,28> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 27> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector28<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            , typename Vector::item24, typename Vector::item25
            , typename Vector::item26
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector27<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            , typename Vector::item25, typename Vector::item26
            , typename Vector::item27
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<28>
{
    template< typename V > struct result_
    {
        typedef typename V::item28 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 28> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,28 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 28> >
    : size_traits< aux::vector_tag< 28> >
{
};

template<>
struct clear_traits< aux::vector_tag< 28> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    , typename T25, typename T26, typename T27, typename T28
    >
struct vector29
{
    typedef aux::vector_tag<29> tag;
    typedef vector29 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    typedef T25 item25;
    typedef T26 item26;
    typedef T27 item27;
    typedef T28 item28;
    

    typedef void_ item29;
    typedef T28 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,29> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 28> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector29<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            , typename Vector::item24, typename Vector::item25
            , typename Vector::item26, typename Vector::item27
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector28<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            , typename Vector::item25, typename Vector::item26
            , typename Vector::item27, typename Vector::item28
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<29>
{
    template< typename V > struct result_
    {
        typedef typename V::item29 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 29> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,29 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 29> >
    : size_traits< aux::vector_tag< 29> >
{
};

template<>
struct clear_traits< aux::vector_tag< 29> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    , typename T10, typename T11, typename T12, typename T13, typename T14
    , typename T15, typename T16, typename T17, typename T18, typename T19
    , typename T20, typename T21, typename T22, typename T23, typename T24
    , typename T25, typename T26, typename T27, typename T28, typename T29
    >
struct vector30
{
    typedef aux::vector_tag<30> tag;
    typedef vector30 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    typedef T9 item9;
    typedef T10 item10;
    typedef T11 item11;
    typedef T12 item12;
    typedef T13 item13;
    typedef T14 item14;
    typedef T15 item15;
    typedef T16 item16;
    typedef T17 item17;
    typedef T18 item18;
    typedef T19 item19;
    typedef T20 item20;
    typedef T21 item21;
    typedef T22 item22;
    typedef T23 item23;
    typedef T24 item24;
    typedef T25 item25;
    typedef T26 item26;
    typedef T27 item27;
    typedef T28 item28;
    typedef T29 item29;
    

    typedef void_ item30;
    typedef T29 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,30> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 29> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector30<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8, typename Vector::item9
            , typename Vector::item10, typename Vector::item11
            , typename Vector::item12, typename Vector::item13
            , typename Vector::item14, typename Vector::item15
            , typename Vector::item16, typename Vector::item17
            , typename Vector::item18, typename Vector::item19
            , typename Vector::item20, typename Vector::item21
            , typename Vector::item22, typename Vector::item23
            , typename Vector::item24, typename Vector::item25
            , typename Vector::item26, typename Vector::item27
            , typename Vector::item28
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector29<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9, typename Vector::item10
            , typename Vector::item11, typename Vector::item12
            , typename Vector::item13, typename Vector::item14
            , typename Vector::item15, typename Vector::item16
            , typename Vector::item17, typename Vector::item18
            , typename Vector::item19, typename Vector::item20
            , typename Vector::item21, typename Vector::item22
            , typename Vector::item23, typename Vector::item24
            , typename Vector::item25, typename Vector::item26
            , typename Vector::item27, typename Vector::item28
            , typename Vector::item29
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<30>
{
    template< typename V > struct result_
    {
        typedef typename V::item30 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 30> >
{
    template< typename V, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,30 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 30> >
    : size_traits< aux::vector_tag< 30> >
{
};

template<>
struct clear_traits< aux::vector_tag< 30> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

} // namespace mpl
} // namespace boost

