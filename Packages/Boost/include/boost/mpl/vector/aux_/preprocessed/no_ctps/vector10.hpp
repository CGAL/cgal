// preprocessed version of 'boost/mpl/vector/vector10.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

namespace aux {
template<> struct vector_item_impl<0>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item0 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 0> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct size_traits< aux::vector_tag< 0> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,0 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 0> >
    : size_traits< aux::vector_tag< 0> >
{
};

template<>
struct clear_traits< aux::vector_tag< 0> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0
    >
struct vector1
{
    typedef aux::vector_tag<1> tag;
    typedef vector1 type;
    typedef T0 item0;
    typedef void_ item1;
    typedef T0 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,1> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 0> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector1<
              T
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<
              
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<1>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item1 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 1> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,1 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 1> >
    : size_traits< aux::vector_tag< 1> >
{
};

template<>
struct clear_traits< aux::vector_tag< 1> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1
    >
struct vector2
{
    typedef aux::vector_tag<2> tag;
    typedef vector2 type;
    typedef T0 item0;
    typedef T1 item1;
    

    typedef void_ item2;
    typedef T1 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,2> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 1> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector2<
              T
              ,
              typename Vector::item0
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector1<
              typename Vector::item1
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<2>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item2 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 2> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,2 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 2> >
    : size_traits< aux::vector_tag< 2> >
{
};

template<>
struct clear_traits< aux::vector_tag< 2> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2
    >
struct vector3
{
    typedef aux::vector_tag<3> tag;
    typedef vector3 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    

    typedef void_ item3;
    typedef T2 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,3> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 2> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector3<
              T
              ,
              typename Vector::item0, typename Vector::item1
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector2<
              typename Vector::item1, typename Vector::item2
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<3>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item3 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 3> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,3 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 3> >
    : size_traits< aux::vector_tag< 3> >
{
};

template<>
struct clear_traits< aux::vector_tag< 3> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3
    >
struct vector4
{
    typedef aux::vector_tag<4> tag;
    typedef vector4 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    

    typedef void_ item4;
    typedef T3 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,4> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 3> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector4<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector3<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<4>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item4 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 4> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,4 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 4> >
    : size_traits< aux::vector_tag< 4> >
{
};

template<>
struct clear_traits< aux::vector_tag< 4> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    >
struct vector5
{
    typedef aux::vector_tag<5> tag;
    typedef vector5 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    

    typedef void_ item5;
    typedef T4 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,5> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 4> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector5<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector4<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<5>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item5 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 5> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,5 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 5> >
    : size_traits< aux::vector_tag< 5> >
{
};

template<>
struct clear_traits< aux::vector_tag< 5> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct vector6
{
    typedef aux::vector_tag<6> tag;
    typedef vector6 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    

    typedef void_ item6;
    typedef T5 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,6> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 5> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector6<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector5<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<6>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item6 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 6> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,6 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 6> >
    : size_traits< aux::vector_tag< 6> >
{
};

template<>
struct clear_traits< aux::vector_tag< 6> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6
    >
struct vector7
{
    typedef aux::vector_tag<7> tag;
    typedef vector7 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    

    typedef void_ item7;
    typedef T6 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,7> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 6> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector7<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector6<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<7>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item7 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 7> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,7 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 7> >
    : size_traits< aux::vector_tag< 7> >
{
};

template<>
struct clear_traits< aux::vector_tag< 7> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7
    >
struct vector8
{
    typedef aux::vector_tag<8> tag;
    typedef vector8 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    

    typedef void_ item8;
    typedef T7 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,8> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 7> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector8<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector7<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<8>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item8 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 8> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,8 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 8> >
    : size_traits< aux::vector_tag< 8> >
{
};

template<>
struct clear_traits< aux::vector_tag< 8> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8
    >
struct vector9
{
    typedef aux::vector_tag<9> tag;
    typedef vector9 type;
    typedef T0 item0;
    typedef T1 item1;
    typedef T2 item2;
    typedef T3 item3;
    typedef T4 item4;
    typedef T5 item5;
    typedef T6 item6;
    typedef T7 item7;
    typedef T8 item8;
    

    typedef void_ item9;
    typedef T8 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,9> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 8> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector9<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector8<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<9>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item9 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 9> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,9 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 9> >
    : size_traits< aux::vector_tag< 9> >
{
};

template<>
struct clear_traits< aux::vector_tag< 9> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    >
struct vector10
{
    typedef aux::vector_tag<10> tag;
    typedef vector10 type;
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
    

    typedef void_ item10;
    typedef T9 back;
    typedef vector_iterator< type,integral_c<long,0> > begin;
    typedef vector_iterator< type,integral_c<long,10> > end;
};

template<>
struct push_front_traits< aux::vector_tag< 9> >
{
    template< typename Vector, typename T > struct algorithm
    {
        typedef vector10<
              T
              ,
              typename Vector::item0, typename Vector::item1
            , typename Vector::item2, typename Vector::item3
            , typename Vector::item4, typename Vector::item5
            , typename Vector::item6, typename Vector::item7
            , typename Vector::item8
            > type;
    };
};

template<>
struct pop_front_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector9<
              typename Vector::item1, typename Vector::item2
            , typename Vector::item3, typename Vector::item4
            , typename Vector::item5, typename Vector::item6
            , typename Vector::item7, typename Vector::item8
            , typename Vector::item9
            > type;
    };
};

namespace aux {
template<> struct vector_item_impl<10>
{
    template< typename V_ > struct result_
    {
        typedef typename V_::item10 type;
    };
};
}

template<>
struct at_traits< aux::vector_tag< 10> >
{
    template< typename V_, typename N > struct algorithm
    {
        typedef typename aux::vector_item_impl<BOOST_MPL_AUX_VALUE_WKND(N)::value>
            ::template result_<V_>::type type;
    };
};

template<>
struct front_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::item0 type;
    };
};

template<>
struct back_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
    {
        typedef typename Vector::back type;
    };
};

template<>
struct empty_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
        : false_
    {
    };
};

template<>
struct size_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
        : integral_c< long,10 >
    {
    };
};

template<>
struct O1_size_traits< aux::vector_tag< 10> >
    : size_traits< aux::vector_tag< 10> >
{
};

template<>
struct clear_traits< aux::vector_tag< 10> >
{
    template< typename Vector > struct algorithm
    {
        typedef vector0<> type;
    };
};

} // namespace mpl
} // namespace boost

