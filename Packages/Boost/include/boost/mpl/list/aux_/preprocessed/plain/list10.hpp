// preprocessed version of 'boost/mpl/list/list10.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T0
    >
struct list1
    : list_node<
          integral_c< long,1 >
        , T0
        , null_node
        >
{
    typedef list1 type;
};

template<
      typename T0, typename T1
    >
struct list2
    : list_node<
          integral_c< long,2 >
        , T0
        , list1<T1>
        >
{
    typedef list2 type;
};

template<
      typename T0, typename T1, typename T2
    >
struct list3
    : list_node<
          integral_c< long,3 >
        , T0
        , list2< T1,T2 >
        >
{
    typedef list3 type;
};

template<
      typename T0, typename T1, typename T2, typename T3
    >
struct list4
    : list_node<
          integral_c< long,4 >
        , T0
        , list3< T1,T2,T3 >
        >
{
    typedef list4 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    >
struct list5
    : list_node<
          integral_c< long,5 >
        , T0
        , list4< T1,T2,T3,T4 >
        >
{
    typedef list5 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct list6
    : list_node<
          integral_c< long,6 >
        , T0
        , list5< T1,T2,T3,T4,T5 >
        >
{
    typedef list6 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6
    >
struct list7
    : list_node<
          integral_c< long,7 >
        , T0
        , list6< T1,T2,T3,T4,T5,T6 >
        >
{
    typedef list7 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7
    >
struct list8
    : list_node<
          integral_c< long,8 >
        , T0
        , list7< T1,T2,T3,T4,T5,T6,T7 >
        >
{
    typedef list8 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8
    >
struct list9
    : list_node<
          integral_c< long,9 >
        , T0
        , list8< T1,T2,T3,T4,T5,T6,T7,T8 >
        >
{
    typedef list9 type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    >
struct list10
    : list_node<
          integral_c< long,10 >
        , T0
        , list9< T1,T2,T3,T4,T5,T6,T7,T8,T9 >
        >
{
    typedef list10 type;
};

} // namespace mpl
} // namespace boost

