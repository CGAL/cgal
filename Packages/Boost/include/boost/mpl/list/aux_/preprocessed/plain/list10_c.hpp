// preprocessed version of 'boost/mpl/list/list10_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T
    , T C0
    >
struct list1_c
    : list_node<
          integral_c< long,1 >
        , integral_c< T,C0 >
        , null_node
        >
{
    typedef list1_c type;
};

template<
      typename T
    , T C0, T C1
    >
struct list2_c
    : list_node<
          integral_c< long,2 >
        , integral_c< T,C0 >
        , list1_c< T,C1 >
        >
{
    typedef list2_c type;
};

template<
      typename T
    , T C0, T C1, T C2
    >
struct list3_c
    : list_node<
          integral_c< long,3 >
        , integral_c< T,C0 >
        , list2_c< T,C1,C2 >
        >
{
    typedef list3_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3
    >
struct list4_c
    : list_node<
          integral_c< long,4 >
        , integral_c< T,C0 >
        , list3_c< T,C1,C2,C3 >
        >
{
    typedef list4_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4
    >
struct list5_c
    : list_node<
          integral_c< long,5 >
        , integral_c< T,C0 >
        , list4_c< T,C1,C2,C3,C4 >
        >
{
    typedef list5_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5
    >
struct list6_c
    : list_node<
          integral_c< long,6 >
        , integral_c< T,C0 >
        , list5_c< T,C1,C2,C3,C4,C5 >
        >
{
    typedef list6_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6
    >
struct list7_c
    : list_node<
          integral_c< long,7 >
        , integral_c< T,C0 >
        , list6_c< T,C1,C2,C3,C4,C5,C6 >
        >
{
    typedef list7_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7
    >
struct list8_c
    : list_node<
          integral_c< long,8 >
        , integral_c< T,C0 >
        , list7_c< T,C1,C2,C3,C4,C5,C6,C7 >
        >
{
    typedef list8_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8
    >
struct list9_c
    : list_node<
          integral_c< long,9 >
        , integral_c< T,C0 >
        , list8_c< T,C1,C2,C3,C4,C5,C6,C7,C8 >
        >
{
    typedef list9_c type;
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9
    >
struct list10_c
    : list_node<
          integral_c< long,10 >
        , integral_c< T,C0 >
        , list9_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9 >
        >
{
    typedef list10_c type;
};

} // namespace mpl
} // namespace boost

