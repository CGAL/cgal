// preprocessed version of 'boost/mpl/list.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T0 = void_, typename T1 = void_, typename T2 = void_
    , typename T3 = void_, typename T4 = void_, typename T5 = void_
    , typename T6 = void_, typename T7 = void_, typename T8 = void_
    , typename T9 = void_
    >
struct list;

template<
     
    >
struct list<
          void_, void_, void_, void_, void_, void_, void_, void_, void_
        , void_
        >
    : list0<  >
{
    typedef list0<  >::type type;
};

template<
      typename T0
    >
struct list<
          T0, void_, void_, void_, void_, void_, void_, void_, void_, void_
        >
    : list1<T0>
{
    typedef typename list1<T0>::type type;
};

template<
      typename T0, typename T1
    >
struct list< T0,T1,void_,void_,void_,void_,void_,void_,void_,void_ >
    : list2< T0,T1 >
{
    typedef typename list2< T0,T1 >::type type;
};

template<
      typename T0, typename T1, typename T2
    >
struct list< T0,T1,T2,void_,void_,void_,void_,void_,void_,void_ >
    : list3< T0,T1,T2 >
{
    typedef typename list3< T0,T1,T2 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3
    >
struct list< T0,T1,T2,T3,void_,void_,void_,void_,void_,void_ >
    : list4< T0,T1,T2,T3 >
{
    typedef typename list4< T0,T1,T2,T3 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    >
struct list< T0,T1,T2,T3,T4,void_,void_,void_,void_,void_ >
    : list5< T0,T1,T2,T3,T4 >
{
    typedef typename list5< T0,T1,T2,T3,T4 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct list< T0,T1,T2,T3,T4,T5,void_,void_,void_,void_ >
    : list6< T0,T1,T2,T3,T4,T5 >
{
    typedef typename list6< T0,T1,T2,T3,T4,T5 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6
    >
struct list< T0,T1,T2,T3,T4,T5,T6,void_,void_,void_ >
    : list7< T0,T1,T2,T3,T4,T5,T6 >
{
    typedef typename list7< T0,T1,T2,T3,T4,T5,T6 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7
    >
struct list< T0,T1,T2,T3,T4,T5,T6,T7,void_,void_ >
    : list8< T0,T1,T2,T3,T4,T5,T6,T7 >
{
    typedef typename list8< T0,T1,T2,T3,T4,T5,T6,T7 >::type type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8
    >
struct list< T0,T1,T2,T3,T4,T5,T6,T7,T8,void_ >
    : list9< T0,T1,T2,T3,T4,T5,T6,T7,T8 >
{
    typedef typename list9< T0,T1,T2,T3,T4,T5,T6,T7,T8 >::type type;
};

// primary template (not a specialization!)

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    >
struct list
    : list10< T0,T1,T2,T3,T4,T5,T6,T7,T8,T9 >
{
    typedef list10< T0,T1,T2,T3,T4,T5,T6,T7,T8,T9 > type;
};

} // namespace mpl
} // namespace boost

