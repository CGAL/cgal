// preprocessed version of 'boost/mpl/vector.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T0 = void_, typename T1 = void_, typename T2 = void_
    , typename T3 = void_, typename T4 = void_, typename T5 = void_
    , typename T6 = void_, typename T7 = void_, typename T8 = void_
    , typename T9 = void_
    >
struct vector;

template<
      
    >
struct vector<
          void_, void_, void_, void_, void_, void_, void_, void_, void_
        , void_
        >
    : vector0<  >
{
    typedef vector0<  > type;
};

template<
      typename T0
    >
struct vector<
          T0, void_, void_, void_, void_, void_, void_, void_, void_, void_
        >
    : vector1<T0>
{
    typedef vector1<T0> type;
};

template<
      typename T0, typename T1
    >
struct vector<
          T0, T1, void_, void_, void_, void_, void_, void_, void_, void_
        >
    : vector2< T0,T1 >
{
    typedef vector2< T0,T1 > type;
};

template<
      typename T0, typename T1, typename T2
    >
struct vector< T0,T1,T2,void_,void_,void_,void_,void_,void_,void_ >
    : vector3< T0,T1,T2 >
{
    typedef vector3< T0,T1,T2 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3
    >
struct vector< T0,T1,T2,T3,void_,void_,void_,void_,void_,void_ >
    : vector4< T0,T1,T2,T3 >
{
    typedef vector4< T0,T1,T2,T3 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    >
struct vector< T0,T1,T2,T3,T4,void_,void_,void_,void_,void_ >
    : vector5< T0,T1,T2,T3,T4 >
{
    typedef vector5< T0,T1,T2,T3,T4 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct vector< T0,T1,T2,T3,T4,T5,void_,void_,void_,void_ >
    : vector6< T0,T1,T2,T3,T4,T5 >
{
    typedef vector6< T0,T1,T2,T3,T4,T5 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6
    >
struct vector< T0,T1,T2,T3,T4,T5,T6,void_,void_,void_ >
    : vector7< T0,T1,T2,T3,T4,T5,T6 >
{
    typedef vector7< T0,T1,T2,T3,T4,T5,T6 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7
    >
struct vector< T0,T1,T2,T3,T4,T5,T6,T7,void_,void_ >
    : vector8< T0,T1,T2,T3,T4,T5,T6,T7 >
{
    typedef vector8< T0,T1,T2,T3,T4,T5,T6,T7 > type;
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8
    >
struct vector< T0,T1,T2,T3,T4,T5,T6,T7,T8,void_ >
    : vector9< T0,T1,T2,T3,T4,T5,T6,T7,T8 >
{
    typedef vector9< T0,T1,T2,T3,T4,T5,T6,T7,T8 > type;
};

// primary template (not a specialization!)
template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    >
struct vector
    : vector10< T0,T1,T2,T3,T4,T5,T6,T7,T8,T9 >
{
    typedef vector10< T0,T1,T2,T3,T4,T5,T6,T7,T8,T9 > type;
};

} // namespace mpl
} // namespace boost

