// preprocessed version of 'boost/mpl/vector/vector10.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T0
    >
struct vector1
    : vector_node<
          1
        , T0
        , vector0<  >
        >
{
};

template<
      typename T0, typename T1
    >
struct vector2
    : vector_node<
          2
        , T0
        , vector1<T1>
        >
{
};

template<
      typename T0, typename T1, typename T2
    >
struct vector3
    : vector_node<
          3
        , T0
        , vector2< T1,T2 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3
    >
struct vector4
    : vector_node<
          4
        , T0
        , vector3< T1,T2,T3 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    >
struct vector5
    : vector_node<
          5
        , T0
        , vector4< T1,T2,T3,T4 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct vector6
    : vector_node<
          6
        , T0
        , vector5< T1,T2,T3,T4,T5 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6
    >
struct vector7
    : vector_node<
          7
        , T0
        , vector6< T1,T2,T3,T4,T5,T6 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7
    >
struct vector8
    : vector_node<
          8
        , T0
        , vector7< T1,T2,T3,T4,T5,T6,T7 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8
    >
struct vector9
    : vector_node<
          9
        , T0
        , vector8< T1,T2,T3,T4,T5,T6,T7,T8 >
        >
{
};

template<
      typename T0, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename T6, typename T7, typename T8, typename T9
    >
struct vector10
    : vector_node<
          10
        , T0
        , vector9< T1,T2,T3,T4,T5,T6,T7,T8,T9 >
        >
{
};

} // namespace mpl
} // namespace boost

