// preprocessed version of 'boost/mpl/vector/vector10_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T
    , T C0
    >
struct vector1_c
    : vector_node<
          1
        , integral_c< T,C0 >
        , vector0_c<T>
        >
{
};

template<
      typename T
    , T C0, T C1
    >
struct vector2_c
    : vector_node<
          2
        , integral_c< T,C0 >
        , vector1_c< T,C1 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2
    >
struct vector3_c
    : vector_node<
          3
        , integral_c< T,C0 >
        , vector2_c< T,C1,C2 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3
    >
struct vector4_c
    : vector_node<
          4
        , integral_c< T,C0 >
        , vector3_c< T,C1,C2,C3 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4
    >
struct vector5_c
    : vector_node<
          5
        , integral_c< T,C0 >
        , vector4_c< T,C1,C2,C3,C4 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5
    >
struct vector6_c
    : vector_node<
          6
        , integral_c< T,C0 >
        , vector5_c< T,C1,C2,C3,C4,C5 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6
    >
struct vector7_c
    : vector_node<
          7
        , integral_c< T,C0 >
        , vector6_c< T,C1,C2,C3,C4,C5,C6 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7
    >
struct vector8_c
    : vector_node<
          8
        , integral_c< T,C0 >
        , vector7_c< T,C1,C2,C3,C4,C5,C6,C7 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8
    >
struct vector9_c
    : vector_node<
          9
        , integral_c< T,C0 >
        , vector8_c< T,C1,C2,C3,C4,C5,C6,C7,C8 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9
    >
struct vector10_c
    : vector_node<
          10
        , integral_c< T,C0 >
        , vector9_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9 >
        >
{
};

} // namespace mpl
} // namespace boost

