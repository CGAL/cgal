// preprocessed version of 'boost/mpl/vector/vector10_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T
    , T C0
    >
struct vector1_c
    : vector1< integral_c<T,C0> >
{
};

template<
      typename T
    , T C0, T C1
    >
struct vector2_c
    : vector2< integral_c<T,C0>,integral_c<T,C1> >
{
};

template<
      typename T
    , T C0, T C1, T C2
    >
struct vector3_c
    : vector3< integral_c<T,C0>,integral_c<T,C1>,integral_c<T,C2> >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3
    >
struct vector4_c
    : vector4<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4
    >
struct vector5_c
    : vector5<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5
    >
struct vector6_c
    : vector6<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4>,integral_c<T,C5 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6
    >
struct vector7_c
    : vector7<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4>,integral_c<T,C5 >
        ,integral_c< T,C6 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7
    >
struct vector8_c
    : vector8<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4>,integral_c<T,C5 >
        ,integral_c< T,C6>,integral_c<T,C7 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8
    >
struct vector9_c
    : vector9<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4>,integral_c<T,C5 >
        ,integral_c< T,C6>,integral_c<T,C7>,integral_c<T,C8 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9
    >
struct vector10_c
    : vector10<
          integral_c< T,C0>,integral_c<T,C1>,integral_c<T,C2 >
        ,integral_c< T,C3>,integral_c<T,C4>,integral_c<T,C5 >
        ,integral_c< T,C6>,integral_c<T,C7>,integral_c<T,C8 >
        ,integral_c< T,C9 >
        >
{
};

} // namespace mpl
} // namespace boost

