// preprocessed version of 'boost/mpl/vector_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T
    , long C0 = LONG_MAX, long C1 = LONG_MAX, long C2 = LONG_MAX
    , long C3 = LONG_MAX, long C4 = LONG_MAX, long C5 = LONG_MAX
    , long C6 = LONG_MAX, long C7 = LONG_MAX, long C8 = LONG_MAX
    , long C9 = LONG_MAX
    >
struct vector_c;

template<
      typename T
     
    >
struct vector_c<
          T, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        >
    : vector0_c<T>
{
    typedef vector0_c<T> type;
};

template<
      typename T
    , long C0
    >
struct vector_c<
          T, C0, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX, LONG_MAX, LONG_MAX
        >
    : vector1_c< T,C0 >
{
    typedef vector1_c< T,C0 > type;
};

template<
      typename T
    , long C0, long C1
    >
struct vector_c<
          T, C0, C1, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX, LONG_MAX, LONG_MAX
        >
    : vector2_c< T,C0,C1 >
{
    typedef vector2_c< T,C0,C1 > type;
};

template<
      typename T
    , long C0, long C1, long C2
    >
struct vector_c<
          T, C0, C1, C2, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX, LONG_MAX
        >
    : vector3_c< T,C0,C1,C2 >
{
    typedef vector3_c< T,C0,C1,C2 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3
    >
struct vector_c<
          T, C0, C1, C2, C3, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX
        >
    : vector4_c< T,C0,C1,C2,C3 >
{
    typedef vector4_c< T,C0,C1,C2,C3 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4
    >
struct vector_c<
          T, C0, C1, C2, C3, C4, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        , LONG_MAX
        >
    : vector5_c< T,C0,C1,C2,C3,C4 >
{
    typedef vector5_c< T,C0,C1,C2,C3,C4 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5
    >
struct vector_c<
          T, C0, C1, C2, C3, C4, C5, LONG_MAX, LONG_MAX, LONG_MAX, LONG_MAX
        >
    : vector6_c< T,C0,C1,C2,C3,C4,C5 >
{
    typedef vector6_c< T,C0,C1,C2,C3,C4,C5 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5, long C6
    >
struct vector_c<
          T, C0, C1, C2, C3, C4, C5, C6, LONG_MAX, LONG_MAX, LONG_MAX
        >
    : vector7_c< T,C0,C1,C2,C3,C4,C5,C6 >
{
    typedef vector7_c< T,C0,C1,C2,C3,C4,C5,C6 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5, long C6, long C7
    >
struct vector_c< T,C0,C1,C2,C3,C4,C5,C6,C7,LONG_MAX,LONG_MAX >
    : vector8_c< T,C0,C1,C2,C3,C4,C5,C6,C7 >
{
    typedef vector8_c< T,C0,C1,C2,C3,C4,C5,C6,C7 > type;
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5, long C6, long C7
    , long C8
    >
struct vector_c< T,C0,C1,C2,C3,C4,C5,C6,C7,C8,LONG_MAX >
    : vector9_c< T,C0,C1,C2,C3,C4,C5,C6,C7,C8 >
{
    typedef vector9_c< T,C0,C1,C2,C3,C4,C5,C6,C7,C8 > type;
};

// primary template (not a specialization!)
template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5, long C6, long C7
    , long C8, long C9
    >
struct vector_c
    : vector10_c< T,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 >
{
    typedef vector10_c< T,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 > type;
};

} // namespace mpl
} // namespace boost

