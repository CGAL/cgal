// preprocessed version of 'boost/mpl/vector/vector20_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    >
struct vector11_c
    : vector_node<
          11
        , integral_c< T,C0 >
        , vector10_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11
    >
struct vector12_c
    : vector_node<
          12
        , integral_c< T,C0 >
        , vector11_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12
    >
struct vector13_c
    : vector_node<
          13
        , integral_c< T,C0 >
        , vector12_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13
    >
struct vector14_c
    : vector_node<
          14
        , integral_c< T,C0 >
        , vector13_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14
    >
struct vector15_c
    : vector_node<
          15
        , integral_c< T,C0 >
        , vector14_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14, T C15
    >
struct vector16_c
    : vector_node<
          16
        , integral_c< T,C0 >
        , vector15_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14, T C15, T C16
    >
struct vector17_c
    : vector_node<
          17
        , integral_c< T,C0 >
        , vector16_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14, T C15, T C16, T C17
    >
struct vector18_c
    : vector_node<
          18
        , integral_c< T,C0 >
        , vector17_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14, T C15, T C16, T C17, T C18
    >
struct vector19_c
    : vector_node<
          19
        , integral_c< T,C0 >
        , vector18_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18 >
        >
{
};

template<
      typename T
    , T C0, T C1, T C2, T C3, T C4, T C5, T C6, T C7, T C8, T C9, T C10
    , T C11, T C12, T C13, T C14, T C15, T C16, T C17, T C18, T C19
    >
struct vector20_c
    : vector_node<
          20
        , integral_c< T,C0 >
        , vector19_c< T,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19 >
        >
{
};

} // namespace mpl
} // namespace boost

