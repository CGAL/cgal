// preprocessed version of 'boost/mpl/list_c.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

namespace aux {
template< int N > struct list_c_impl_chooser;
}

namespace aux {

template<>
struct list_c_impl_chooser<0>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list0_c<
              T  
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<1>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list1_c<
              T, C0
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<2>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list2_c<
              T, C0, C1
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<3>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list3_c<
              T, C0, C1, C2
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<4>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list4_c<
              T, C0, C1, C2, C3
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<5>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list5_c<
              T, C0, C1, C2, C3, C4
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<6>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list6_c<
              T, C0, C1, C2, C3, C4, C5
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<7>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list7_c<
              T, C0, C1, C2, C3, C4, C5, C6
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<8>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list8_c<
              T, C0, C1, C2, C3, C4, C5, C6, C7
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<9>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list9_c<
              T, C0, C1, C2, C3, C4, C5, C6, C7, C8
            >::type type;
    };
};

} // namespace aux

namespace aux {

template<>
struct list_c_impl_chooser<10>
{
    template<
          typename T
        , long C0, long C1, long C2, long C3, long C4, long C5, long C6
        , long C7, long C8, long C9
        >
    struct result_
    {
        typedef typename list10_c<
              T, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9
            >::type type;
    };
};

} // namespace aux

namespace aux {

template< long T >
struct is_list_c_arg
{
    enum { value = true };
};

template<>
struct is_list_c_arg<LONG_MAX>
{
    enum { value = false };
};

template<
      long T1, long T2, long T3, long T4, long T5, long T6, long T7, long T8
    , long T9, long T10
    >
struct list_c_count_args
{
    enum { value =
          is_list_c_arg<T1>::value + is_list_c_arg<T2>::value 
        + is_list_c_arg<T3>::value + is_list_c_arg<T4>::value 
        + is_list_c_arg<T5>::value + is_list_c_arg<T6>::value 
        + is_list_c_arg<T7>::value + is_list_c_arg<T8>::value 
        + is_list_c_arg<T9>::value + is_list_c_arg<T10>::value
        };
};

template<
      typename T
    , long C0, long C1, long C2, long C3, long C4, long C5, long C6, long C7
    , long C8, long C9
    >
struct list_c_impl
{
    typedef aux::list_c_count_args< C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 > arg_num_;
    typedef typename aux::list_c_impl_chooser< arg_num_::value >
        ::template result_< T,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 >::type type;
};

} // namespace aux

template<
      typename T
    , long C0 = LONG_MAX, long C1 = LONG_MAX, long C2 = LONG_MAX
    , long C3 = LONG_MAX, long C4 = LONG_MAX, long C5 = LONG_MAX
    , long C6 = LONG_MAX, long C7 = LONG_MAX, long C8 = LONG_MAX
    , long C9 = LONG_MAX
    >
struct list_c
    : aux::list_c_impl< T,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9 >::type
{
    typedef typename aux::list_c_impl<
          T, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9
        >::type type;
};

} // namespace mpl
} // namespace boost

