// preprocessed version of 'boost/mpl/arg.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<> struct arg<-1>
{
    enum { value = -1 };
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U1 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

template<> struct arg<1>
{
    enum { value = 1 };
    typedef arg<2> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U1 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

template<> struct arg<2>
{
    enum { value = 2 };
    typedef arg<3> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U2 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

template<> struct arg<3>
{
    enum { value = 3 };
    typedef arg<4> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U3 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

template<> struct arg<4>
{
    enum { value = 4 };
    typedef arg<5> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U4 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

template<> struct arg<5>
{
    enum { value = 5 };
    typedef arg<6> next;
    BOOST_MPL_AUX_ARG_TYPEDEF(void_, tag)

    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
        typedef U5 type;

     private:
        enum { nv = !is_void_<type>::value };
        BOOST_STATIC_ASSERT(nv);

    };
};

BOOST_MPL_AUX_NONTYPE_ARITY_SPEC(1, int, arg)

} // namespace mpl
} // namespace boost

