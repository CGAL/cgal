// preprocessed version of 'boost/mpl/aux_/config/use_preprocessed.hpp' header
// see the original for copyright information

namespace boost { namespace mpl {

namespace aux {

template< bool C_ > struct and_impl
{
    template<
          typename T1, typename T2, typename T3, typename T4
        >
    struct result_
        : false_
    {
    };
};

template<> struct and_impl<true>
{
    template<
          typename T1, typename T2, typename T3, typename T4
        >
    struct result_
        : and_impl< 
              BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
            >::template result_< T2,T3,T4,true_ >
    {
    };

    template<> struct result_<true_, true_, true_, true_>
        : true_
    {
    };
};

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    , typename T3 = true_, typename T4 = true_, typename T5 = true_
    >
struct and_

    : aux::and_impl< 
          BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
        >::template result_< T2,T3,T4,T5 >

{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          5
        , and_
        , (T1, T2, T3, T4, T5)
        )
};

BOOST_MPL_AUX_VOID_SPEC_EXT(
      2
    , 5
    , and_
    )

}} // namespace boost::mpl

