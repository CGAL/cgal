// preprocessed version of 'boost/mpl/aux_/config/use_preprocessed.hpp' header
// see the original for copyright information

namespace boost { namespace mpl {

namespace aux {

template< bool C_ > struct or_impl
{
    template<
          typename T1, typename T2, typename T3, typename T4
        >
    struct result_
        : true_
    {
    };
};

template<> struct or_impl<false>
{
    template<
          typename T1, typename T2, typename T3, typename T4
        >
    struct result_
        : or_impl< 
              BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
            >::template result_< T2,T3,T4,false_ >
    {
    };

};

template<>
struct or_impl<false>
    ::result_< false_,false_,false_,false_ >
        : false_
{
};

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    , typename T3 = false_, typename T4 = false_, typename T5 = false_
    >
struct or_

    : aux::or_impl< 
          BOOST_MPL_AUX_NESTED_TYPE_WKND(T1)::value
        >::template result_< T2,T3,T4,T5 >

{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          5
        , or_
        , (T1, T2, T3, T4, T5)
        )
};

BOOST_MPL_AUX_VOID_SPEC_EXT(
      2
    , 5
    , or_
    )

}} // namespace boost::mpl

