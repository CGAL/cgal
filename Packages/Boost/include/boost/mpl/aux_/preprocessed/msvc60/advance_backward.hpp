// preprocessed version of 'boost/mpl/aux_/advance_backward.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {
namespace aux {

template< nttp_long N > struct advance_backward;
template<>
struct advance_backward<0>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef iter0 type;
    };

    // ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };

};

template<>
struct advance_backward<1>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::prior iter1;
        typedef iter1 type;
    };

    // ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };

};

template<>
struct advance_backward<2>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::prior iter1;
        typedef typename iter1::prior iter2;
        typedef iter2 type;
    };

    // ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };

};

template<>
struct advance_backward<3>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::prior iter1;
        typedef typename iter1::prior iter2;
        typedef typename iter2::prior iter3;
        typedef iter3 type;
    };

    // ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };

};

template<>
struct advance_backward<4>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::prior iter1;
        typedef typename iter1::prior iter2;
        typedef typename iter2::prior iter3;
        typedef typename iter3::prior iter4;
        typedef iter4 type;
    };

    // ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };

};

template< nttp_long N >
struct advance_backward
{
    template< typename Iterator > struct apply
    {
        typedef typename apply1< advance_backward<4>,Iterator >::type chunk_result_;
        typedef typename apply1<advance_backward<( (N - 4) < 0 ? 0 : N - 4 )>,chunk_result_>::type type;
    };
};

} // namespace aux
} // namespace mpl
} // namespace boost

