// preprocessed version of 'boost/mpl/aux_/advance_forward.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {
namespace aux {

template< long N > struct advance_forward;
template<>
struct advance_forward<0>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef iter0 type;
    };
};

template<>
struct advance_forward<1>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::next iter1;
        typedef iter1 type;
    };
};

template<>
struct advance_forward<2>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::next iter1;
        typedef typename iter1::next iter2;
        typedef iter2 type;
    };
};

template<>
struct advance_forward<3>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::next iter1;
        typedef typename iter1::next iter2;
        typedef typename iter2::next iter3;
        typedef iter3 type;
    };
};

template<>
struct advance_forward<4>
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;
        typedef typename iter0::next iter1;
        typedef typename iter1::next iter2;
        typedef typename iter2::next iter3;
        typedef typename iter3::next iter4;
        typedef iter4 type;
    };
};

template< long N >
struct advance_forward
{
    template< typename Iterator > struct apply
    {
        typedef typename advance_forward<4>::template apply<Iterator>::type chunk_result_;
        typedef typename advance_forward<( (N - 4) < 0 ? 0 : N - 4 )>::template apply<chunk_result_>::type type;
    };
};

} // namespace aux
} // namespace mpl
} // namespace boost

