// preprocessed version of 'boost/mpl/aux_/template_arity.hpp' header
// see the original for copyright information

namespace boost { namespace mpl { namespace aux {

template< bool >
struct template_arity_impl
{
    template< typename F > struct result_
    {
        enum { value = -1 };
    };
};

template<>
struct template_arity_impl<true>
{
    template< typename F > struct result_
    {
        enum { value = F::arity };

    };
};

template< typename F >
struct template_arity
    : template_arity_impl< ::boost::mpl::aux::has_rebind<F>::value >
        ::template result_<F>
{
};

template<>
struct template_arity<int>
{
    enum { value = -1 };
};

}}} // namespace boost::mpl::aux

