// preprocessed version of 'boost/mpl/lambda_helper.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      template< typename P1 > class F
    , typename T1
    >
struct lambda_helper1
{
    struct rebind
    {
        enum { arity = 1 };
        typedef T1 arg1;

        template< typename U1 > struct apply
            : F<U1>
        {
        };
   };
};

template<
      template< typename P1, typename P2 > class F
    , typename T1, typename T2
    >
struct lambda_helper2
{
    struct rebind
    {
        enum { arity = 2 };
        typedef T1 arg1;
        typedef T2 arg2;
        
        template< typename U1, typename U2 > struct apply
            : F< U1,U2 >
        {
        };
   };
};

template<
      template< typename P1, typename P2, typename P3 > class F
    , typename T1, typename T2, typename T3
    >
struct lambda_helper3
{
    struct rebind
    {
        enum { arity = 3 };
        typedef T1 arg1;
        typedef T2 arg2;
        typedef T3 arg3;
        
        template< typename U1, typename U2, typename U3 > struct apply
            : F< U1,U2,U3 >
        {
        };
   };
};

template<
      template< typename P1, typename P2, typename P3, typename P4 > class F
    , typename T1, typename T2, typename T3, typename T4
    >
struct lambda_helper4
{
    struct rebind
    {
        enum { arity = 4 };
        typedef T1 arg1;
        typedef T2 arg2;
        typedef T3 arg3;
        typedef T4 arg4;
        
        template<
              typename U1, typename U2, typename U3, typename U4
            >
        struct apply
            : F< U1,U2,U3,U4 >
        {
        };
   };
};

template<
      template<
          typename P1, typename P2, typename P3, typename P4
        , typename P5
        >
      class F
    , typename T1, typename T2, typename T3, typename T4, typename T5
    >
struct lambda_helper5
{
    struct rebind
    {
        enum { arity = 5 };
        typedef T1 arg1;
        typedef T2 arg2;
        typedef T3 arg3;
        typedef T4 arg4;
        typedef T5 arg5;
        
        template<
              typename U1, typename U2, typename U3, typename U4
            , typename U5
            >
        struct apply
            : F< U1,U2,U3,U4,U5 >
        {
        };
   };
};

} // namespace mpl
} // namespace boost

