/* supplementing compose function objects
 * Fri Jul 16 21:01:58 MEST 1999
 */
/* The following code example is taken from the book
 * "The C++ Standard Library - A Tutorial and Reference"
 * by Nicolai M. Josuttis, Addison-Wesley, 1999
 *
 * (C) Copyright Nicolai M. Josuttis 1999.
 * Permission to copy, use, modify, sell and distribute this software
 * is granted provided this copyright notice appears in all copies.
 * This software is provided "as is" without express or implied
 * warranty, and with no claim as to its suitability for any purpose.
 */

// See http://www.boost.org/libs/compose for Documentation.

#ifndef BOOST_DEPRECATED
#  error Boost.Compose has been deprecated in favor of Boost.Bind or Boost.Lambda, and will be removed in a future release. You may define the macro BOOST_DEPRECATED to suppress this warning.
#endif

#ifndef BOOST_COMPOSE_HPP
#define BOOST_COMPOSE_HPP

#include <functional>

namespace boost {

/**********************************************************
 * type nullary_function
 * - as supplement to unary_function and binary_function
 **********************************************************/
template <class Result>
struct nullary_function {
    typedef Result result_type;
};

/**********************************************************
 * ptr_fun for functions with no argument
 **********************************************************/
template <class Result>
class pointer_to_nullary_function : public nullary_function<Result>
{
  protected:
    Result (*ptr)();
  public:
    pointer_to_nullary_function() {
    }
    explicit pointer_to_nullary_function(Result (*x)()) : ptr(x) {
    }
    Result operator()() const {
        return ptr();
    }
};

template <class Result>
inline pointer_to_nullary_function<Result> ptr_fun(Result (*x)())
{
  return pointer_to_nullary_function<Result>(x);
}

/*********** compose_f_gx_t and compose_f_gx **********************/

/* class for the compose_f_gx adapter
 */
template <class OP1, class OP2>
class compose_f_gx_t
 : public std::unary_function<typename OP2::argument_type,
                              typename OP1::result_type>
{
  private:
    OP1 op1;    // process: op1(op2(x))
    OP2 op2;
  public:
    // constructor
    compose_f_gx_t(const OP1& o1, const OP2& o2)
     : op1(o1), op2(o2) {
    }

    // function call
    typename OP1::result_type
    operator()(const typename OP2::argument_type& x) const {
        return op1(op2(x));
    }
};

/* convenience functions for the compose_f_gx adapter
 */
template <class OP1, class OP2>
inline compose_f_gx_t<OP1,OP2>
compose_f_gx (const OP1& o1, const OP2& o2) {
    return compose_f_gx_t<OP1,OP2>(o1,o2);
}

/*********** compose_f_gx_hx_t and compose_f_gx_hx **********************/

/* class for the compose_f_gx_hx adapter
 */
template <class OP1, class OP2, class OP3>
class compose_f_gx_hx_t
 : public std::unary_function<typename OP2::argument_type,
                              typename OP1::result_type>
{
  private:
    OP1 op1;    // process: op1(op2(x),op3(x))
    OP2 op2;
    OP3 op3;
  public:
    // constructor
    compose_f_gx_hx_t (const OP1& o1, const OP2& o2, const OP3& o3)
     : op1(o1), op2(o2), op3(o3) {
    }

    // function call
    typename OP1::result_type
    operator()(const typename OP2::argument_type& x) const {
        return op1(op2(x),op3(x));
    }
};

/* convenience functions for the compose_f_gx_hx adapter
 */
template <class OP1, class OP2, class OP3>
inline compose_f_gx_hx_t<OP1,OP2,OP3>
compose_f_gx_hx (const OP1& o1, const OP2& o2, const OP3& o3) {
    return compose_f_gx_hx_t<OP1,OP2,OP3>(o1,o2,o3);
}

/*********** compose_f_gxy_t and compose_f_gxy **********************/

/* class for the compose_f_gxy adapter
 */
template <class OP1, class OP2>
class compose_f_gxy_t
 : public std::binary_function<typename OP2::first_argument_type,
                               typename OP2::second_argument_type,
                               typename OP1::result_type>
{
  private:
    OP1 op1;    // process: op1(op2(x,y))
    OP2 op2;
  public:
    // constructor
    compose_f_gxy_t (const OP1& o1, const OP2& o2)
     : op1(o1), op2(o2) {
    }

    // function call
    typename OP1::result_type
    operator()(const typename OP2::first_argument_type& x,
               const typename OP2::second_argument_type& y) const {
        return op1(op2(x,y));
    }
};

/* convenience function for the compose_f_gxy adapter
 */
template <class OP1, class OP2>
inline compose_f_gxy_t<OP1,OP2>
compose_f_gxy (const OP1& o1, const OP2& o2) {
    return compose_f_gxy_t<OP1,OP2>(o1,o2);
}

/*********** compose_f_gx_hy_t and compose_f_gx_hy **********************/

/* class for the compose_f_gx_hy adapter
 */
template <class OP1, class OP2, class OP3>
class compose_f_gx_hy_t
 : public std::binary_function<typename OP2::argument_type,
                               typename OP3::argument_type,
                               typename OP1::result_type>
{
  private:
    OP1 op1;    // process: op1(op2(x),op3(y))
    OP2 op2;
    OP3 op3;
  public:
    // constructor
    compose_f_gx_hy_t (const OP1& o1, const OP2& o2, const OP3& o3)
     : op1(o1), op2(o2), op3(o3) {
    }

    // function call
    typename OP1::result_type
    operator()(const typename OP2::argument_type& x,
               const typename OP3::argument_type& y) const {
        return op1(op2(x),op3(y));
    }
};

/* convenience function for the compose_f_gx_hy adapter
 */
template <class OP1, class OP2, class OP3>
inline compose_f_gx_hy_t<OP1,OP2,OP3>
compose_f_gx_hy (const OP1& o1, const OP2& o2, const OP3& o3) {
    return compose_f_gx_hy_t<OP1,OP2,OP3>(o1,o2,o3);
}

/*********** compose_f_g_t and compose_f_g **********************/

/* class for the compose_f_g adapter
 */
template <class OP1, class OP2>
class compose_f_g_t
 : public boost::nullary_function<typename OP1::result_type>
{
  private:
    OP1 op1;    // process: op1(op2())
    OP2 op2;
  public:
    // constructor
    compose_f_g_t(const OP1& o1, const OP2& o2)
     : op1(o1), op2(o2) {
    }

    // function call
    typename OP1::result_type
    operator()() const {
        return op1(op2());
    }
};

/* convenience functions for the compose_f_g adapter
 */
template <class OP1, class OP2>
inline compose_f_g_t<OP1,OP2>
compose_f_g (const OP1& o1, const OP2& o2) {
    return compose_f_g_t<OP1,OP2>(o1,o2);
}

} /* namespace boost */

#endif /*BOOST_COMPOSE_HPP*/
