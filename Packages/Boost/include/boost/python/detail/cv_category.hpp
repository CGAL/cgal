// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CV_CATEGORY_DWA200222_HPP
# define CV_CATEGORY_DWA200222_HPP
# include <boost/type_traits/cv_traits.hpp>

namespace boost { namespace python { namespace detail { 

template <bool is_const_, bool is_volatile_>
struct cv_tag
{
    BOOST_STATIC_CONSTANT(bool, is_const = is_const_);
    BOOST_STATIC_CONSTANT(bool, is_volatile = is_const_);
};

typedef cv_tag<false,false> cv_unqualified;
typedef cv_tag<true,false> const_;
typedef cv_tag<false,true> volatile_;
typedef cv_tag<true,true> const_volatile_;

template <class T>
struct cv_category
{
//    BOOST_STATIC_CONSTANT(bool, c = is_const<T>::value);
//    BOOST_STATIC_CONSTANT(bool, v = is_volatile<T>::value);
    typedef cv_tag<
        ::boost::is_const<T>::value
      , ::boost::is_volatile<T>::value
    > type;
};

}}} // namespace boost::python::detail

#endif // CV_CATEGORY_DWA200222_HPP
