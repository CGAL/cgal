// Copyright Daniel Wallin, David Abrahams 2005. Use, modification and
// distribution is subject to the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef TAGGED_ARGUMENT_050328_HPP
#define TAGGED_ARGUMENT_050328_HPP

#include <boost/parameter/aux_/void.hpp>
#include <boost/parameter/aux_/arg_list.hpp>
#include <boost/detail/is_xxx.hpp>

namespace boost { namespace parameter { namespace aux {

// Holds a reference to an argument of type Arg associated with
// keyword Keyword
    
template <class Keyword, class Arg>
struct tagged_argument
{
    typedef Keyword key_type;
    typedef Arg value_type;
    typedef Arg& reference;

    tagged_argument(reference x) : value(x) {}

    // Comma operator to compose argument list without using parameters<>.
    // Useful for argument lists with undetermined length.
    template <class Keyword2, class Arg2>
    arg_list<
        tagged_argument<Keyword, Arg>
      , arg_list<tagged_argument<Keyword2, Arg2> > 
    >
    operator,(tagged_argument<Keyword2, Arg2> x) const
    {
        return arg_list<
            tagged_argument<Keyword, Arg>
          , arg_list<tagged_argument<Keyword2, Arg2> > 
        >(
            *this
          , arg_list<tagged_argument<Keyword2, Arg2> >(x, empty_arg_list())
        );
    }
    
    reference value;
};

// Defines a metafunction, is_tagged_argument, that identifies
// tagged_argument specializations.
BOOST_DETAIL_IS_XXX_DEF(tagged_argument,tagged_argument,2)

}}} // namespace boost::parameter::aux

#endif // TAGGED_ARGUMENT_050328_HPP

