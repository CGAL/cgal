// Boost Lambda Library -- control_constructs_common.hpp -------------------

// Copyright (C) 1999, 2000 Jaakko Järvi (jaakko.jarvi@cs.utu.fi)
// Copyright (C) 2000 Gary Powell (powellg@amazon.com)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies. 
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice 
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty, 
// and with no claim as to its suitability for any purpose.
//
// For more information, see www.boost.org

// --------------------------------------------------------------------------

#if !defined(BOOST_CONTROL_CONSTRUCTS_COMMON_HPP)
#define BOOST_CONTROL_CONSTRUCTS_COMMON_HPP

namespace boost { 
namespace lambda {

  // special types of lambda functors, used with control structures
  // to guarantee that they are composed correctly.

template<class Tag, class LambdaFunctor>
class tagged_lambda_functor;

template<class Tag, class Args>
class tagged_lambda_functor<Tag, lambda_functor<Args> > 
  : public lambda_functor<Args> 
{
public:
  tagged_lambda_functor(const Args& a) : lambda_functor<Args>(a) {}

  tagged_lambda_functor(const lambda_functor<Args>& a) 
    : lambda_functor<Args>(a) {}

  // for the no body cases in control structures.
  tagged_lambda_functor() : lambda_functor<Args>() {}
};

} // lambda
} // boost

#endif // BOOST_CONTROL_CONSTRUCTS_COMMON_HPP







