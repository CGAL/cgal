//  © Copyright Fernando Luis Cacciola Carballal 2000-2004
//  Use, modification, and distribution is subject to the Boost Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See library home page at http://www.boost.org/libs/numeric/conversion
//
// Contact the author at: fernando_cacciola@hotmail.com
// 
//
//  Revision History
//
//    19 Nov 2001 Syntatic changes as suggested by Darin Adler (Fernando Cacciola)
//    08 Nov 2001 Fixes to accommodate MSVC (Fernando Cacciola)
//    04 Nov 2001 Fixes to accommodate gcc2.92 (Fernando Cacciola)
//    30 Oct 2001 Some fixes suggested by Daryle Walker (Fernando Cacciola)
//    25 Oct 2001 Initial boostification (Fernando Cacciola)
//    23 Jan 2004 Inital add to cvs (post review)s
//
#ifndef BOOST_NUMERIC_CONVERSION_CAST_25OCT2001_HPP
#define BOOST_NUMERIC_CONVERSION_CAST_25OCT2001_HPP

#include <boost/type.hpp>
#include <boost/numeric/conversion/converter.hpp>

#ifndef BOOST_EXPLICIT_DEFAULT_TARGET
# if defined(BOOST_MSVC) && BOOST_MSVC <= 1200 // 1200 = VC6
#  define BOOST_EXPLICIT_DEFAULT_TARGET , ::boost::type<Target>* = 0
# else
#  define BOOST_EXPLICIT_DEFAULT_TARGET
# endif
#endif

namespace boost
{
  template<typename Target, typename Source>
  inline
  typename boost::numeric::converter<Target,Source>::result_type
    numeric_cast ( Source arg BOOST_EXPLICIT_DEFAULT_TARGET )
  {
    typedef boost::numeric::converter<Target,Source> Converter ;
    return Converter::convert(arg);
  }

} // namespace boost

#endif
