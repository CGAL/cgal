// -*- C++ -*-
//  Boost general library 'format'   ---------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2001
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ------------------------------------------------------------------------------
// config_macros.hpp : configuration macros for the format library
//   only BOOST_IO_STD is absolutely needed. other are just used to trigger workaround
//   codes here and there.
// ------------------------------------------------------------------------------

#ifndef BOOST_FORMAT_CONFIG_MACROS_HPP
#define BOOST_FORMAT_CONFIG_MACROS_HPP

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

// make sure our local macros wont override something :
#if defined(BOOST_NO_LOCALE_ISDIGIT) || defined(BOOST_OVERLOAD_FOR_NON_CONST) \
  || defined(BOOST_IO_STD) || defined( BOOST_IO_NEEDS_USING_DECLARATION )
#error "boost::format defines a local macro that would overwrite a previously defined macro."
#endif

// specific workarounds. each header can define BOOS_IO_STD if it 
// needs. (e.g. because of IO_NEEDS_USING_DECLARATION)
#include <boost/format/detail/workarounds_gcc-2.95.hpp>
#include <boost/format/detail/workarounds_stlport.hpp>  // stlport workarounds

#ifndef BOOST_IO_STD
#  define BOOST_IO_STD std::
#endif

#if defined(BOOST_NO_STD_LOCALE) || \
 ( BOOST_WORKAROUND(__BORLANDC__, <= 0x564) \
   || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT( 0x570 ) )  )
// some future __BORLANDC__ >0x564  versions might not need this
// 0x570 is Borland's kylix branch
#define BOOST_NO_LOCALE_ISIDIGIT
#endif

#if  BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x570) ) || BOOST_WORKAROUND( BOOST_MSVC, BOOST_TESTED_AT(1300))
#define BOOST_NO_OVERLOAD_FOR_NON_CONST
#endif

// gcc-2.95's stringstream is not usable, unless it's the one from STLPORT :
#if BOOST_WORKAROUND(__GNUC__, < 3) && !(defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)) 
#define BOOST_FORMAT_IGNORE_STRINGSTREAM  
#endif


// **** Workaround for io streams, stlport and msvc.
#ifdef BOOST_IO_NEEDS_USING_DECLARATION
namespace boost {
  using std::char_traits;
  using std::basic_ostream;
  using std::basic_ostringstream;
  namespace io {
    using std::basic_ostream;
    namespace detail {
      using std::basic_ios;
      using std::basic_ostream;
      using std::basic_ostringstream;
    }
  }
}
#endif

// ------------------------------------------------------------------------------

#endif // BOOST_FORMAT_MACROS_DEFAULT_HPP
