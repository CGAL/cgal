// (C) Copyright Gennaro Prota 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.


#ifndef BOOST_NON_TYPE_HPP_GP_20030417
#define BOOST_NON_TYPE_HPP_GP_20030417


namespace boost {

  // Just a simple "envelope" for non-type template parameters. Useful
  // to work around some MSVC deficiencies.

 template <typename T, T n>
 struct non_type { };


}


#endif // include guard
