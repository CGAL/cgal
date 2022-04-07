#ifndef LEXICAL_CAST_HPP
#define LEXICAL_CAST_HPP

/*! This files provides lexical casts from std::string to any one of the number
 * types we intend to use in the benchmark, and a lexical cast does not exist.
 * It is inspired by boost::lexical_cast
 */

#include <CGAL/basic.h>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>

#include "number_type.hpp"

#if 0
template<typename Target, typename Source>
Target lexical_cast(Source arg)
{
  return static_cast<Target>(arg);
}
#endif

template<typename Target>
Target lexical_cast(std::string & arg)
{
  Target result;
  arg >> result;
  return result;
}

#ifdef CGAL_MP_FLOAT_H

template <>
CGAL::MP_Float lexical_cast<CGAL::MP_Float>(std::string & str)
{
  std::stringstream ss(str);
  CGAL::MP_Float result;
  ss >> result;
  return result;
}

template<>
CGAL::Quotient<CGAL::MP_Float>
lexical_cast<CGAL::Quotient<CGAL::MP_Float> >(std::string & str)
{
  typedef CGAL::MP_Float        RT;

  RT rt = lexical_cast<RT>(str);
  CGAL::Quotient<RT> result(rt, 1);
  return result;
}

#endif

#ifdef LEDA_INTEGER_H

template <> leda::integer lexical_cast<leda::integer>(std::string & str)
{
  leda::integer result(str.c_str());
  return result;
}

#endif

#ifdef LEDA_RATIONAL_H

template <> leda::rational lexical_cast<leda::rational>(std::string & str)
{
  leda::rational result(str.c_str());
  return result;
}

#endif

#ifdef CGAL_GMPZ_H

template <> CGAL::Gmpz lexical_cast<CGAL::Gmpz>(std::string & str)
{
  CGAL::Gmpz result(str);
  return result;
}

template<>
CGAL::Quotient<CGAL::Gmpz>
lexical_cast<CGAL::Quotient<CGAL::Gmpz> >(std::string & str)
{
  typedef CGAL::Gmpz            RT;

  RT rt = lexical_cast<RT>(str);
  CGAL::Quotient<RT> result(rt, 1);
  return result;
}

#endif

#ifdef CGAL_GMPQ_H

template <> CGAL::Gmpq lexical_cast<CGAL::Gmpq>(std::string & str)
{
  CGAL::Gmpq result(str);
  return result;
}

#endif

#ifdef DOUBLE_HPP

template <> CGAL::Double lexical_cast<CGAL::Double>(std::string & str)
{
  double d = boost::lexical_cast<double>(str);
  return CGAL::Double(d);
}

#endif

#endif
