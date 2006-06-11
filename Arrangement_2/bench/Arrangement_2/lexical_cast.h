#ifndef LEXICAL_CAST_H
#define LEXICAL_CAST_H

#include <CGAL/basic.h>
#include <string>
#include <assert.h>
#include <boost/lexical_cast.hpp>

#include "number_type.h"

template<typename Target, typename Source>
Target lexical_cast(Source arg)
{
  Target result;
  std::cerr << "throw an exeption!\n";
  assert(0);
  return result;
}


#ifdef CGAL_MP_FLOAT_H

template<>
CGAL::MP_Float lexical_cast<CGAL::MP_Float, std::string&>(std::string & str)
{
  float f = boost::lexical_cast<float>(str);
  CGAL::MP_Float result = f;
  return result;
}

template<>
CGAL::Quotient<CGAL::MP_Float> lexical_cast<CGAL::Quotient<CGAL::MP_Float>, std::string&>(std::string & str)
{
  float f = boost::lexical_cast<float>(str);
  CGAL::MP_Float mpf = f;
  CGAL::Quotient<CGAL::MP_Float> result(mpf, 1);
  return result;
}

#endif

#endif
