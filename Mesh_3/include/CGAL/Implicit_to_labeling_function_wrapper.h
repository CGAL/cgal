// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stéphane Tayeb, Aymeric PELLE
//
//******************************************************************************
// File Description :
// Implicit_to_labeling_function_wrapper and
// Implicit_vector_to_labeling_function_wrapper class declaration
// and implementation.
//
// See classes description to have more information.
//******************************************************************************

#ifndef CGAL_IMPLICIT_TO_LABELING_FUNCTION_WRAPPER_H
#define CGAL_IMPLICIT_TO_LABELING_FUNCTION_WRAPPER_H


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

#include <boost/dynamic_bitset.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>

#include <CGAL/config.h>
#include <CGAL/assertions.h>

namespace CGAL {

/**
 * @class Implicit_to_labeling_function_wrapper
 *
 * This class is designed to wrap an implicit function which describes a domain
 * by [p is inside if f(p)<0] to a function which takes its values into {0,1}.
 * f(p)=0 means that p is outside the domain.
 */
template<class Function_, class BGT>
class Implicit_to_labeling_function_wrapper
{
public:
  // Types
  typedef int                     return_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Implicit_to_labeling_function_wrapper(const Function_& f)
    : r_f_(f) {}

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Implicit_to_labeling_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, const bool = true) const
  {
    return ( (r_f_(p)<0) ? 1 : 0 );
  }

private:
  /// Function to wrap
  const Function_& r_f_;

};  // end class Implicit_to_labeling_function_wrapper



/**
 * \deprecated
 *
 * @class Implicit_vector_to_labeling_function_wrapper
 *
 * Wraps a set of implicit function [f1,f2,...] to one function F which
 * takes its values into N.
 *
 * Let p be a point.
 * F(p) = 0b000000(f2(p)<0)(f1(p)<0)
 *
 * It can handle at most 8 functions.
 */
template<class Function_, class BGT>
class Implicit_vector_to_labeling_function_wrapper
{
public:
  // Types
  typedef int                       return_type;
  typedef std::vector<Function_*>   Function_vector;
  typedef typename BGT::Point_3     Point_3;

  /// Constructor
  Implicit_vector_to_labeling_function_wrapper(const std::vector<Function_*>& v)
    : function_vector_(v) {}

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Implicit_vector_to_labeling_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, const bool = true) const
  {
    int nb_func = function_vector_.size();
    if ( nb_func > 8 )
    {
      CGAL_error_msg("We support at most 8 functions !");
    }

    char bits = 0;
    for ( int i = 0 ; i < nb_func ; ++i )
    {
      // Insert value into bits : we compute fi(p) and insert result at
      // bit i of bits
      bits |= ( ((*function_vector_[i])(p) < 0) << i );
    }

    return ( static_cast<return_type>(bits) );
  }

private:
  /// Functions to wrap
  const Function_vector function_vector_;

};  // end class Implicit_to_labeling_function_wrapper

template <class ImplicitFunction>
class Implicit_multi_domain_to_labeling_function_wrapper
{
  template <class T_>
  class Implicit_function_traits
  {
  public:
    typedef typename T_::Point Point;
  };

  template <class RT_, class Point_>
  class Implicit_function_traits<RT_ (*)(Point_)>
  {
  public:
    typedef typename boost::remove_reference<
            typename boost::remove_cv< Point_ >::type>::type Point;
  };

public:
  typedef int                     return_type;
  typedef ImplicitFunction        Function;
  typedef typename Implicit_function_traits<ImplicitFunction>::Point   Point_3;
  typedef std::vector<Function>   Function_vector;

private:
  std::vector<Function> funcs;
  typedef boost::dynamic_bitset<std::size_t> Bmask;
  std::vector<Bmask> bmasks;

public:
  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& vf, const std::vector<std::vector<Sign> >& vps)
  : funcs(vf), bmasks(vps.size(), Bmask(funcs.size() * 2, false))
  {
    CGAL_assertion(funcs.size() != 0);

    std::size_t mask_index = 0;
    for (std::vector<std::vector<Sign> >::const_iterator mask_iter = vps.begin(), mask_end_iter = vps.end();
         mask_iter != mask_end_iter;
         ++mask_iter)
    {
      const std::vector<Sign>& mask = *mask_iter;
      CGAL_assertion(funcs.size() == mask.size());
      Bmask& bmask = bmasks[mask_index++];

      typename Bmask::size_type bit_index = 0;
      for (std::vector<Sign>::const_iterator iter = mask.begin(), endIter = mask.end(); iter != endIter; ++iter)
      {
        std::string::value_type character = *iter;
        CGAL_assertion(character == POSITIVE || character == NEGATIVE);

        bmask[bit_index] = (character == POSITIVE);
        ++bit_index;
        bmask[bit_index] = (character == NEGATIVE);
        ++bit_index;
      }
    }
    std::sort(bmasks.begin(), bmasks.end());
  }

  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& vf)
  : funcs(vf)
  {
    CGAL_assertion(funcs.size() != 0);

    bmasks.reserve((1 << funcs.size()) - 1);
    bmasks.push_back(Bmask(std::string("10")));
    bmasks.push_back(Bmask(std::string("01")));

    for (std::size_t i = 0; i < funcs.size()-1; ++i)
    {
      std::size_t c_size = bmasks.size();
      for (std::size_t index = 0; index < c_size; ++index)
      {
        Bmask aux = bmasks[index];
        aux.push_back(true);
        aux.push_back(false);
        bmasks.push_back(aux);
        bmasks[index].push_back(false);
        bmasks[index].push_back(true);
      }
    }
    bmasks.pop_back();
    std::sort(bmasks.begin(), bmasks.end());
  }

  Implicit_multi_domain_to_labeling_function_wrapper (const Function_vector& vf, const std::vector<std::string>& vps)
  : funcs(vf), bmasks(vps.size(), Bmask(funcs.size() * 2, false))
  {
    CGAL_assertion(funcs.size() != 0);

    std::size_t mask_index = 0;
    for (std::vector<std::string>::const_iterator mask_iter = vps.begin(), mask_end_iter = vps.end();
         mask_iter != mask_end_iter;
         ++mask_iter)
    {
      const std::string& mask_str = *mask_iter;
      CGAL_assertion(funcs.size() == mask_str.length());
      Bmask& bmask = bmasks[mask_index++];

      typename Bmask::size_type bit_index = 0;
      for (std::string::const_iterator iter = mask_str.begin(), endIter = mask_str.end(); iter != endIter; ++iter)
      {
        std::string::value_type character = *iter;
        CGAL_assertion(character == '+' || character == '-');

        bmask[bit_index] = (character == '+');
        ++bit_index;
        bmask[bit_index] = (character == '-');
        ++bit_index;
      }
    }
    std::sort(bmasks.begin(), bmasks.end());
  }

  return_type operator() (const Point_3& p, const bool = true) const
  {
    Bmask bmask(funcs.size() * 2, false);

    std::size_t i = 0;
    for (typename std::vector<Function>::const_iterator iter = funcs.begin(), endIter = funcs.end();
         iter != endIter;
         ++iter)
    {
      const Function& function = *iter;

      double fres = function(p);
      bmask[i] = fres > 0;
      ++i;
      bmask[i] = fres < 0;
      ++i;
    }

    std::vector<Bmask>::const_iterator iter = std::lower_bound(bmasks.begin(), bmasks.end(), bmask);
    if (iter != bmasks.end() && *iter == bmask)
      return static_cast<return_type>(1 + (iter - bmasks.begin()));
    return 0;
  }
};

}  // end namespace CGAL



#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_IMPLICIT_TO_LABELING_FUNCTION_WRAPPER_H
