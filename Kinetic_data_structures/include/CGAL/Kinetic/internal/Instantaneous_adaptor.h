// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_INSTANTANEOUS_ADAPTOR_H
#define CGAL_INSTANTANEOUS_ADAPTOR_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/functional_base.h>

CGAL_KINETIC_BEGIN_NAMESPACE;

//! An object to help convert between moving objects and their static representations to wrap a predicate.
/*!  A ref counted pointer is stored to the part of the
  Instantaneous_kernel which actually performs the conversions.

  Look at the source to figure out how it works.
*/
template <class Predicate, class Curcoord, class Object>
class Instantaneous_adaptor
{
public:
  Instantaneous_adaptor(Curcoord ik,
			Predicate pred=Predicate()): cc_(ik), pred_(pred) {
  }

  typedef typename Predicate::result_type result_type;
  typedef Object argument_type;
  typedef argument_type first_argument_type;
  typedef argument_type second_argument_type;
  typedef argument_type third_argument_type;
  typedef argument_type fourth_argument_type;
  typedef argument_type fifth_argument_type;
  typedef typename Arity_traits<Predicate>::Arity Arity;

  result_type operator()(const first_argument_type &arg0) const
  {
    return pred_(cc_(arg0));
  }

  result_type operator()(const first_argument_type &arg0,
			 const second_argument_type &arg1) const
  {
    /*std::cout << "Args " << cc_(arg0) <<", "
      << cc_(arg1) << " result " << pred_(cc_(arg0), cc_(arg1))
      << " antiresult " << pred_(cc_(arg1), cc_(arg0)) << std::endl;*/

    return pred_(cc_(arg0), cc_(arg1));
  }

  result_type operator()(const first_argument_type &arg0,
			 const second_argument_type &arg1,
			 const third_argument_type &arg2) const
  {
    return pred_(cc_(arg0), cc_(arg1),
		 cc_(arg2));
  }

  result_type operator()(const first_argument_type &arg0,
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3) const
  {
    return pred_(cc_(arg0), cc_(arg1),
		 cc_(arg2), cc_(arg3));
  }

  result_type operator()(const first_argument_type &arg0,
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3,
			 const fifth_argument_type &arg4) const
  {
    return pred_(cc_(arg0), cc_(arg1),
		 cc_(arg2), cc_(arg3),
		 cc_(arg4));
  }

protected:

  Curcoord cc_;
  Predicate pred_;
};

CGAL_KINETIC_END_NAMESPACE;
#endif
