// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_DESCARTES_HAS_ROOT_H
#define CGAL_POLYNOMIAL_DESCARTES_HAS_ROOT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Tools/Alternation_counter.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE;



template <class Kernel>
class Descartes_has_root {
public:
  Descartes_has_root(const typename Kernel::Function &f,
			 Kernel k= Kernel()): map_(k.map_rational_interval_to_positive_object(f)){
    
  }

  typedef bool result_type;
  typedef typename Kernel::Function::NT first_argument_type;
  typedef typename Kernel::Function::NT second_argument_type;

  //! Note, the result is an upper bound
  /*!
    Bezier root counter needs the sign version of the function. This is added to make the two compatible.
  */
  template <class NTT>
  result_type operator()(const NTT &lb, const NTT &ub,
			 POLYNOMIAL_NS::Sign, POLYNOMIAL_NS::Sign) const {
    return operator()(lb,ub);
  }
  template <class NTT>
  result_type operator()(const NTT &lb, const NTT &ub) const {
    typename Kernel::Function mf= map_(lb, ub);
    
    typename POLYNOMIAL_NS::Alternation_counter<first_argument_type> ac;
    for (int i=0; i<= mf.degree(); ++i){
      ac.push_back(mf[i]);
      if (ac.number_of_alternations() >0) return true;
      
    }
    //std::cout << "Num alternations is " << ac.number_of_alternations() << std::endl;
    return false;
  }
protected:
  typename Kernel::Map_rational_interval_to_positive map_;
};
CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE;

#endif
