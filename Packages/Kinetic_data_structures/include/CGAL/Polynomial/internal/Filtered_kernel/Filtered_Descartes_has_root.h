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

#ifndef CGAL_POLYNOMIAL_FILTERED_DESCARTES_HAS_ROOT_H
#define CGAL_POLYNOMIAL_FILTERED_DESCARTES_HAS_ROOT_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_Descartes_root_counter.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE;





template <class Kernel>
class Filtered_Descartes_has_root {
public:
  Filtered_Descartes_has_root(){}

  Filtered_Descartes_has_root(const typename Kernel::Function &fh, Kernel k= Kernel()): h_(fh), kernel_(k) {
  }

  ~Filtered_Descartes_has_root(){
  }

  typedef bool result_type;

  template <class NTT>
  result_type operator()(const NTT &begin, const NTT &end,
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO, 
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO) const {
    return filtered_Descartes_root_counter(h_, begin, end, false, kernel_)!= 0;
  }

protected:
  typename Kernel::Function h_;
  Kernel kernel_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE;

#endif
