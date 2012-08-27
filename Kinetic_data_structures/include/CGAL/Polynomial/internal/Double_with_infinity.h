// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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


#ifndef CGAL_POLYNOMIAL_INTERNAL_WRAPPED_DOUBLE_H
#define CGAL_POLYNOMIAL_INTERNAL_WRAPPED_DOUBLE_H
#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

 template <bool S, bool Q>
  struct Double_with_infinity_default {
    static double value(){return (std::numeric_limits<double>::max)();};
  };
  template<bool O>
  struct Double_with_infinity_default<true,O> {
    static double value(){return std::numeric_limits<double>::signaling_NaN();};
  };
  template<>
  struct Double_with_infinity_default<false,true> {
    static double value(){return std::numeric_limits<double>::quiet_NaN();};
  };

struct Double_with_infinity {
  static double get_default() {
    return Double_with_infinity_default<std::numeric_limits<double>::has_signaling_NaN,
      std::numeric_limits<double>::has_quiet_NaN>::value();
      
  }
  
  Double_with_infinity(): d_(get_default()){}
  Double_with_infinity(double d): d_(d){}
  operator const double&() const {
    return d_;
  }

  bool is_rational() const {
    return true;
  }
  double to_rational() const {return d_;}
  bool is_even_multiplicity() const {
    return false;
  }

  std::pair<double, double> isolating_interval() const {
    return std::make_pair(d_, d_);
  }
  
protected:
  double d_;
};

} } } //namespace CGAL::POLYNOMIAL::internal

namespace std {
  template <>
  class numeric_limits<CGAL_POLYNOMIAL_NS::internal::Double_with_infinity >: public numeric_limits<double>
  {
  public:
    static const bool is_specialized = true;
    static const bool has_infinity=true;
    static double infinity() throw() {return (std::numeric_limits<double>::max)();}
  };
}


namespace CGAL {

inline double to_double(CGAL_POLYNOMIAL_NS::internal::Double_with_infinity d) {
  return to_double(static_cast<double>(d));
}

} //namespace CGAL
  
#endif
