// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_SHEAR_CONTROLLER
#define CGAL_ACK_SHEAR_CONTROLLER 1

#include <CGAL/basic.h>
#include <CGAL/tss.h>

#include<set>

namespace CGAL {

  namespace internal {

    /*!
     * \brief A class that controls the used shear factors
     *
     * The objects returns positive integers that are used as shear factors.
     * It choses integers from the range \c 1..max at random , 
     * where \c c is a positive integer
     * initially set in the constructor (8 by default). No integer is given
     * twice by \c get_shear_factor(), at least if the failed ones are reported
     * with the \c report_failure() method.
     * If more than half of the integers in the range were bad, the range is
     * enlarged to \c 1..2*max.
     */
    template<typename Int,int InitialMax=8>
      class Shear_controller {
      
    public:
      
      //! Constructor, getting the maximal absolute value of the shear factor
      Shear_controller()
	: m_max(InitialMax)
        , pos_next_factor(0)  {
          CGAL_assertion(m_max>=1);
#if CGAL_ACK_STATIC_SEED
	  #warning Warning, uses static seed!
          srand(CGAL_ACK_STATIC_SEED);
#else
          srand(time(NULL));
#endif
	}
	
        //! Reports that the shear factor \c factor was bad.
	void report_failure(Int factor) {
          this->bad_shears.insert(factor);
	  long failures=static_cast<long>(this->bad_shears.size())+1;
	  if(2*failures>this->m_max) {
	    this->m_max*=2;
	  }
	  
	}

        //! Gets a shear factor
        Int get_shear_factor() {
          if(pos_next_factor==static_cast<int>(value_order().size())) {
            value_order().push_back(get_new_shear_factor());
          }
          return value_order()[pos_next_factor++];
        }

    private:

        //! Gets a new shear factor
	Int get_new_shear_factor() {
	  CGAL_assertion(int(this->bad_shears.size())<m_max);
	  while(true) {
	    Int s = Int(rand()%m_max)+1;
	    if(bad_shears.find(s)==bad_shears.end()) {
	      return s;
	    }
	  }
	}

	// Maximal absolute value
	Int m_max;

        static std::vector<Int>& value_order() {
          CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(std::vector<Int>, value_order_);
          return value_order_;
        }

        int pos_next_factor;
	
	// Unsuccessfull shear factors
	std::set<Int> bad_shears;
	
    };

  } // namespace internal
} //namespace CGAL

#endif // CGAL_ACK_SHEAR_CONTROLLER
