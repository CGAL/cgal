// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_SHEAR_CONTROLLER
#define CGAL_ACK_SHEAR_CONTROLLER 1

#include <CGAL/basic.h>

#include<set>

CGAL_BEGIN_NAMESPACE

  namespace CGALi {

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
	: max(InitialMax)
#if !AcX_INDEPENDENT_SHEAR_ORDER
           , pos_next_factor(0)
#endif
	{
          CGAL_assertion(max>=1);
#if AcX_STATIC_SEED
	  #warning Warning, uses static seed!
          srand(AcX_STATIC_SEED);
#else
          srand(time(NULL));
#endif
	}
	
        //! Reports that the shear factor \c factor was bad.
	void report_failure(Int factor) {
          this->bad_shears.insert(factor);
	  long failures=static_cast<long>(this->bad_shears.size())+1;
	  if(2*failures>this->max) {
	    this->max*=2;
	  }
	  
	}

        //! Gets a shear factor
        Int get_shear_factor() {
#if !AcX_INDEPENDENT_SHEAR_ORDER
          if(pos_next_factor==static_cast<int>(value_order().size())) {
            value_order().push_back(get_new_shear_factor());
          }
          return value_order()[pos_next_factor++];
#else 
          return get_new_shear_factor();
#endif
        }

    private:

        //! Gets a new shear factor
	Int get_new_shear_factor() {
	  CGAL_assertion(this->bad_shears.size()<max);
	  while(true) {
	    Int s = Int(rand()%max)+1;
	    if(bad_shears.find(s)==bad_shears.end()) {
	      return s;
	    }
	  }
	}

	// Maximal absolute value
	Int max;

#if !AcX_INDEPENDENT_SHEAR_ORDER
        static std::vector<Int>& value_order() {
          static std::vector<Int> value_order_;
          return value_order_;
        }

        int pos_next_factor;
#endif

	
	// Unsuccesfull shear factors
	std::set<Int> bad_shears;
	
    };

  } // namespace CGALi
CGAL_END_NAMESPACE

#endif // CGAL_ACK_SHEAR_CONTROLLER
