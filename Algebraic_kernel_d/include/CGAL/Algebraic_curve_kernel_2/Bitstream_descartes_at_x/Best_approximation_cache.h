// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id$
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_BEST_APPROXIMATION_CACHE
#define CGAL_BEST_APPROXIMATION_CACHE 1


#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <map>
#include <vector>

CGAL_BEGIN_NAMESPACE


  namespace CGALi {

    /*!
     * \brief Stores computed approximations for coefficients.
     *
     * During the execution of the Bitstream Descartes method, approximations
     * of the coefficients are computed independently by different branches of
     * the Descartes tree. To avoid computing them several times, the best
     * konwn approximation for a special coefficient is stored in a \c std::map
     * structure. The data contains the used precision and the integer that
     * represents the coefficient for this precision, according to the 
     * \c CGAL::CGALi::Bitstream_descartes_rndl_tree_traits definition.
     */
    template<typename Coefficient_,typename Integer_>
      class Best_approximation_cache 
      : public ::CGAL::Handle_with_policy<std::map<Coefficient_,std::pair<long,Integer_> > > {

      public:

      //! The Coefficient type
      typedef Coefficient_ Coefficient;

      //! The integer type
      typedef Integer_ Integer;

      typedef std::pair<long,Integer> Data_pair;

      typedef std::map<Coefficient,Data_pair> Map;

      typedef ::CGAL::Handle_with_policy<Map> Base;

      typedef Best_approximation_cache<Coefficient,Integer> Self;

      //! Default constructor
      Best_approximation_cache() {
      }
      
      //! Copy constructor
      Best_approximation_cache(const Self& s)
	: Base(static_cast<const Base&>(s))
	{}
      
      /*! 
       * \brief Checks whether the coefficient \c c already has already been
       * approximated.
       */
      bool approximation_exists(Coefficient c) const {
	typename Map::const_iterator it = this->ptr()->find(c);
	return it!=this->ptr()->end();
      }

      /*!
       * \brief The best approximation known for the coefficient <tt>c</tt>
       * is returned.
       *
       * It is necessary to check whether an approximation of <tt>c</tt> exists
       * using the \c approximation_exists method beforehand!
       */
      void get_best_approximation(Coefficient c,long& prec, Integer& val) const{
	typename Map::const_iterator it = this->ptr()->find(c);
	CGAL_assertion(it!=this->ptr()->end());
	prec = it->second.first;
	val = it->second.second;
      }

      /*!
       * \brief Updates the approximation memory
       *
       * If an approximation for <tt>c</tt> is known, this method can be used
       * to store it for later usage. The update is only performed if the
       * precision is indeed higher than the current precision in the map.
       */
      void update_approximation(Coefficient c,long prec, const Integer& val) {
	typename Map::iterator it = this->ptr()->find(c);
	if(it!=this->ptr()->end()) {
	  if(it->second.first >= prec) {
	    return;
	  }
	  else {
	    Data_pair better_pair = std::make_pair(prec,val);
	    it->second=better_pair;
	  }
	}
	else {
	  Data_pair new_pair = std::make_pair(prec,val);
	  this->ptr()->operator[](c)=new_pair;
	}
	
      }
    };
    

  }// namespace CGALi

CGAL_END_NAMESPACE

#endif // AcX_BEST_APPROXIMATION_CACHE
