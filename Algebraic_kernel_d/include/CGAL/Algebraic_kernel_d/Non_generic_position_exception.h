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



#ifndef CGAL_NON_GENERIC_POSITION_EXCEPTION
#define CGAL_NON_GENERIC_POSITION_EXCEPTION


CGAL_BEGIN_NAMESPACE

  namespace CGALi {

    /*! 
     * \brief Exception class for not sufficiently generic positions.
     *
     * Must be thrown whenever a curve cannot be analysed because its position
     * is not "good enough".
     */
    class Non_generic_position_exception {
      
    public:

      //! Default constructible
      Non_generic_position_exception() {}
      
    };

  } // namespace CGALi

CGAL_END_NAMESPACE


#endif
