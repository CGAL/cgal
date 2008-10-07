// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

/*!\file include/CGAL/Arrangement_2l/Z_stack_helpers.h
 * \brief definition of Z_stack helpers
 */

#ifndef CGAL_ARRANGEMENT_2l_Z_STACK_HELPERS_H
#define CGAL_ARRANGEMENT_2l_Z_STACK_HELPERS_H 1

#include <CGAL/config.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief
 * Functor to determine the order and overlaps of two realrootisolators
 */
template < class RealRootIsolator >
class Refineable_interval_helper {
public:
    //! this instance's template parameter
    typedef RealRootIsolator Real_root_isolator;
    
    /*!\brief
     * Compute order of two refineable intervals given by
     * pair \c isolator1 and \c i1 and the pair
     * \c isolator2 and \c i2
     *
     * Returns \c CGAL::SMALLER if \c i1 -th interval of \c isolator1
     * is strictly less then the other, CGAL::LARGER, if strictly
     * greater, and CGAL::EQUAL in case of an overlap of the two
     * given intervals.
     */
    ::CGAL::Comparison_result overlap_or_order(
            const Real_root_isolator& isolator1, int i1,
            const Real_root_isolator& isolator2, int i2) const {
        
        if (isolator1.right_boundary(i1) < isolator2.left_boundary(i2)) {
            return CGAL::SMALLER;
        }
        if (isolator2.right_boundary(i2) < isolator1.left_boundary(i1)) {
            return CGAL::LARGER;
        }
        return CGAL::EQUAL;
    }

    /*!\brief
     * Return whether interval defined by  the pair \c isolator1 and \c i1 
     * is included in intervals defined by pair \c isolator2 and \c i2
     *
     * \c strict = true indicates a inclusion in interior
     */
    bool is_included(
            const Real_root_isolator& isolator1, int i1,
            const Real_root_isolator& isolator2, int i2,
            bool strict = false) const {
        
        if (strict) {
            return 
                (isolator1.left_boundary(i1) > isolator2.left_boundary(i2)) &&
                (isolator1.right_boundary(i1) < isolator2.right_boundary(i2));
        } 
        // else 
        return 
            (isolator1.left_boundary(i1) >= isolator2.left_boundary(i2)) &&
            (isolator1.right_boundary(i1) <= isolator2.right_boundary(i2));
    }
    
    
    /*!\brief
     * computs all pairwise refined overlapping id
     *
     * value_type of OutputIterator is std::pair< int, int >
     */
    template < class OutputIterator >
    OutputIterator refined_overlaps(
            const Real_root_isolator& isolator1,
            const Real_root_isolator& isolator2,
            OutputIterator oi
    ) const {
        
        // obtain minimal number of overlaps, i.e., all each interval
        // of one isolator overlaps with at most one interval of the
        // other isolator

        const int num1 = isolator1.number_of_real_roots();
        const int num2 = isolator2.number_of_real_roots();
        CGAL_precondition(num1 > 0);
        CGAL_precondition(num2 > 0);
        
        // TASK cache result for two isolators?

        // compute all overlaps

        std::list< int > empty;

        std::vector< std::list< int > > ovl1(num1, empty);
        std::vector< std::list< int > > ovl2(num2, empty);

        int i1 = 0;
        int i2 = 0;
        while (true) {
            if (i1 == num1) {
                // no further overlaps possible
                break;
            }
            if (i2 == num2) {
                // no further overlaps possible
                break;
            }
            CGAL::Comparison_result res = 
                this->overlap_or_order(isolator1, i1, isolator2, i2);
            
            switch (res) {
            case CGAL::SMALLER:
                //std::cout << "SM" << std::endl;
                i1++;
                break;
            case CGAL::EQUAL:
                //std::cout << "EQUAL " << i1 << " " << i2 << std::endl;
                ovl1[i1].push_back(i2);
                ovl2[i2].push_back(i1);
                if (isolator1.right_boundary(i1) <=
                    isolator2.right_boundary(i2)) {
                    i1++;
                } else if (isolator1.right_boundary(i1) >=
                           isolator2.right_boundary(i2)) {
                    i2++;
                }
                break;
            case CGAL::LARGER:
                //std::cout << "LA" << std::endl;
                i2++;
                break;
            }
        }
        
        for (i1 = 0; i1 < num1; i1++) {
            switch (ovl1[i1].size()) {
            case 0:
                // nothing to do
            case 1:
                // is also fine
                break;
            default:
                // size is greater than 1, refine against all found interval
                do {
                    std::list< typename std::list< int >::iterator > erase;
                    for (typename std::list< int >::iterator it = 
                             ovl1[i1].begin();
                         it != ovl1[i1].end(); it++) {
                        isolator1.refine_interval(i1);
                        isolator2.refine_interval(*it);
                        if (this->overlap_or_order(isolator1, i1, 
                                                   isolator2, *it) != 
                            CGAL::EQUAL) {
                            erase.push_back(it);
                        }
                    }
                    for (typename 
                             std::list< 
                             typename std::list< int >::iterator >::iterator
                             eit = erase.begin(); eit != erase.end(); eit++) {
                        ovl2[*(*eit)].remove(i1);
                        ovl1[i1].erase(*eit);
                    }
                } while (static_cast< int >(ovl1[i1].size()) > 1);
            }
        }

        for (i2 = 0; i2 < num2; i2++) {
            switch (ovl2[i2].size()) {
            case 0:
                // nothing to do
            case 1:
                // is also fine
                break;
            default: 
                // size is greater than 1, refine against all found interval
                do {
                    std::list< typename std::list< int >::iterator > erase;
                    for (typename std::list< int >::iterator it = 
                             ovl2[i2].begin();
                         it != ovl2[i2].end(); it++) {
                        isolator1.refine_interval(*it);
                        isolator2.refine_interval(i2);
                        if (this->overlap_or_order(isolator1, *it, 
                                                   isolator2, i2) != 
                            CGAL::EQUAL) {
                            erase.push_back(it);
                        }
                    }
                    for (typename 
                             std::list< typename std::list< int >::iterator >::
                             iterator
                             eit = erase.begin(); eit != erase.end(); eit++) {
                        ovl1[*(*eit)].remove(i2);
                        ovl2[i2].erase(*eit);
                    }
                } while (static_cast< int >(ovl2[i2].size()) > 1);
            }
        }

        for (i1 = 0; i1 < num1; i1++) {
            switch (ovl1[i1].size()) {
            case 0:
                // nothing to do
                break;
            case 1: {
                // check whether it exist in other
                int i2 = ovl1[i1].front();
                CGAL_assertion(static_cast< int >(ovl2[i2].size()) == 1);
                CGAL_assertion(ovl2[i2].front() == i1);
                *oi++ = std::make_pair(i1,i2);
                break;
            }
            default:
                CGAL_error_msg("should not occur");
            }
        }

        return oi;
    }
    
    //! returns the unique overlapping pair of intervals of two isolators
    template< class InputIterator >
    std::pair< int, int > unique_overlap(
            const Real_root_isolator& isolator1,
            const Real_root_isolator& isolator2,
            InputIterator begin,
            InputIterator end
    ) const {
        
        // obtain the unique overlap by further refinement, as we know
        // from the external that is is exactly one overlap
        std::list< std::pair< int, int > > overlaps;
        std::copy(begin,end, std::back_inserter(overlaps));
        
        CGAL_precondition(static_cast< int >(overlaps.size()) > 0);
        
        while (static_cast< int >(overlaps.size()) > 1) {
            std::list< typename std::list< std::pair< int, int > >::iterator > 
                erase;
            
            // for each overlapping pair
            for (typename std::list< std::pair< int, int > >::iterator 
                     it = overlaps.begin(); it != overlaps.end(); it++) {
                // refine overlapping intervals
                isolator1.refine_interval(it->first);
                isolator2.refine_interval(it->second);

                // and remove non-overlapping pairs
                if (this->overlap_or_order(
                            isolator1, it->first, isolator2, it->second
                    ) != CGAL::EQUAL) {
                    erase.push_back(it);
                }
            }
            for (typename   
                         std::list< typename std::list< 
                         std::pair< int, int > >::iterator >::iterator it =
                         erase.begin(); it != erase.end(); it++) {
                overlaps.erase(*it);
            }
        }
        CGAL_assertion(static_cast< int >(overlaps.size()) == 1);

        CGAL_postcondition(overlaps.begin()->first >= 0);
        CGAL_postcondition(overlaps.begin()->second >= 0);
        return std::make_pair(
                overlaps.begin()->first, overlaps.begin()->second
        );
    }
    
}; // Refineable_interval_helper

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_Z_STACK_HELPERS_H
// EOF
