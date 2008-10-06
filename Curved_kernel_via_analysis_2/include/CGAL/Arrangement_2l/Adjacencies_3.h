// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
// File          : include/SoX/GAPS/Adjacencies_3.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef SoX_GAPS_ADJACENCIES_3_H
#define SoX_GAPS_ADJACENCIES_3_H 1

/*! \file SoX/GAPS/Adjacencies_3
    \brief Contains Adjacencies_3 class
*/

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

#include <iostream>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

class Adjacencies_3_rep {

public:

    //! the class itself
    typedef Adjacencies_3_rep Self;

    //! a pair of adjacencies
    typedef std::pair< int, int > Adjacency_pair;
    
    //! an interval 
    typedef std::pair< int, int > Adjacency_interval;
    
    //! an adjacency vector
    typedef std::vector < Adjacency_pair > Adjacency_vector;

    //! Default constructor
    Adjacencies_3_rep() {
    }
    
    //! Constructor with data
    Adjacencies_3_rep(Adjacency_vector data) : 
        _m_data(data) {
    }
    
    Adjacency_vector _m_data;

};

} // namespace CGALi

//! The Adjacency output object
class Adjacencies_3 : 
    public CGAL::Handle_with_policy< CGAL::CGALi::Adjacencies_3_rep > {
    
public:
    //! the class itself
    typedef Adjacencies_3 Self;
    
    //! the rep type
    typedef CGALi::Adjacencies_3_rep Rep;

    //! the base type
    typedef CGAL::Handle_with_policy< Rep > Base;
    
    //! a pair of adjacencies
    typedef std::pair< int, int > Adjacency_pair;
    
    //! an interval 
    typedef std::pair< int, int > Adjacency_interval;
    
    //! a vector of adjacency
    typedef std::vector < Adjacency_pair > Adjacency_vector;
    
    //! const version 
    typedef Adjacency_vector::const_iterator Const_adjacency_iterator;
    
public:

    //! Default constructor
    Adjacencies_3() : 
        Base(Rep()) {
    }

    //! Constructor with data
    Adjacencies_3(Adjacency_vector data) :
        Base(Rep(data)) {
    }
    
    //! Returns an interval [i,j] such that all sheets between 
    //! i and j are adjacenct to index
    Adjacency_interval interval(int index) const {
        
        std::vector<int> adjacent_features;
        
        for (int i = 0; 
             i < static_cast<int>(this->ptr()->_m_data.size()); 
             i++) {
            if (this->ptr()->_m_data[i].first == index) {
                adjacent_features.push_back(this->ptr()->_m_data[i].second);
            }
        }
        
        std::sort(adjacent_features.begin(), adjacent_features.end());
        
        int n = static_cast<int>(adjacent_features.size());
        
        if (n == 0) {
            return std::make_pair(1,0);
        }
        // Make sure that the adjacent features form an interval
        for (int i = 0; i < n-1; i++) {
            
            CGAL_assertion( adjacent_features[i] == 
                            adjacent_features[i+1]-1 );
            
        }
        
        return std::make_pair(adjacent_features[0],
                              adjacent_features[n-1]);
    }
    
    //! beginning of adjacencies
    Const_adjacency_iterator begin() const {
        return this->ptr()->_m_data.begin();
    }
    
    //! end of adjacencies
    Const_adjacency_iterator end() const {
        return this->ptr()->_m_data.end();
    }

    //! The number of adjacency pairs
    unsigned int size() const {
        return this->ptr()->_m_data.size();
    }
    
    //! print (for testing only)
    void print() const {
        for(int i = 0; 
            i < static_cast<int>(this->ptr()->_m_data.size()); 
            i++) {
            std::cout << "(" << this->ptr()->_m_data[i].first << ", " 
                      << this->ptr()->_m_data[i].second << ")" << std::endl;
        }
    }
    
    //! swap two instances
    Self swap() const {
        Self copy;
        for(Const_adjacency_iterator it = begin(); it!=end(); it++) {
            copy.ptr()->_m_data.push_back
                (std::make_pair(it->second,it->first));
        }
        return copy;
    }

    //! output operator
    void pretty_print(std::ostream& os) const {
        for (int i = 0; 
             i < static_cast<int>(this->ptr()->_m_data.size()); 
             i++) {
            os << "<" << this->ptr()->_m_data[i].first
               << "," << this->ptr()->_m_data[i].second << ">";
            if (i < static_cast<int>(this->ptr()->_m_data.size()) - 1) {
                os << ",";
            }
        }
    }
        
}; // end of Adjacencies data structure

/*!\relates Adjacencies_3
 * \brief
 * outputs Adjacencies_3 object to stream \c os
 */
inline
std::ostream& operator<<(
        std::ostream& os, 
        const Adjacencies_3& adj) {
    adj.pretty_print(os);
    return os;
}

CGAL_END_NAMESPACE

#endif // SoX_GAPS_ADJACENCIES_3
