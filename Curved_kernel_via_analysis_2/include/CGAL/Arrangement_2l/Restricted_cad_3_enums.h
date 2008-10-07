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

#ifndef CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ENUMS_H
#define CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ENUMS_H 1

/*!\file include/CGAL/Arrangement_2l/Restricted_cad_3_enums.h
 * \brief Contains enumeration related to  
 * \link Restricted_cad_3 \endlink
 */

#include <CGAL/config.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE

//! distinguishes between the three possible dcel-features
enum Dcel_feature {
    VERTEX = 0, //!< identifies a vertex
    EDGE = 1, //!< identifies an edge
    FACE = 2 //!< identifies a face
};

/*!\relates Dcel_feature
   \brief output operator
   
   The output is intended for debugging purposes only and hence always
   is in a human-readable pretty-print format.
   
   
*/
// \pre ::CGAL::is_pretty(os)
inline std::ostream& operator<< (std::ostream& os, Dcel_feature feat) {
    //CGAL_precondition(::CGAL::is_pretty(os));
    static const char* names[] = { "VERTEX", "EDGE", "FACE" };
    
    CGAL_assertion(feat >= 0 && 
                   feat < static_cast<int>(sizeof names / sizeof *names));
    return os << names[feat];
}


//! container for integers
class Nk {
public:

    enum Value_type {
        MULT = 1,
        N = 2,
        K = 3
    };
    
    Nk() :
        _m_feat1(CGAL::VERTEX),
        _m_feat2(CGAL::VERTEX),
        _m_mult(-1),
        // TODO other initial value?
        _m_n(-2),
        _m_k(-2),
        _m_fixed(false) {
    }
    
public:
    
    CGAL::Dcel_feature feature1() const {
        return _m_feat1;
    }

    CGAL::Dcel_feature feature2() const {
        return _m_feat2;
    }

    void set_feature1(CGAL::Dcel_feature feature) const {
        CGAL_precondition(!_m_fixed);
        _m_feat1 = feature;
    }
    
    void set_feature2(CGAL::Dcel_feature feature) const {
        CGAL_precondition(!_m_fixed);
        _m_feat2 = feature;
    }
    
    int mult() const {
        return _m_mult;
    }
    
    int n() const {
        return _m_n;
    }
    
    int k() const {
        return _m_k;
    }
    
    void set_mult(int mult) const {
        CGAL_precondition(mult >= 0);
        CGAL_precondition(!_m_fixed);
        _m_mult = mult;
    }
    
    void set_n(int n) const {
        CGAL_precondition(n >= -1);
        CGAL_precondition(!_m_fixed);
        _m_n = n;
    }
    
    void set_k(int k) const {
        CGAL_precondition(k >= -1);
        CGAL_precondition(!_m_fixed);
        _m_k = k;
    }
    
    void fix() const {
        _m_fixed = true;
    }

    bool same_nk(const Nk& nk) const {
        return 
            (this->_m_n == nk._m_n) &&
            (this->_m_k == nk._m_k);
    }

private:
    // members
    mutable CGAL::Dcel_feature _m_feat1;
    mutable CGAL::Dcel_feature _m_feat2;


    mutable int _m_mult;
    mutable int _m_n;
    mutable int _m_k;
    
    mutable bool _m_fixed;
};


inline
std::ostream& operator<<(std::ostream& out, const CGAL::Nk& nk) {
    out << "NK(mult=" << nk.mult() << ",n=" << nk.n() 
        << ",k=" << nk.k() << ")";
    return out;
}

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_RESTRICTED_CAD_3_ENUMS_H
// EOF
