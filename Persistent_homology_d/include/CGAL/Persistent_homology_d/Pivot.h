// Copyright (c) 2003,2004  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Lutz Kettner, Afra Zomorodian

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_PIVOT_H
#define CGAL_PERSISTENT_HOMOLOGY_D_PIVOT_H

#include <iostream>   // std::cout
#include <vector>
#include <functional> // std::less


namespace CGAL {
namespace Persistent_homology_d {


// Class:  Pivot
// *************
// Note for IsPositive to make any sense, the pivots must
// be declared in an array
template <class SimplexHandle, class OutputPolicy, class CoefficientPolicy>
class Pivot {
public:
  // members
  typedef SimplexHandle Simplex_handle;
  typedef OutputPolicy  Output_policy;
  typedef typename std::vector< Pivot< Simplex_handle, OutputPolicy, 
	  CoefficientPolicy> >::iterator               Pivot_handle;
  typedef std::less< Pivot_handle>                     CompareWrtoFiltration;
  // define chain, boundary (and extended versions)
  typedef typename CoefficientPolicy::template 
    Bind< Pivot_handle, CompareWrtoFiltration>::Chain     Chain;
  typedef typename CoefficientPolicy::template 
    Bind< Pivot_handle, CompareWrtoFiltration>::Ext_Chain Ext_Chain;
  typedef typename Output_policy::template 
    Bind< Chain>::Boundary                                Boundary;
  typedef typename Output_policy::template 
    Bind< Ext_Chain>::Boundary                            Ext_Boundary;

protected:
  static Pivot_handle base_;          // needed for printing purposes
  Pivot_handle   paired_with_;
  Boundary       boundary_;
  Simplex_handle simplex_;            // need to store to return answer  
public:
  const Pivot_handle&   paired_with() const { return paired_with_; }
  Pivot_handle&         paired_with()       { return paired_with_; }
  const Boundary&       boundary()    const { return boundary_;    }
  Boundary&             boundary()          { return boundary_;    }
  const Simplex_handle& simplex()     const { return simplex_;     }
  Simplex_handle&       simplex()           { return simplex_;     }

  bool is_positive()                  const { return this < &(*paired_with_);}
  void set_base(Pivot_handle base)          { base_ = base;        }
  void print()                        const { 
    std::cout << this - &(*base_);
  }
};

// Definition of static variable
template< class SimplexHandle, class OutputPolicy, class CoefficientPolicy>
typename Pivot< SimplexHandle, OutputPolicy, CoefficientPolicy>::Pivot_handle
         Pivot< SimplexHandle, OutputPolicy, CoefficientPolicy>::base_;

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_PIVOT_H
