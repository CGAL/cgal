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

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_CHAIN_H
#define CGAL_PERSISTENT_HOMOLOGY_D_CHAIN_H

#include <vector>
#include <functional>

namespace CGAL {
namespace Persistent_homology_d {


// Forward declaration (necessary for operator=)
template <class PivotHandle, class CompareWrtoFiltration>
class Ext_Chain_Z2;

// Class:  Chain_Z2
// ****************
// A chain with Z2 coefficients.  The coefficients are not
// stored.  Rather, a pivothandle being in the chain means
// that it has coefficient 1.

template <class PivotHandle, 
          class CompareWrtoFiltration = std::less< PivotHandle> >
class Chain_Z2 :
  public std::vector< PivotHandle> 
{
public:
  typedef Chain_Z2< PivotHandle, CompareWrtoFiltration> Self;
  typedef std::vector< PivotHandle>       Vector;
  typedef typename Vector::size_type      size_type;
  typedef int                             Field;
  typedef PivotHandle                     Pivot_handle;
  typedef CompareWrtoFiltration           Compare_wrto_filtration;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;

  Chain_Z2() {}
  Chain_Z2(size_type n) : Vector( n) {}
  Self& operator=(const Ext_Chain_Z2
		  < Pivot_handle, CompareWrtoFiltration> &c) {
    *(static_cast<Vector*>(this)) = Vector( c.begin(), c.end());
    return *this;
  }
  // coeff, handle accessors
  // no reference type for coefficient as you cannot set it!
  Pivot_handle handle        (iterator i)        const { return *i;   }
  const Field& coeff         (const_iterator i)  const { return one;  }
  const Pivot_handle handle  (const_iterator i)  const { return *i;   }

  void print() const
  {
    for(const_iterator i = begin(); i != end(); i++) {
      handle(i)->print();
      std::cout << ' ';
    }
  }
};

// Specialization so that vectors aren't copied
template <class PivotHandle>
inline void swap(Chain_Z2< PivotHandle> &x, Chain_Z2< PivotHandle> &y)
{
  x.swap(y);
}

// Incomplete Definition
template <class PivotHandle,
          class Field_,
          class CompareWrtoFiltration>
class Ext_Chain_field;

// Class:  Chain_field
// *******************
// Chains of PivotHandles with coefficients in Field.
// Coefficients are stored explicitly.

template <class PivotHandle,
          class Field_,
          class CompareWrtoFiltration = std::less< PivotHandle> >
class Chain_field : 
  public std::vector< std::pair< PivotHandle, Field_ > >
{
public:
  typedef Chain_field< PivotHandle, Field_, CompareWrtoFiltration> Self;
  typedef std::vector< std::pair< PivotHandle, Field_ > >          Vector;
  typedef typename Vector::size_type      size_type;
  typedef Field_                          Field;
  typedef PivotHandle                     Pivot_handle;
  typedef CompareWrtoFiltration           Compare_wrto_filtration;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;

  Chain_field() {}
  Chain_field(size_type n) : Vector( n) {}
  Self& operator=(const Ext_Chain_field< Pivot_handle, Field_, 
		  CompareWrtoFiltration> &c) {
    *(static_cast<Vector*>(this)) = Vector( c.begin(), c.end());
    return *this;
  }

  // coeff, handle accessors
  Field& coeff              (iterator i)       const { return i->second; }
  PivotHandle handle        (iterator i)       const { return i->first;  }
  const Field& coeff        (const_iterator i) const { return i->second; }
  const PivotHandle handle  (const_iterator i) const { return i->first;  }

  // print
  void print() const
  {
    for(const_iterator i = begin(); i != end(); i++) {
      std::cout << "(" << coeff(i) << ") ";
      handle(i)->print();
      std::cout << ' ';
    }
  }
};

// Specialization so that vectors aren't copied
template <class PivotHandle, class Field_>
inline void swap(Chain_field< PivotHandle, Field_> &x, 
		 Chain_field< PivotHandle, Field_> &y)
{
  x.swap(y);
}

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_CHAIN_H
