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

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_BOUNDARY_H
#define CGAL_PERSISTENT_HOMOLOGY_D_BOUNDARY_H

#include <iostream> // std::cout
#include <utility>  // std::pair

namespace CGAL {
namespace Persistent_homology_d {


// Class:  Boundary
// ****************
// Contains a boundary and allows for pivot elimination
// using second boundary.  The latter must have the same
// pivot and have pivot coefficient equal to 1.

template <class Chain>
class Boundary {
protected:
  Chain boundary_;
public:
  typedef Boundary< Chain>                         Self;
  typedef typename Chain::Field                    Field;
  typedef typename Chain::Pivot_handle             Pivot_handle;
  typedef typename Chain::const_iterator           const_iterator;
  typedef typename Chain::iterator                 iterator;
  typedef typename Chain::size_type                size_type;
  typedef typename Chain::Compare_wrto_filtration  Compare_wrto_filtration;
  typedef typename Chain::reference                reference;
  typedef typename Chain::const_reference          const_reference;

  // Constructors
  Boundary() {}
  Boundary(size_type n) : boundary_(n) {}
  Boundary(Chain& c)    : boundary_(c) {}

  // operator=
  template< class Chain2>
  Self& operator=(const Boundary< Chain2> &b)
  {
    boundary_ = b.boundary();
    return *this;
  }

  // extend some chain operations out
  const Chain&    boundary()  const { return boundary_;           }  
  void            clear()           { boundary_.clear();          }
  bool            empty()     const { return boundary_.empty();   }
  size_type       size()      const { return boundary_.size();    }
  iterator        pivot()           { return boundary_.end()-1;   }
  const_iterator  pivot()     const { return boundary_.end()-1;   }
  void            print()     const { boundary_.print();          }
  void            push_back(const std::pair< Pivot_handle, Field> &p)
  { 
      boundary_.push_back( p); 
  }
  void            pop_pivot()       { boundary_.pop_back();       }
  
  // eliminate_pivot:  
  // boundary2 must have same pivot and coefficient equal to 1.
  // we assume the pivot is not stored in boundary2, but is
  // rather implied.  therefore, we eliminate the pivot in 
  // boundary at the end to maintain this invariant.
  template< class Boundary2>
  void eliminate_pivot( const Boundary2& boundary2,
			Compare_wrto_filtration handle_compare = 
			Compare_wrto_filtration())
  {
    // me = me - mypivotcoeff * boundary2 (which has pivotcoeff = 1)
    const typename Chain::Field& pivot_coeff = boundary_.coeff(pivot());
    boundary_.multiply_and_add( boundary2.boundary(), -pivot_coeff, 
				handle_compare);
    // delete pivot
    // note:  this would've happened in multiply_and_add if we had
    // explicit pivots
    pop_pivot();
  }
  
  void normalize() 
  {
    if( empty()) return;
    const Field& pivot_coeff = boundary_.coeff(pivot());
    if( pivot_coeff != Chain::one) {
      Field inverse = Chain::one / pivot_coeff;
      boundary_.multiply( inverse);
    }
  }
};

// Class:  Boundary_with_manifold
// ******************************
// Also contains the immersed manifold whose boundary
// it stores.

template <class Chain>
class Boundary_with_manifold : public Boundary< Chain> {
 protected:
  Chain manifold_;
public:
  typedef typename Chain::Field                    Field;
  typedef typename Chain::Pivot_handle             Pivot_handle;
  typedef typename Chain::const_iterator           const_iterator;
  typedef typename Chain::iterator                 iterator;
  typedef typename Chain::size_type                size_type;
  typedef typename Chain::Compare_wrto_filtration  Compare_wrto_filtration;
  typedef typename Chain::reference                reference;
  typedef typename Chain::const_reference          const_reference;

  // Constructor from chains
  Boundary_with_manifold( ) {}
  Boundary_with_manifold( size_type n) : 
    Boundary< Chain>( n), manifold_( n) {}
  Boundary_with_manifold( Chain& b, Chain &m) : 
    Boundary< Chain>( b), manifold_( m) {}
  
  // Access manifold
  const Chain& manifold() const     { return manifold_; }

  void clear()
  { 
    boundary_.clear();  
    manifold_.clear();
  }
  
  void print() const 
  { 
    boundary_.print(); 
    std::cout << "| ";
    manifold_.print();
  }
  void push_back_manifold( const std::pair< Pivot_handle, Field> &p)
  { 
    manifold_.push_back( p); 
  }
  
  // eliminate_pivot:  
  // boundary2 must have same pivot and coefficient equal to 1.
  // also does same operations to manifold chain
  template< class BoundaryWithManifold>
  void eliminate_pivot( const BoundaryWithManifold &boundary2,
			Compare_wrto_filtration handle_compare =
			Compare_wrto_filtration())
  {
    const typename Chain::Field& pivot_coeff = boundary_.coeff(pivot());
    boundary_.multiply_and_add( boundary2.boundary(), -pivot_coeff, 
				handle_compare);
    manifold_.multiply_and_add( boundary2.manifold(), -pivot_coeff,
				handle_compare);
    // delete pivot 
    // note:  this would've happened in multiply_and_add if we had
    // explicit pivots
    pop_pivot();
  }

  void normalize() 
  {
    if( empty()) return;
    const Field& pivot_coeff = boundary_.coeff(pivot());
    if( pivot_coeff != Chain::one) {
      Field inverse = Chain::one / pivot_coeff;
      boundary_.multiply( inverse);
      manifold_.multiply( inverse);
    }
  }
};

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_BOUNDARY_H
