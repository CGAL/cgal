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

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_PERSISTENCE_H
#define CGAL_PERSISTENT_HOMOLOGY_D_PERSISTENCE_H

#include <CGAL/basic.h>
#include <vector>
#include <iterator>
#include <iostream>

#include <CGAL/Unique_hash_map.h>

#include <CGAL/Persistent_homology_d/chain.h>
#include <CGAL/Persistent_homology_d/ext_chain.h>
#include <CGAL/Persistent_homology_d/Boundary.h>
#include <CGAL/Persistent_homology_d/Pivot.h>

namespace CGAL {
namespace Persistent_homology_d {


// **********************************************************************
// Coefficient Policy                                                   *
// **********************************************************************
// Note:  policy allows chains to be defined

class Z2_coefficients {
public:
  typedef int Field;
  template <class Pivot_handle, class CompareWrtoFiltration>
    struct Bind {
      typedef ::CGAL::Persistent_homology_d::
              Chain_Z2< Pivot_handle, CompareWrtoFiltration>     Chain;
      typedef ::CGAL::Persistent_homology_d::
              Ext_Chain_Z2< Pivot_handle, CompareWrtoFiltration> Ext_Chain;
    };
};

template < class Field_>
class Field_coefficients {
public:
  typedef Field_ Field;
  template <class Pivot_handle, class CompareWrtoFiltration>
  struct Bind {
    typedef ::CGAL::Persistent_homology_d::Chain_field< Pivot_handle, Field, 
			       CompareWrtoFiltration> Chain;
    typedef ::CGAL::Persistent_homology_d::Ext_Chain_field<Pivot_handle, Field,
			       CompareWrtoFiltration> Ext_Chain;
  };
};

// **********************************************************************
// Storage Policies                                                     *
// **********************************************************************

// Our storage, provided by us
class Our_storage {
  public:
  template < class SimplexHandle, class OutputPolicy, 
	     class CoefficientPolicy>
  class Pivot_space {
  public:
    typedef SimplexHandle                                Simplex_handle;
    typedef OutputPolicy                                 Output_policy;
    typedef ::CGAL::Persistent_homology_d::Pivot < SimplexHandle, OutputPolicy,
		      CoefficientPolicy>                 Pivot;
    typedef typename Pivot::Pivot_handle                 Pivot_handle;
    typedef typename Pivot::CompareWrtoFiltration        CompareWrtoFiltration;
    typedef std::size_t                                  size_type;
    typedef typename Pivot::Boundary                     Boundary;
    typedef typename Pivot::Ext_Boundary                 Ext_Boundary;
    // defined positive and negative infinities
    Pivot_handle                                         pos_inf;
    Pivot_handle                                         neg_inf;
  private:
    // private members
    std::vector< Pivot> pivots;
    CGAL::Unique_hash_map< Simplex_handle, Pivot_handle> hash;
    typename std::vector< Pivot>::iterator               next;
  public:    
    Pivot_space( size_type n) {
      pivots.reserve( n+1);
      pivots.insert( pivots.begin(), n+1, Pivot());
      pos_inf = pivots.end();
      neg_inf = pivots.begin();
      next = neg_inf + 1;
      hash = CGAL::Unique_hash_map< Simplex_handle, Pivot_handle>(neg_inf);
      pivots[0].set_base(next);  // for printing capability
    }

    Pivot_handle pivot_handle( Simplex_handle s) 
    { 
      Pivot_handle& pivot = hash[s];
      if ( pivot == neg_inf) {
 	pivot = next++;
 	pivot->simplex() = s;
      }
      return pivot;
    }

    Pivot_handle pivot_handle_from_index( size_type i)
    {
      return neg_inf + i + 1;
    }

    size_type pivot_index_from_handle( Pivot_handle p)
    {
      return p - neg_inf - 1;
    }

    Boundary& boundary( Pivot_handle p)        { return p->boundary();       }
    Pivot_handle& paired_with( Pivot_handle p) { return p->paired_with();    }
    bool is_paired( Pivot_handle p) { 
      return 
	p->paired_with() != pos_inf &&
	p->paired_with() != neg_inf;
    }
    bool is_positive( Pivot_handle p)          { return p->is_positive();    }
    void set_positive( Pivot_handle p)         { p->paired_with() = pos_inf; }
    void set_negative( Pivot_handle p)         { p->paired_with() = neg_inf; }
  };
};

// Example user storage, with mostly empty functions
class User_storage {
  template < class SimplexHandle, class Boundary_>
  class Pivot_space {
  public:
    typedef SimplexHandle Pivot_handle;
    typedef SimplexHandle Simplex_handle;
    typedef Boundary_     Boundary;
    typedef typename std::iterator_traits< Simplex_handle>
      ::size_type size_type;

    Pivot_space( size_type n)                           {}
    Pivot_handle pivot_handle( Simplex_handle s)        { return s; }
    Pivot_handle pivot_handle_from_index( size_type i)  {}
    Boundary& boundary( Pivot_handle p)                 {}
    Pivot_handle& paired_with( Pivot_handle p)          {}
    bool is_positive ( Pivot_handle p)                  {}
    typedef std::less< Pivot_handle>                    CompareWrtoFiltration;
  };
};

// **********************************************************************
// Output Policies                                                      *
// **********************************************************************

class Pairs_only {
public:
  template < class Chain>
  struct Bind {
    typedef ::CGAL::Persistent_homology_d::Boundary< Chain> Boundary;
  };
  template< class PivotHandle, class Boundary, class Field>
  void Initialize_manifold( PivotHandle p, Boundary& b) {}
};

class Pairs_and_boundary {
public:
  template < class Chain>
  struct Bind {
    typedef ::CGAL::Persistent_homology_d::Boundary< Chain> Boundary;
  };
  template< class PivotHandle, class Boundary, class Field>
  void Initialize_manifold( PivotHandle p, Boundary& b) {}
};

class Pairs_and_manifold {
public:
  template <class Chain>
  struct Bind {
    typedef ::CGAL::Persistent_homology_d::Boundary_with_manifold <Chain> 
        Boundary;
  };
  template< class PivotHandle, class Boundary, class Field>
  void Initialize_manifold( PivotHandle p, Boundary& b)
  {
    typedef PivotHandle Pivot_handle;
    b.push_back_manifold( std::pair< Pivot_handle, Field>
			  (p, b.boundary().one));
  }
};

// **********************************************************************
// Boundary Operator                                                    *
// **********************************************************************

class Boundary_handle_tag {};
class Boundary_index_tag  {};

// Example Boundary Operator
// MUST return a pair type, the second type being the coefficient type.
// class BoundaryOperator {
// public:
//   typedef Boundary_handle_tag rep_tag;
//   typedef int (*Simplex_handle)[4];
//   typedef std::pair< int, int> value_type;
//   template <class OutputIterator>
//   OutputIterator operator()(SimplexHandle, OutputIterator out) {}
//};

// Boundary_output_transformer
template <class BoundaryOperator, class PivotSpace, class OutputPolicy>
class Boundary_output_transformer {
public:
  // types
  typedef Boundary_output_transformer<BoundaryOperator, 
    PivotSpace, OutputPolicy>                   Self;
  typedef std::output_iterator_tag              iterator_category;
  typedef void                                  difference_type;
  typedef void                                  pointer;
  typedef void                                  reference;
  typedef typename PivotSpace::Pivot_handle     Pivot_handle;
  typedef typename PivotSpace::Ext_Boundary     Ext_Boundary;
  typedef typename BoundaryOperator::rep_tag    rep_type;
  typedef typename BoundaryOperator::value_type value_type;
protected:
  Ext_Boundary*                    boundary;
  PivotSpace*                      pivot_space;
public:
  void action( value_type v, Boundary_handle_tag, Pairs_only)
  {
    Pivot_handle handle = pivot_space->pivot_handle( v.first);
    if(pivot_space->is_positive( handle)) {
      boundary->push_back( std::make_pair( handle, v.second));
    }
  }

  void action( value_type v, Boundary_index_tag, Pairs_only)
  {
    Pivot_handle handle = pivot_space->pivot_handle_from_index( v.first);
    if(pivot_space->is_positive( handle)) {
      boundary->push_back( std::make_pair( handle, v.second));
    }
  }

  template< class OutputPolicy_>
  void action( value_type v, Boundary_handle_tag, OutputPolicy_)
  {
    Pivot_handle handle = pivot_space->pivot_handle( v.first);
    boundary->push_back( std::make_pair( handle, v.second));
  }

  template< class OutputPolicy_>
  void action( value_type v, Boundary_index_tag, OutputPolicy_)
  {
    Pivot_handle handle = pivot_space->pivot_handle_from_index( v.first);
    boundary->push_back( std::make_pair( handle, v.second));
  }
  
  explicit Boundary_output_transformer(Ext_Boundary *b, PivotSpace *p)
    : boundary( b), pivot_space( p) {}
  Self& operator=( value_type p) { 
    action(p, rep_type(), OutputPolicy());
    return *this;
  }
  Self& operator*()     { return *this; }
  Self& operator++()    { return *this; }
  Self& operator++(int) { return *this; }
};


// **********************************************************************
// Persistence Function                                                 *
// **********************************************************************

template <class ForwardIterator,
          class BoundaryOperator,
	  class CoefficientPolicy, 
          class StoragePolicy,
	  class OutputPolicy>
typename StoragePolicy::template 
  Pivot_space< ForwardIterator, OutputPolicy, CoefficientPolicy>*
persistence( ForwardIterator   first, 
             ForwardIterator   last,
             CoefficientPolicy coefficient_policy,
             StoragePolicy     storage_policy,
             OutputPolicy      output_policy,
             BoundaryOperator  boundary_operator)
{
  // types
  typedef ForwardIterator                      Simplex_handle;
  typedef std::size_t                          size_type;
  typedef typename StoragePolicy::template
    Pivot_space<Simplex_handle, OutputPolicy, CoefficientPolicy>
                                               Pivot_space;
  typedef typename Pivot_space::Pivot_handle   Pivot_handle;
  typedef typename Pivot_space::Boundary       Boundary;
  typedef typename Pivot_space::Ext_Boundary   Ext_Boundary;
  typedef typename Boundary::iterator          Output_iterator;
  typedef Boundary_output_transformer 
    < BoundaryOperator, Pivot_space, OutputPolicy> Transformer;
  typedef typename CoefficientPolicy::Field    Field;
 
  // allocate space
  size_type n = std::distance( first, last); 
  Pivot_space *pivots = new Pivot_space( n);  
  Ext_Boundary b( n);
  Transformer transformer( &b, pivots);

  // process each boundary in turn 
  Pivot_handle p = pivots->pivot_handle_from_index( 0);
  for( Simplex_handle simplex = first; simplex != last; ++simplex, ++p) {
    // get simplex into our world:
    // 1.  hash table way:
    // const Pivot_handle p = pivots->pivot_handle(simplex);

    // 2.  non-hash table way.  we need to set simplex
    // we already have the right pivot_handle from for loop
    p->simplex() = simplex;
    
    // get boundary
    boundary_operator( simplex, transformer);
    
    // set up manifold, if needed
    output_policy.template 
      Initialize_manifold< Pivot_handle, Ext_Boundary, Field>( p, b);
    
    // process boundary
    process_boundary( pivots, p, b);
    // get rid of it
    b.clear();    
  }
  
  // output
  return pivots;
}

// process_boundary
// note on testing positive's boundary:
// we test if positive is already paired with somebody
// while in C, the test is something like 
// (positive's boundary == NULL), running empty() on it does not work
// as the it may still be empty for Pairs_only even though
// positive _is_ paired.  hence the current test
template < class PivotSpace>
void process_boundary( PivotSpace* pivots, 
		       typename PivotSpace::Pivot_handle simplex,
		       typename PivotSpace::Ext_Boundary &boundary)
{
  typedef PivotSpace                         Pivot_space;
  typedef typename Pivot_space::Boundary     Boundary;
  typedef typename Pivot_space::Pivot_handle Pivot_handle;

  while( !boundary.empty()) {
    //boundary.boundary().print(); std::cout << std::endl;
    const Pivot_handle positive = 
      boundary.boundary().handle( boundary.pivot());
    if( !pivots->is_paired( positive)) {
      // normalize _first_ and remove pivot (coeff 1 is implied)  
      boundary.normalize();
      boundary.pop_pivot();
      //std::cout << "storing ";
      //boundary.boundary().print(); std::cout << std::endl;
      pivots->boundary( positive) = boundary;
      pivots->paired_with( positive) = simplex;
      pivots->paired_with( simplex)  = positive;
      return;
    }
    // Collision -- eliminate pivot!
    const Boundary& stored_boundary = pivots->boundary( positive);
    //std::cout << "stored is ";
    //stored_boundary.boundary().print(); std::cout << std::endl;
    boundary.eliminate_pivot( stored_boundary);
  }
  // empty boundary means positive simplex
  pivots->set_positive( simplex);
}

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_PERSISTENCE_H
