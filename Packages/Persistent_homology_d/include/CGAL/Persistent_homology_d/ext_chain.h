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

#ifndef CGAL_PERSISTENT_HOMOLOGY_D_EXT_CHAIN_H
#define CGAL_PERSISTENT_HOMOLOGY_D_EXT_CHAIN_H

#include <vector>
#include <functional>
#include <iostream>
#include <iterator>
#include <algorithm> // for set_symmetric_difference, etc.
#include <CGAL/Persistent_homology_d/chain.h>


namespace CGAL {
namespace Persistent_homology_d {

// Class:  Ext_Chain_Z2
// ********************
// Extended class for Chain_Z2 with functionality needed for
// persistence computation.  Also has static temporary space
// it allocates and uses.
template <class PivotHandle, 
          class CompareWrtoFiltration = std::less< PivotHandle> >
class Ext_Chain_Z2 {
  public:
  // types
  typedef PivotHandle                     Pivot_handle;
  typedef Ext_Chain_Z2<Pivot_handle>      Self;
  typedef Chain_Z2< PivotHandle, CompareWrtoFiltration> 
                                          Chain;
  typedef typename Chain::iterator        iterator;
  typedef typename Chain::const_iterator  const_iterator;
  typedef typename Chain::reference       reference;
  typedef typename Chain::const_reference const_reference;
  typedef typename Chain::size_type       size_type;
  typedef CompareWrtoFiltration           Compare_wrto_filtration;
  typedef int                             Field;

  // how we store 1 in the field for this chain (which is Z_2, so we use 1)
  static const Field one;

  protected:
  // members
  Chain             chain_;
  iterator          end_;

  // static temp
  static Chain      temp;  
  static size_type  temp_size;
  static iterator   temp_end;
 
  public:
  // constructors
  Ext_Chain_Z2() {}
  Ext_Chain_Z2( size_type n) : chain_( n) {
    // allocate and initialize chain
    chain_.reserve( n);
    typename Chain::value_type tmp;
    for( size_type i = 0; i < n; ++i) {
      chain_.push_back( tmp);
    }
    end_ = chain_.begin();
    // temp reallocation, if necessary
    if (temp_size < n) {
      temp.reserve( n);
      for( size_type i = temp_size; i < n; ++i) {
	temp.push_back( tmp);
      }
      temp_size = n;
      temp_end  = temp.begin();
    }
  }
  
  // extend vector operations out
  const Chain&    chain()     const { return chain_;                   }
  void            clear()           { end_ = chain_.begin();           }
  bool            empty()     const { return end_ == chain_.begin();   }
  size_type       size()      const { return end_ - chain_.begin();    }
  iterator        begin()           { return chain_.begin();           }
  iterator        end()             { return end_;                     }
  const_iterator  begin()     const { return chain_.begin();           }
  const_iterator  end()       const { return end_;                     }

  // coeff, handle accessors
  // no reference type for coefficient as you cannot set it!
  Pivot_handle handle        (iterator i)        const { return *i;   }
  const Field& coeff         (const_iterator i)  const { return one;  }
  const Pivot_handle handle  (const_iterator i)  const { return *i;   }

  // disregard the coefficient, have overwrite semantics
  void            push_back(const std::pair< Pivot_handle, int> &p)
  {
    *end_++ = p.first;
  }

  void            pop_back()        { --end_; }

  // multiply_and_add:  me += scalar * chain2
  template< class Chain2>
  void multiply_and_add( const Chain2 &chain2, 
			 int scalar = 1, 
			 CompareWrtoFiltration handle_compare =
			 CompareWrtoFiltration()) 
  {
    // overwrite semantics
    end_ = std::set_symmetric_difference( chain_.begin(), end_,
					  chain2.begin(), chain2.end(),
					  temp.begin(),
					  handle_compare);
    swap(chain_, temp);
  }

  // multiply:  me *= scalar
  void multiply( const Field& scalar) {
    // should never be called
  }  
  // print
  void print() const
  {
    for(const_iterator i = chain_.begin(); i != end_; i++) {
      handle(i)->print();
      std::cout << ' ';
    }
  }
};

// Definition of static variables
template <class PivotHandle, 
          class CompareWrtoFiltration>
typename Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::Chain
         Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::temp;
template <class PivotHandle, 
          class CompareWrtoFiltration>
typename Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::size_type
         Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::temp_size = 0;
template <class PivotHandle, 
          class CompareWrtoFiltration>
typename Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::iterator
         Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::temp_end;
template <class PivotHandle, 
          class CompareWrtoFiltration>
const typename Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::Field 
         Ext_Chain_Z2<PivotHandle, CompareWrtoFiltration>::one = 1;

// Class:  Ext_Chain_field
// ***********************
// Chains of PivotHandles with coefficients in Field.
// Coefficients are stored explicitly.

template <class PivotHandle,
          class Field_,
          class CompareWrtoFiltration = std::less< PivotHandle> >
class Ext_Chain_field {
  public:
  // types
  typedef Field_                                           Field;
  typedef Ext_Chain_field<PivotHandle, Field>              Self;
  typedef Chain_field< PivotHandle, Field_, CompareWrtoFiltration> 
                                                           Chain;
  typedef typename Chain::iterator                         iterator;
  typedef typename Chain::const_iterator                   const_iterator;
  typedef PivotHandle                                      Pivot_handle;
  typedef CompareWrtoFiltration                   Compare_wrto_filtration;
  typedef typename Chain::size_type                        size_type;
  typedef typename Chain::reference                        reference;
  typedef typename Chain::const_reference                  const_reference;
  // how we store 1 in the field for this chain
  static const Field one;
  static const Field zero;

  protected:
  // members
  Chain             chain_;
  iterator          end_;

  // static temp
  static Chain      temp;  
  static size_type  temp_size;
  static iterator   temp_end;
 
  public:
  // constructors
  Ext_Chain_field() {}
  Ext_Chain_field( size_type n) {
    // allocate and initialize chain
    chain_.reserve( n);
    typename Chain::value_type tmp;
    for( size_type i = 0; i < n; ++i) {
      chain_.push_back( tmp);
    }
    end_ = chain_.begin();
    // temp reallocation, if necessary
    if (temp_size < n) {
      temp.reserve( n);
      for( size_type i = temp_size; i < n; ++i) {
	temp.push_back( tmp);
      }
      temp_size = n;
      temp_end  = temp.begin();
    }
  }

  const Chain&    chain()     const { return chain_;                   }
  void            clear()           { end_ = chain_.begin();           }
  bool            empty()     const { return end_ == chain_.begin();   }
  size_type       size()      const { return end_ - chain_.begin();    }
  iterator        begin()           { return chain_.begin();           }
  iterator        end()             { return end_;                     }
  const_iterator  begin()     const { return chain_.begin();           }
  const_iterator  end()       const { return end_;                     }

  // coeff, handle accessors
  Field& coeff              (iterator i)       const { return i->second; }
  PivotHandle handle        (iterator i)       const { return i->first;  }
  const Field& coeff        (const_iterator i) const { return i->second; }
  const PivotHandle handle  (const_iterator i) const { return i->first;  }

  // disregard the coefficient, overwrite semantics
  void            push_back(const std::pair< Pivot_handle, Field> &p)
  {
    *end_++ = p;
  }

  void            pop_back()        { --end_; }

  // multiply_and_add:  me += scalar * chain2
  template< class Chain2>
  void multiply_and_add( const Chain2 &chain2, 
			 const Field& scalar = Field(1), 
			 CompareWrtoFiltration handle_compare =
			 CompareWrtoFiltration())
  {
    //std::cout << "scalar is " << scalar << std::endl;
    // "clear" temp
    temp_end = temp.begin();
    const_iterator i1 = chain_.begin();
    const_iterator i2 = chain2.begin();
    while ( i1 != end_ && i2 != chain2.end()) {
      if ( handle( i1) == handle( i2)) {
	Field c = coeff( i1) + scalar * coeff( i2);
	if ( c != zero) {
	  *temp_end++ = std::make_pair( handle( i1), c);
	}
	// either way increment both!
	i1++; i2++;
      } else if( handle_compare( handle( i1), handle( i2))) {
	*temp_end++ = *i1++;
      } else {
	*temp_end++ = std::make_pair( handle( i2), scalar * coeff( i2));
	i2++;
      }
    }
    // Copy rest of whatever remains:  i1 remains same,
    // but i2 is multiplied by scalar
    while(i1 != end_) {
      *temp_end++ = *i1++;
    }
    while(i2 != chain2.end()) {
      *temp_end++ = std::make_pair( handle( i2), scalar * coeff( i2));
      i2++;
    }
    // put the result back
    swap(chain_, temp);
    end_ = temp_end;
  }
  
  // multiply:  me *= scalar
  void multiply( const Field& scalar) 
  {
    for(iterator i = chain_.begin(); i != end_; i++) {
      coeff(i) *= scalar;
    }    
  }
  
  // print
  void print() const
  {
    for(const_iterator i = chain_.begin(); i != end_; i++) {
      std::cout << "(" << coeff(i) << ") ";
      handle(i)->print();
      std::cout << ' ';
    }
  }
};

// Definition of static variable
template< class PivotHandle,
          class Field_,
          class CompareWrtoFiltration
typename Ext_Chain_field< PivotHandle, Field_, CompareWrtoFiltration>::Chain
         Ext_Chain_field< PivotHandle, Field_, CompareWrtoFiltration>::temp;
template< class PivotHandle, 
	  class Field_,
          class CompareWrtoFiltration>
typename Ext_Chain_field< PivotHandle, Field_, 
			  CompareWrtoFiltration>::size_type
	 Ext_Chain_field< PivotHandle, Field_,
			  CompareWrtoFiltration>::temp_size = 0;
template< class PivotHandle, 
	  class Field_, 
          class CompareWrtoFiltration>
typename Ext_Chain_field< PivotHandle, Field_, 
			  CompareWrtoFiltration>::iterator
         Ext_Chain_field< PivotHandle, Field_, 
			  CompareWrtoFiltration>::temp_end;
template< class PivotHandle,
          class Field_,
          class CompareWrtoFiltration>
const Field_ Ext_Chain_field< PivotHandle, Field_, 
			      CompareWrtoFiltration>::zero = Field_(0);
template< class PivotHandle,
          class Field_,
          class CompareWrtoFiltration>
const Field_ Ext_Chain_field< PivotHandle, Field_, 
			      CompareWrtoFiltration>::one = Field_(1);

} // namespace Persistent_homology_d
} // namespace CGAL


#endif // CGAL_PERSISTENT_HOMOLOGY_D_EXT_CHAIN_H

