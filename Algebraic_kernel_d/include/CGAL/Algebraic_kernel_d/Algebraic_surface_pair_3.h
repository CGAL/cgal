// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_PAIR_3_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_PAIR_3_H 1

/*!\file include/CGAL/Algebraic_kernel_d/Algebraic_surface_pair_3.h
 * \brief definition of \c Algebraic_surface_pair_3<>
 */

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arrangement_2l/z_stack_predeclarations.h>
#include <CGAL/Arrangement_2l/Surface_pair_3.h>

CGAL_BEGIN_NAMESPACE

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits, class Rep_ >
class Algebraic_surface_pair_3;

namespace CGALi {

template < class SurfaceZAtXyIsolatorTraits >
class Algebraic_surface_pair_3_rep : 
        public CGAL::CGALi::Surface_pair_3_rep< SurfaceZAtXyIsolatorTraits > {
    
public:
    //! this instance's template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! base class
    typedef CGAL::CGALi::Surface_pair_3_rep< Surface_z_at_xy_isolator_traits > 
    Base;

    //! the instance itself
    typedef Algebraic_surface_pair_3_rep< Surface_z_at_xy_isolator_traits >
    Self;

    //! Container for trivariate polynomials
    typedef std::vector< Polynomial_3 > Polynomial_3_container;

public:
    //!\name Constructors
    //!@{
    
    Algebraic_surface_pair_3_rep(
            const Surface_3& surface1, 
            const Surface_3& surface2) :
        Base(surface1, surface2) {
    }
    
    //!@}
    
public:

    // The Subresultant sequence of both surfaces
    mutable boost::optional<Polynomial_3_container> _m_subresultants;
    
    friend class Algebraic_surface_pair_3< Surface_z_at_xy_isolator_traits,
    Self >;
};

} // namespace CGALi

template < 
class SurfaceZAtXyIsolatorTraits, 
class Rep_ = 
CGAL::CGALi::Algebraic_surface_pair_3_rep < SurfaceZAtXyIsolatorTraits >
>
class Algebraic_surface_pair_3 : 
        public CGAL::Surface_pair_3<  SurfaceZAtXyIsolatorTraits, Rep_ > {
    
public:
    //! this instance's first template parameter
    typedef SurfaceZAtXyIsolatorTraits Surface_z_at_xy_isolator_traits;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;
    
    CGAL_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(
            Surface_z_at_xy_isolator_traits
    );
    
    //! type of Base
    typedef CGAL::Surface_pair_3< Surface_z_at_xy_isolator_traits, Rep > Base;

    //! this instance itself
    typedef Algebraic_surface_pair_3< Surface_z_at_xy_isolator_traits, Rep > 
    Self;

    //! type of restricted cad
    typedef typename Base::Restricted_cad_3 Restricted_cad_3;
    
    //!\name Constructors
    //!@{
    
    //! standard constructor
    Algebraic_surface_pair_3(
            const Surface_3& surface1, const Surface_3& surface2) :
        Base(surface1, surface2) {
    }
    
    //!@}

    //!\name Caching
    //!@{
private:
    
    typedef typename Surface_3::Surface_less_than Surface_less_than;
    
    //! type of Surface_pair
    typedef std::pair< Surface_3, Surface_3 > Surface_pair;
    
    //! type of Less of surface pair
    typedef CGAL::Pair_lexicographical_less_than< 
                            Surface_3, Surface_3, 
                            Surface_less_than, Surface_less_than > 
    Surface_pair_less;
    
    struct Pair_creator {
        Self operator()(Surface_pair pair) {
            return Self(pair.first, pair.second);
        }
    };

    //! pair canonicalizer
    struct Canonicalizer {
        
        //! first in pair in always less than second
        Surface_pair operator()(Surface_pair pair) {
            Surface_less_than less;
            if (!less(pair.first, pair.second)) {
                std::swap(pair.first, pair.second);
            }
            return pair;
        }
    };

public:
    //! type of surface pair cache
    typedef CGAL::Cache< Surface_pair, Self, 
    Pair_creator, Canonicalizer, 
    Surface_pair_less > Surface_pair_cache;
    
    //! instance of surface pair cache
    static
    Surface_pair_cache& surface_pair_cache() {
        static Surface_pair_cache _m_cache;
        return _m_cache;
    }
    
    //!@}
    
    //!\name (Sub)Resultants
    //!@{
private:
    
     //! Computes the subresultant sequence of the polynomials
    void _compute_subresultants() const {
        if (!this->ptr()->_m_subresultants) {
            typename Rep::Polynomial_3_container sres;
            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Polynomial_subresultants()(this->surface1().f(),
                                             this->surface2().f(),
                                             std::back_inserter(sres));
            this->ptr()->_m_subresultants = sres;
        }
    }

public:
    
    //! Compute the ith polynomial subresultant of the surfaces
    Polynomial_3 get_polynomial_subresultant(int i) {

      if (!this->ptr()->_m_subresultants ) {
          _compute_subresultants();
      }
      CGAL_assertion(i >= 0 
                     && i < static_cast<int>
                     (this->ptr()->_m_subresultants->size()));
      return this->ptr()->_m_subresultants.get()[i];
      
    }

    //! Compute the ith principal subresultant of the surfaces
    Polynomial_2 get_principal_subresultant(int i) {

      return get_polynomial_subresultant(i).lcoeff();
    
    }
    
    //!@}
};

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_PAIR_3_H
// EOF
