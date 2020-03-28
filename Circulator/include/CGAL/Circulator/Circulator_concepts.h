// Copyright (c) 2012
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philipp Möller

#ifndef CGAL_CIRCULATOR_CONCEPTS_H
#define CGAL_CIRCULATOR_CONCEPTS_H

#include <boost/concept/assert.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/detail/concept_def.hpp>

#include <iterator>

#include <CGAL/circulator_bases.h>

namespace CGAL { namespace Concepts {

template <typename T>
void deref(T& t) {
  t.operator->();
}

template <typename T>
void deref(T*) {
  // valid anyway
}

// N.B.: there is no such concept as Circulator as it is immaterial
template <typename C>
struct ForwardCirculator
  : boost::Assignable<C>, boost::DefaultConstructible<C>
  , boost::CopyConstructible<C>
{
  // for some odd reason circulators have no associated traits
  typedef typename C::value_type value_type;
  typedef typename C::reference reference;
  typedef typename C::pointer pointer;
  typedef typename C::difference_type difference_type;
  // odd name but useful for compatibility
  typedef typename C::iterator_category iterator_category;

  // this requirement ought to be a mistake because it breaks
  // compatibility with iterators
  // typedef typename C::size_type size_type;

  BOOST_CONCEPT_USAGE(ForwardCirculator)
  {
    // CGAL_UNUSED is needed with g++ 4.8, otherwise it warns:
    // "typedef ‘boost_concept_check68’ locally defined but not used
    // [-Wunused-local-typedefs]"
    BOOST_CONCEPT_ASSERT((boost::SignedInteger<difference_type>)) CGAL_UNUSED;
    BOOST_CONCEPT_ASSERT((boost::Convertible<iterator_category, CGAL::Forward_circulator_tag>)) CGAL_UNUSED;

    boost::require_boolean_expr(a == nullptr);
    boost::require_boolean_expr(a != nullptr);
    ++a;
    a++;
    (void)*a; // suppress unused warning, don't check the return type
              // to allow for proxies

    deref(a);
  }
private:
  C a;
};

template<typename C>
struct BidirectionalCirculator
  : ForwardCirculator<C>
{
  BOOST_CONCEPT_USAGE(BidirectionalCirculator)
  {
    BOOST_CONCEPT_ASSERT((boost::Convertible<typename C::iterator_category, CGAL::Bidirectional_circulator_tag>)) CGAL_UNUSED;
    --a;
    a--;
  }
private:
  C a;
};

template<typename C>
struct RandomAccessCirculator
  : BidirectionalCirculator<C>
{
  BOOST_CONCEPT_USAGE(RandomAccessCirculator)
  {
    BOOST_CONCEPT_ASSERT((boost::Convertible<typename C::iterator_category, CGAL::Random_access_circulator_tag>)) CGAL_UNUSED;
    c += n; // addition
    c = c + n; c = n + c;
    c -= n; // subtraction
    c = c - n;
    n = c - b; // difference
    (void)c[n]; // operator[]

    c.min_circulator(); // minimum
  }
private:
  C c, b;
  C i;
  typename C::difference_type n;
};

} // Concept
} // CGAL



#include <boost/concept/detail/concept_undef.hpp>


#endif /* CGAL_CIRCULATOR_CONCEPTS_H */
