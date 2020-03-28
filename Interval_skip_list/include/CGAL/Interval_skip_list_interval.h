// Copyright (c) 2003 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERVAL_SKIP_LIST_INTERVAL_H
#define CGAL_INTERVAL_SKIP_LIST_INTERVAL_H

#include <CGAL/license/Interval_skip_list.h>


#include <CGAL/basic.h>
#include <cstdlib>
#include <iostream>


namespace CGAL {

  template <class Value_>
  class Interval_skip_list_interval
  {
  public:
    typedef Value_ Value;

  private:
    bool lbound_;
    bool rbound_;
    Value inf_;
    Value sup_;
  public:

    Interval_skip_list_interval(){}
    Interval_skip_list_interval(const Value& inf_,
                                const Value& sup_,
                                bool lb = true,
                                bool rb = true);

    const Value& inf() const {return inf_;}

    const Value& sup() const {return sup_;}

    bool inf_closed() const {return lbound_;}

    bool sup_closed() const {return rbound_;}

    bool contains(const Value& V) const;

    // true iff this contains (l,r)
    bool contains_interval(const Value& l, const Value& r) const;

    bool operator==(const Interval_skip_list_interval& I) const
    {
      return ( (inf() == I.inf()) && (sup() == I.sup()) &&
               (inf_closed() == I.inf_closed()) && (sup_closed() == I.sup_closed()) );
    }

    bool operator!=(const Interval_skip_list_interval& I) const
    {
      return ! (*this == I);
    }
  };



  template <class V>
  std::ostream& operator<<(std::ostream& os,
                           const Interval_skip_list_interval<V>& i)
  {
    os << (i.inf_closed()?"[":"(") << i.inf() << ", " << i.sup() << (i.sup_closed()?"]":")");
    return os;
  }


  template <class V>
  Interval_skip_list_interval<V>::Interval_skip_list_interval(
                                                              const Value& i,
                                                              const Value& s,
                                                              bool lb, bool rb)
    : lbound_(lb), rbound_(rb), inf_(i), sup_(s)
  {
    CGAL_precondition( !(inf_ > sup_) );
  }


  template <class V>
  bool
  Interval_skip_list_interval<V>::contains_interval(const Value& i,
                                                    const Value& s) const
    // true iff this contains (l,r)
  {
    return( (inf() <= i) && (sup() >= s) );
  }


  template <class V>
  bool
  Interval_skip_list_interval<V>::contains(const Value& v) const
  {
    // return true if this contains V, false otherwise
    if((v > inf()) && (v < sup()))
      return true;
    else if ((v == inf()) && inf_closed())
      return true;
    else if ((v == sup()) && sup_closed())
      return true;
    else
      return false;
  }




} // namespace CGAL

#endif // CGAL_INTERVAL_SKIP_LIST_INTERVAL_H
