// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_STANDARD_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_CRITERIA_H

#include <CGAL/license/Mesh_3.h>



#include <optional>
#include <boost/ptr_container/ptr_vector.hpp>

namespace CGAL {

namespace Mesh_3 {


/**
 * @class Abstract_criterion
 */
template <typename C3T3, typename Visitor_>
class Abstract_criterion
{
  typedef Abstract_criterion<C3T3, Visitor_> Self;

  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Geom_traits::FT FT;

public:
  typedef FT Quality;
  typedef std::optional<Quality>  Is_bad;
  typedef typename Visitor_::Handle Handle;

  /// Destructor
  virtual ~Abstract_criterion() {}

  void accept(Visitor_& v) const { do_accept(v); }
  Is_bad is_bad(const C3T3& c3t3, const Handle& h) const { return do_is_bad(c3t3, h); }
  Self* clone() const { return do_clone(); }

protected:
  virtual void do_accept(Visitor_& v) const = 0;
  virtual Is_bad do_is_bad(const C3T3& c3t3, const Handle& h) const = 0;
  virtual Self* do_clone() const = 0;

};  // end class Abstract_criterion


template <typename C3T3, typename Visitor_>
inline
Abstract_criterion<C3T3,Visitor_>*
new_clone(const Abstract_criterion<C3T3,Visitor_>& criterion)
{
  return criterion.clone();
}



template <typename C3T3, typename Handle_>
class Criterion_visitor
{
  typedef Criterion_visitor<C3T3, Handle_> Self;
public:
  typedef Handle_ Handle;

protected:
  typedef Abstract_criterion<C3T3, Self> Criterion;

public:
  typedef std::pair<int, typename Criterion::Quality> Quality;
  typedef std::optional<Quality>                    Is_bad;


  // Constructor
  Criterion_visitor(const C3T3& c3t3, const Handle_& h)
    : c3t3_(c3t3)
    , handle_(h)
    , is_bad_()
    , criterion_counter_(0) {}

  // Destructor
  ~Criterion_visitor() {}

  Is_bad is_bad() const
  {
    return is_bad_;
  }

  bool go_further() const
  {
    return !is_bad_;
  }

protected:
  void increment_counter()
  {
    ++criterion_counter_;
  }

  void set_is_bad(const Is_bad& is_bad)
  {
    is_bad_ = is_bad;
  }

  template<typename Derived>
  void do_visit(const Abstract_criterion<C3T3, Derived>& criterion)
  {
    typedef typename Abstract_criterion<C3T3, Derived>::Is_bad Is_bad;

    const Is_bad is_bad = criterion.is_bad(c3t3_, handle_);
    if ( is_bad )
      is_bad_ = std::make_pair(criterion_counter_, *is_bad);

    increment_counter();
  }

private:
  const C3T3& c3t3_;
  Handle_ handle_;
  Is_bad is_bad_;
  int criterion_counter_;
};  // end class Criterion_visitor




template <typename C3T3, typename Visitor_>
class Criteria
{
  typedef Criteria<C3T3, Visitor_> Self;

  typedef Abstract_criterion<C3T3, Visitor_> Criterion;
  typedef boost::ptr_vector<Criterion> Criterion_vector;
  typedef typename Visitor_::Quality Quality;
  typedef typename Visitor_::Is_bad  Is_bad;

public:
  /// Constructor
  Criteria() {}

  /// Copy constructor
  Criteria(const Self& rhs)
    : criterion_vector_(rhs.criterion_vector_.clone())    { }

  /// Destructor
  ~Criteria() {} // ptr_vector do the job of criterion deleting

  /// Add a criterion
  void add(Criterion* criterion)
  {
    criterion_vector_.push_back(criterion);
  }

  Is_bad operator()(const C3T3& c3t3,
                    const typename Visitor_::Handle& h) const
  {
    Visitor_ visitor(c3t3, h);

    typename Criterion_vector::const_iterator it = criterion_vector_.begin();
    for (  ; it != criterion_vector_.end() ; ++it )
    {
      it->accept(visitor);
      if ( ! visitor.go_further() )
        return visitor.is_bad();
    }

    return visitor.is_bad();
  }

private:
  Criterion_vector criterion_vector_;
};  // end class Criteria


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_STANDARD_CRITERIA_H
