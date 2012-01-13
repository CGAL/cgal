// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_STANDARD_CRITERIA_H
#define CGAL_MESH_3_MESH_STANDARD_CRITERIA_H


#include <boost/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace CGAL {

namespace Mesh_3 {


/**
 * @class Abstract_criterion
 */
template <typename Tr, typename Visitor_>
class Abstract_criterion
{
  typedef typename Tr::Geom_traits::FT FT;
  typedef Abstract_criterion<Tr,Visitor_> Self;

public:
  typedef FT Quality;
  typedef boost::optional<Quality> Badness;
  typedef typename Visitor_::Handle Handle;

  /// Destructor
  virtual ~Abstract_criterion() {}

  void accept(Visitor_& v) const { do_accept(v); }
  Badness is_bad(const Handle& h) const { return do_is_bad(h); }
  Self* clone() const { return do_clone(); }

protected:
  virtual void do_accept(Visitor_& v) const = 0;
  virtual Badness do_is_bad(const Handle& h) const = 0;
  virtual Self* do_clone() const = 0;

};  // end class Abstract_criterion


template <typename Tr, typename Visitor_>
inline
Abstract_criterion<Tr,Visitor_>*
new_clone(const Abstract_criterion<Tr,Visitor_>& criterion)
{
  return criterion.clone();
}



template <typename Tr, typename Handle_>
class Criterion_visitor
{
  typedef Criterion_visitor<Tr, Handle_> Self;
public:
  typedef Handle_ Handle;

protected:
  typedef Abstract_criterion<Tr, Self> Criterion;

public:
  typedef std::pair<int, typename Criterion::Quality> Quality;
  typedef boost::optional<Quality> Badness;


  // Constructor
  Criterion_visitor(const Handle_& h)
    : handle_(h)
    , badness_()
    , criterion_counter_(0) {}

  // Destructor
  ~Criterion_visitor() {}

  Badness badness() const
  {
    return badness_;
  }

  bool go_further() const
  {
    return !badness_;
  }

protected:
  void increment_counter()
  {
    ++criterion_counter_;
  }

  void set_badness(const Badness& badness)
  {
    badness_ = badness;
  }

  template<typename Derived>
  void do_visit(const Abstract_criterion<Tr, Derived>& criterion)
  {
    typedef typename Abstract_criterion<Tr, Derived>::Badness Badness;

    const Badness badness = criterion.is_bad(handle_);
    if ( badness )
      badness_ = std::make_pair(criterion_counter_, *badness);

    increment_counter();
  }

private:
  Handle_ handle_;
  Badness badness_;
  int criterion_counter_;

};  // end class Criterion_visitor




template <typename Tr, typename Visitor_>
class Criteria
{
  typedef Criteria<Tr, Visitor_> Self;

  typedef Abstract_criterion<Tr,Visitor_> Criterion;
  typedef boost::ptr_vector<Criterion> Criterion_vector;
  typedef typename Visitor_::Quality Quality;
  typedef typename Visitor_::Badness Badness;

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

  Badness operator()(const typename Visitor_::Handle& h) const
  {
    Visitor_ visitor(h);

    typename Criterion_vector::const_iterator it = criterion_vector_.begin();
    for (  ; it != criterion_vector_.end() ; ++it )
    {
      it->accept(visitor);
      if ( ! visitor.go_further() ) return visitor.badness();
    }

    return visitor.badness();
  }

private:
  Criterion_vector criterion_vector_;
};  // end class Criteria


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_STANDARD_CRITERIA_H
