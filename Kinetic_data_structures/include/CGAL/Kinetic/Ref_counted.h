// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_REF_COUNTED_H
#define CGAL_REF_COUNTED_H
#include <CGAL/Kinetic/basic.h>
#include <boost/intrusive_ptr.hpp>
#include <boost/utility.hpp>

//#define NEW_REF_COUNTED

namespace CGAL { namespace Kinetic { namespace internal {

class Ref_counted_base;

void intrusive_ptr_add_ref(const Ref_counted_base *t);

void intrusive_ptr_release(const Ref_counted_base *t);


class Ref_counted_base: boost::noncopyable
{
  typedef Ref_counted_base This;

public:
  Ref_counted_base() : reference_count_(0) {}

  void write(std::ostream &out) const
  {
    out << "(" << reference_count_ << ")";
  }
  bool is_referenced() const
  {
    return reference_count_ != 0;
  }

  virtual ~Ref_counted_base() CGAL_NOEXCEPT(CGAL_NO_ASSERTIONS_BOOL)
  {
    CGAL_destructor_assertion(reference_count_==0);
  }

  friend void intrusive_ptr_release(const This *t);
  friend void intrusive_ptr_add_ref(const This *t);

  unsigned int reference_count() const {return reference_count_;}

  void  new_ref() const { ++reference_count_; }
  void  delete_ref() const
  {
    CGAL_precondition(reference_count_!=0);
    --reference_count_;
  }

  mutable unsigned int reference_count_;
};

inline void intrusive_ptr_add_ref(const Ref_counted_base *t)
{
  t->new_ref();
}


inline void intrusive_ptr_release(const Ref_counted_base *t)
{
  t->delete_ref();
  if (t->reference_count() == 0) {
    delete t;
  }
}



struct Non_ref_counted_base{};

inline void intrusive_ptr_add_ref(const Non_ref_counted_base *)
{
}


inline void intrusive_ptr_release(const Non_ref_counted_base *)
{
}



} } } //namespace CGAL::Kinetic::internal

namespace CGAL { namespace Kinetic {

template <class T>
class Ref_counted: public internal::Ref_counted_base

{
  typedef internal::Ref_counted_base P;
 public:
  typedef T This;
  typedef typename boost::intrusive_ptr<T> Handle;

  typedef typename boost::intrusive_ptr<const T> Const_handle;
};

template <class T>
struct Non_ref_counted: public internal::Non_ref_counted_base
{
  typedef T This;
  typedef typename boost::intrusive_ptr<T> Handle;
  typedef typename boost::intrusive_ptr<const T> Const_handle;
};



} } //namespace CGAL::Kinetic
#endif
