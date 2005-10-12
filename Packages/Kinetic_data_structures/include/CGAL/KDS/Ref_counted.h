// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_REF_COUNTED_H
#define CGAL_REF_COUNTED_H
#include <CGAL/KDS/basic.h>
#include <boost/intrusive_ptr.hpp>

//#define NEW_REF_COUNTED

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE;

class Ref_counted_base;

void intrusive_ptr_add_ref(const Ref_counted_base *t);

void intrusive_ptr_release(const Ref_counted_base *t);

class Ref_counted_base {
  typedef Ref_counted_base This;
  Ref_counted_base(const Ref_counted_base&) : reference_count_(0) {
    std::cerr << "Copy constructor called. Why?" << std::endl;
  }
  Ref_counted_base operator=(const Ref_counted_base&) const {
    // preserve current reference count
    std::cerr << "Assignment called. Why?" << std::endl;
    return *this;
  }
public:

  //! Initialize the count to 0.
  Ref_counted_base() : reference_count_(0) {}

  void write(std::ostream &out) const {
    out << "(" << reference_count_ << ")";
  }

  //! Use this to verify that an object is allocated on the heap.
  bool is_referenced() const {
    return reference_count_ != 0;
  }

  virtual ~Ref_counted_base(){
    CGAL_assertion(reference_count_==0);
  }

  //protected:

  friend void intrusive_ptr_release(const This *t);
  friend void intrusive_ptr_add_ref(const This *t);

  unsigned int reference_count() const {return reference_count_;}

  void  new_ref() const { ++reference_count_; }
  void  delete_ref() const { 
    CGAL_precondition(reference_count_!=0);
    --reference_count_;
  }

  mutable unsigned int reference_count_;
};


inline void intrusive_ptr_add_ref(const Ref_counted_base *t) {
  t->new_ref();
}

inline void intrusive_ptr_release(const Ref_counted_base *t) {
  t->delete_ref();
  if (t->reference_count() == 0){
    delete t;
  }
}

CGAL_KDS_END_INTERNAL_NAMESPACE;

CGAL_KDS_BEGIN_NAMESPACE;
//#ifndef NEW_REF_COUNTED
template <class T>
class Ref_counted;

/*template <class T>
void intrusive_ptr_add_ref(const T *t) {
  //t->new_ref();
  internal::intrusive_ptr_add_ref(t);
}
template <class T>
void intrusive_ptr_release(const T *t) {
  internal::intrusive_ptr_release(t);
}*/


/*inline void intrusive_ptr_add_ref(const internal::Ref_counted_base *t) {
  assert(0);
  t->new_ref();
}

inline void intrusive_ptr_release(const internal::Ref_counted_base *t) {
  assert(0);
  t->delete_ref();
  if (t->reference_count() == 0){
    delete t;
  }
  }*/
//#endif

//! The base class for ref counted objects
/*!
  The class T is there to avoid objects having a common base
  class. The best value for T is the type of the thing actually
  being pointed to. Then the conventional Pointer type is
  automatically defined.

  As a convention, any ref counted object should define a type
  Pointer which is the type of the reference counted pointer.
*/
template <class T=int>
class Ref_counted
//#ifdef NEW_REF_COUNTED
  : public internal::Ref_counted_base 
    //#endif

{
  typedef internal::Ref_counted_base P;
  typedef Ref_counted<T> This;
  //! This is necessary
  /*!  I need this constructor since the reference count needs to be
    reset on a copy.
   */
  Ref_counted(const Ref_counted&) {
    this_should_not_compile(T());
    assert(0);
  }

 
  
  This operator=(const This &) {
    this_should_not_compile(T());
    assert(0);
    return *this;
  }
public:

  //! Initialize the count to 0.
  //#ifdef NEW_REF_COUNTED
  Ref_counted() {} /*reference_count_(0) {}*/
  //#else
  //Ref_counted(): reference_count_(0) {}
  //#endif
  
  //! The pointer for an object of type T.
  /*!
    If T is the type inheriting from Ref_counted_base,
    then all is good. If not, the type should define its
    own Pointer type.
  */
  typedef typename boost::intrusive_ptr<T> Pointer;
  //! Constant pointer
  typedef typename boost::intrusive_ptr<const T> Const_pointer;
#ifndef NEW_REF_COUNTED
  /*void write(std::ostream &out) const {
    out << "(" << reference_count_ << ")";
    }*/

  //! Use this to verify that an object is allocated on the heap.
  /*bool is_referenced() const {
    return P::reference_count_ != 0;
    }*/

protected:

  //friend void intrusive_ptr_release<T>(const T *t);
  //friend void intrusive_ptr_add_ref<T>(const T *t);

  //unsigned int reference_count() const {return P::reference_count_;}

  /*void  new_ref() const { ++P::reference_count_; }
  void  delete_ref() const { 
    CGAL_precondition(P::reference_count_!=0);
    --P::reference_count_;
    }*/

private:
  //mutable unsigned int reference_count_;
#endif
};



CGAL_KDS_END_NAMESPACE;
#endif
