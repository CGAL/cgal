// Copyright (c) 2001-2004  ENS of Paris (France).
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
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_BITANGENT_2_H
#define CGAL_VISIBILITY_COMPLEX_2_BITANGENT_2_H

#include <string>
#include <iostream>

CGAL_BEGIN_NAMESPACE

namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------
// -------------------- Base Bitangent class -----------------------------------

struct Bitangent_type_wrapper {
  enum Type  {
    LL,LR,RL,RR
  };
  static Type reverse(Type t) {
    switch (t) {
    case Bitangent_type_wrapper::LL: return Bitangent_type_wrapper::RR;
    case Bitangent_type_wrapper::LR: return Bitangent_type_wrapper::LR;
    case Bitangent_type_wrapper::RL: return Bitangent_type_wrapper::RL;
    case Bitangent_type_wrapper::RR: return Bitangent_type_wrapper::LL;
    default: CGAL_assertion(false); return Bitangent_type_wrapper::LL;
    }
  }
  class Type_util {
  public:
    Type operator()(bool b1,bool b2) const {
      if ( b1 &&  b2) return LL;
      if (!b1 &&  b2) return RL;
      if (!b1 && !b2) return RR;
      return LR;
    }
  };
};

typedef Bitangent_type_wrapper::Type Bitangent_type;

std::ostream &
operator<<(std::ostream &os, Bitangent_type t) {
  switch (t) {
  case Bitangent_type_wrapper::LL : return os<<"LL"; 
  case Bitangent_type_wrapper::LR : return os<<"LR";
  case Bitangent_type_wrapper::RL : return os<<"RL";
  default : return os<<"RR";
  }
}

std::istream &
operator>>(std::istream &is, Bitangent_type& t) {
  std::string s;
  is>>s;
  if (s=="LL") {
    t=Bitangent_type_wrapper::LL; return is;
  }
  if (s=="LR") {
    t=Bitangent_type_wrapper::LR; return is; 
  }
  if (s=="RL") {
    t=Bitangent_type_wrapper::RL; return is;
  }
  if (s=="RR") {
    t=Bitangent_type_wrapper::RR; return is; 
  }
  if (s=="") return is;
  CGAL_assertion(false);
  return is;
}

class Constraint_input :public Bitangent_type_wrapper {
public:
  Constraint_input() :t_(LL), source_(0), target_(0) {};
  Constraint_input(Type t,size_t source,size_t target):
    t_(t),source_(source),target_(target) {}
  Type type() const {return t_;}
  size_t source() const {return source_;}
  size_t target() const {return target_;}
private:
  Type t_;
  size_t source_,target_;
};

std::ostream &
operator<<(std::ostream &os, const Constraint_input& c) {
  return os<<c.type()<<" "<<c.source()<<" "<<c.target();
}

std::istream &
operator>>(std::istream &is, Constraint_input&c) {
  Bitangent_type t;
  size_t source,target;
  is>>t>>source>>target;
  if (is) c=Constraint_input(t,source,target);
  return is;
}



template< class D_ >
class Bitangent_base :public Bitangent_type_wrapper
{
public:
    // -------------------------------------------------------------------------
    typedef D_                                        Disk;
    typedef const Disk*                           Disk_handle;
    // -------------------------------------------------------------------------
    typedef Bitangent_type_wrapper::Type Type;
    using Bitangent_type_wrapper::LL;
    using Bitangent_type_wrapper::LR;
    using Bitangent_type_wrapper::RL;
    using Bitangent_type_wrapper::RR;
public:
    // Constructeurs -----------------------------------------------------------
    Bitangent_base() 
	: type_(LL) , source_object_(0)     , target_object_(0) { }
    Bitangent_base(Type t, Disk_handle o1, Disk_handle o2)  
	: type_(t)  , source_object_(o1) , target_object_(o2)   { }
    Bitangent_base(const Bitangent_base&sibling,bool reverse,Type t) {
      if (reverse) {
        type_=Bitangent_type_wrapper::reverse(t);
        source_object_=sibling.target_object();
        target_object_=sibling.source_object();
      } else {
        target_object_=sibling.target_object();
        source_object_=sibling.source_object();
        type_=t;        
      }
    }
    // ----------------------- Operators ---------------------------------------
    bool operator==(const Bitangent_base& b) const{
	return (b.type() == type() && source_object() == b.source_object() && 
				      target_object() == b.target_object());
    }
    bool operator!=(const Bitangent_base& b) const{ return !(*this == b); }
    // -------- return the opposite oriented bitangent -------------------------
    Type type()                    const { return type_;          }
    // ---- accesing the two objects defining the bitangent --------------------
    Disk_handle source_object() const { return source_object_; }
    Disk_handle target_object() const { return target_object_; }
    // ---------- informations on the type of the bitangent -------------------- 
    bool is_left_right()   const { return (type() == LR);                 }
    bool is_left_left()    const { return (type() == LL);                 }
    bool is_right_right()  const { return (type() == RR);                 }
    bool is_right_left()   const { return (type() == RL);                 }
    bool is_left_xx()      const { return (type() == LL || type() == LR); }
    bool is_xx_left()      const { return (type() == LL || type() == RL); }
    bool is_right_xx()     const { return (type() == RR || type() == RL); }
    bool is_xx_right()     const { return (type() == RR || type() == LR); }
    bool is_internal()     const { return (type() == RL || type() == LR); }
    bool is_external()     const { return (type() == LL || type() == RR); }
    // -------------------------------------------------------------------------
  std::ostream& print(std::ostream& os) const {
    switch (type()) {
	case LL : os<<"LL"; break;
	case LR : os<<"LR"; break;
	case RL : os<<"RL"; break;
	case RR : os<<"RR"; break;
	default: os<<"Unknown type ";
    }
    os<<" {"<<*source_object()<<"}{"<<*target_object()<<"}";
    return os;
  }


private :
    Type              type_;
    Disk_handle   source_object_ , target_object_;
};

// General definition of Bitangent_2.
// Most traits classes define a partial specialisation.

template< class Gtr_ >
struct Bitangent_2 
  : public Bitangent_base<typename Gtr_::Disk> 
{

    typedef Gtr_ Gt;
    typedef typename Gt::Disk Disk;
    typedef typename Gt::Point_2        Point_2;
    typedef typename Gt::Arc_2     Arc_2;
    typedef Bitangent_base<Disk>    Base;
    typedef typename Base::Disk_handle        Disk_handle;
    typedef typename Base::Type                  Type;

    Bitangent_2() : Base() { }
    Bitangent_2(const Point_2& v1 , const Point_2& v2 , Type t ,
		Disk_handle start, Disk_handle finish)
	: Base(v1,v2,t,start,finish) { }
    Bitangent_2(Type t, const Arc_2& source,const Arc_2& target) 
	: Base(t,source.object(),target.object()) { }
    Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) 
	: Base(t,o1,o2) { }
    Bitangent_2(const Bitangent_2&sibling,bool reverse,Type t)
      : Base(sibling,reverse,t) {}

    bool operator==(const Bitangent_2& b) const 
    { return Base::operator==(b); }
    bool operator!=(const Bitangent_2& b) const 
    { return Base::operator!=(b); }
};


template < class D>
std::ostream &
operator<<(std::ostream &os, const Bitangent_2<D> &b) {
  return b.print(os);
}

}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_BITANGENT_2_H
