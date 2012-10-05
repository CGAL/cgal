// Copyright (c) 2003,2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_FUNCTION_OBJECTS_H
#define CGAL_FUNCTION_OBJECTS_H 1

#include <functional>

#include <CGAL/enum.h>

namespace CGAL {

template < class Value>
struct Identity {
  typedef Value argument_type;
  typedef Value result_type;
  Value&       operator()( Value& x)       const { return x; }
  const Value& operator()( const Value& x) const { return x; }
};

template < class Value>
struct Dereference {
  typedef Value* argument_type;
  typedef Value  result_type;
  Value&         operator()( Value* x)       const { return *x;}
  const Value&   operator()( const Value* x) const { return *x;}
};

template < class Value>
struct Get_address {
  typedef Value  argument_type;
  typedef Value* result_type;
  Value*         operator()( Value& x)       const { return &x;}
  const Value*   operator()( const Value& x) const { return &x;}
};

template < class Arg, class Result>
struct Cast_function_object {
  typedef Arg    argument_type;
  typedef Result result_type;
  Result&       operator()( Arg& x)       const { return (Result&)(x); }
  const Result& operator()( const Arg& x) const {
    return (const Result&)(x);
  }
};

template < class Node>
struct Project_vertex {
  typedef Node                  argument_type;
  typedef typename Node::Vertex Vertex;
  typedef Vertex                result_type;
  Vertex&       operator()( Node& x)       const { return x.vertex(); }
  const Vertex& operator()( const Node& x) const { return x.vertex(); }
};

template < class Node>
struct Project_facet {
  typedef Node                  argument_type;
  typedef typename Node::Facet  Facet;
  typedef Facet                 result_type;
  Facet&       operator()( Node& x)       const { return x.facet(); }
  const Facet& operator()( const Node& x) const { return x.facet(); }
};

template < class Node>
struct Project_point {
  typedef Node                  argument_type;
  typedef typename Node::Point  Point;
  typedef Point                 result_type;
  Point&       operator()( Node& x)       const { return x.point(); }
  const Point& operator()( const Node& x) const { return x.point(); }
};

template < class Node>
struct Project_normal {
  typedef Node                  argument_type;
  typedef typename Node::Normal Normal;
  typedef Normal                result_type;
  Normal&       operator()( Node& x)       const { return x.normal(); }
  const Normal& operator()( const Node& x) const { return x.normal(); }
};

template < class Node>
struct Project_plane {
  typedef Node                  argument_type;
  typedef typename Node::Plane  Plane;
  typedef Plane                 result_type;
  Plane&       operator()( Node& x)       const { return x.plane(); }
  const Plane& operator()( const Node& x) const { return x.plane(); }
};

// The following four adaptors are used to create the circulators
// for polyhedral surfaces.
template < class Node>
struct Project_next {
  typedef Node*   argument_type;
  typedef Node*   result_type;
  Node*       operator()( Node* x)       const { return x->next(); }
  const Node* operator()( const Node* x) const { return x->next(); }
};

template < class Node>
struct Project_prev {
  typedef Node*   argument_type;
  typedef Node*   result_type;
  Node*       operator()( Node* x)       const { return x->prev(); }
  const Node* operator()( const Node* x) const { return x->prev(); }
};

template < class Node>
struct Project_next_opposite {
  typedef Node*   argument_type;
  typedef Node*   result_type;
  Node*       operator()( Node* x)       const {
    return x->next()->opposite();
  }
  const Node* operator()( const Node* x) const {
    return x->next()->opposite();
  }
};

template < class Node>
struct Project_opposite_prev {
  typedef Node*   argument_type;
  typedef Node*   result_type;
  Node*       operator()( Node* x)       const {
    return x->opposite()->prev();
  }
  const Node* operator()( const Node* x) const {
    return x->opposite()->prev();
  }
};
template < class Arg, class  Result >
class Creator_1 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Result result_type;
  Result operator()(Arg a) const { return Result(a);}
};

template < class Arg1, class Arg2, class  Result >
class Creator_2 {
public:
  typedef Arg1   argument1_type;
  typedef Arg2   argument2_type;
  typedef Result result_type;
  Result operator()(Arg1 a1, Arg2 a2) const { return Result(a1,a2);}
};

template < class Arg1, class Arg2, class Arg3, class  Result >
class Creator_3 {
public:
  typedef Arg1   argument1_type;
  typedef Arg2   argument2_type;
  typedef Arg3   argument3_type;
  typedef Result result_type;
  Result operator()(Arg1 a1, Arg2 a2, Arg3 a3) const {
    return Result(a1,a2,a3);
  }
};

template < class Arg1, class Arg2, class Arg3, class Arg4, class  Result >
class Creator_4 {
public:
  typedef Arg1   argument1_type;
  typedef Arg2   argument2_type;
  typedef Arg3   argument3_type;
  typedef Arg4   argument4_type;
  typedef Result result_type;
  Result operator()(Arg1 a1, Arg2 a2, Arg3 a3, Arg4 a4) const {
    return Result(a1,a2,a3,a4);
  }
};

template < class Arg1, class Arg2, class Arg3, class Arg4, class Arg5,
           class  Result >
class Creator_5 {
public:
  typedef Arg1   argument1_type;
  typedef Arg2   argument2_type;
  typedef Arg3   argument3_type;
  typedef Arg4   argument4_type;
  typedef Arg5   argument5_type;
  typedef Result result_type;
  Result operator()(Arg1 a1, Arg2 a2, Arg3 a3, Arg4 a4, Arg5 a5) const {
    return Result(a1,a2,a3,a4,a5);
  }
};

template < class Arg, class  Result >
class Creator_uniform_2 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2) const { return Result(a1,a2);}
};

template < class Arg, class  Result >
class Creator_uniform_3 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3) const {
    return Result(a1,a2,a3);
  }
};

template < class Arg, class  Result >
class Creator_uniform_4 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4) const {
    return Result(a1,a2,a3,a4);
  }
};

template < class Arg, class  Result >
class Creator_uniform_5 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Arg    argument5_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4, Arg a5) const {
    return Result(a1,a2,a3,a4,a5);
  }
};

template < class Arg, class  Result >
class Creator_uniform_6 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Arg    argument5_type;
  typedef Arg    argument6_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4, Arg a5, Arg a6
  ) const {
    return Result(a1,a2,a3,a4,a5,a6);
  }
};

template < class Arg, class  Result >
class Creator_uniform_7 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Arg    argument5_type;
  typedef Arg    argument6_type;
  typedef Arg    argument7_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4, Arg a5, Arg a6,
                    Arg a7) const {
                      return Result(a1,a2,a3,a4,a5,a6,a7);
                    }
};

template < class Arg, class  Result >
class Creator_uniform_8 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Arg    argument5_type;
  typedef Arg    argument6_type;
  typedef Arg    argument7_type;
  typedef Arg    argument8_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4, Arg a5, Arg a6,
                    Arg a7, Arg a8) const {
                      return Result(a1,a2,a3,a4,a5,a6,a7,a8);
                    }
};

template < class Arg, class  Result >
class Creator_uniform_9 {
public:
  typedef Arg    argument_type;
  typedef Arg    argument1_type;
  typedef Arg    argument2_type;
  typedef Arg    argument3_type;
  typedef Arg    argument4_type;
  typedef Arg    argument5_type;
  typedef Arg    argument6_type;
  typedef Arg    argument7_type;
  typedef Arg    argument8_type;
  typedef Arg    argument9_type;
  typedef Result result_type;
  Result operator()(Arg a1, Arg a2, Arg a3, Arg a4, Arg a5, Arg a6,
                    Arg a7, Arg a8, Arg a9) const {
                      return Result(a1,a2,a3,a4,a5,a6,a7,a8,a9);
                    }
};

template < class Arg, class  Result >
class Creator_uniform_d {
  int d;

 private:
  Creator_uniform_d(){}

 public:
  typedef Arg   argument1_type;
  typedef Result result_type;

  Creator_uniform_d(int dim)
    : d(dim)
    {}

  Result operator()(Arg a1, Arg a2) const { return Result(d, a1,a2);}
};

template < class Op1, class Op2 >
class Unary_compose_1
  : public std::unary_function< typename Op2::argument_type,
                                typename Op1::result_type >
{
protected:
  Op1 op1;
  Op2 op2;
public:
  typedef typename Op2::argument_type   argument_type;
  typedef typename Op1::result_type     result_type;

  Unary_compose_1(const Op1& x, const Op2& y) : op1(x), op2(y) {}

  result_type
  operator()(const argument_type& x) const
  { return op1(op2(x)); }
};

template < class Op1, class Op2 >
inline Unary_compose_1< Op1, Op2 >
compose1_1(const Op1& op1, const Op2& op2)
{ return Unary_compose_1< Op1, Op2 >(op1, op2); }

template < class Op1, class Op2, class Op3 >
class Binary_compose_1
  : public std::unary_function< typename Op2::argument_type,
                                typename Op1::result_type >
{
protected:
  Op1 op1;
  Op2 op2;
  Op3 op3;
public:
  typedef typename Op2::argument_type  argument_type;
  typedef typename Op1::result_type    result_type;

  Binary_compose_1(const Op1& x, const Op2& y, const Op3& z)
  : op1(x), op2(y), op3(z) {}

  result_type
  operator()(const argument_type& x) const
  { return op1(op2(x), op3(x)); }
};

template < class Op1, class Op2, class Op3 >
inline Binary_compose_1< Op1, Op2, Op3 >
compose2_1(const Op1& op1, const Op2& op2, const Op3& op3)
{ return Binary_compose_1< Op1, Op2, Op3 >(op1, op2, op3); }

template < class Op1, class Op2 >
class Unary_compose_2
  : public std::binary_function< typename Op2::first_argument_type,
                                 typename Op2::second_argument_type,
                                 typename Op1::result_type >
{
protected:
  Op1 op1;
  Op2 op2;
public:
  typedef typename Op2::first_argument_type   first_argument_type;
  typedef typename Op2::second_argument_type  second_argument_type;
  typedef typename Op1::result_type           result_type;

  Unary_compose_2(const Op1& x, const Op2& y) : op1(x), op2(y) {}

  result_type
  operator()(const first_argument_type& x,
             const second_argument_type& y) const
  { return op1(op2(x, y)); }
};

template < class Op1, class Op2 >
inline Unary_compose_2< Op1, Op2 >
compose1_2(const Op1& op1, const Op2& op2)
{ return Unary_compose_2< Op1, Op2 >(op1, op2); }

template < class Op1, class Op2, class Op3 >
class Binary_compose_2
  : public std::binary_function< typename Op2::argument_type,
                                 typename Op3::argument_type,
                                 typename Op1::result_type >
{
protected:
  Op1 op1;
  Op2 op2;
  Op3 op3;
public:
  typedef typename Op2::argument_type  first_argument_type;
  typedef typename Op3::argument_type  second_argument_type;
  typedef typename Op1::result_type    result_type;

  Binary_compose_2(const Op1& x, const Op2& y, const Op3& z)
  : op1(x), op2(y), op3(z) {}

  result_type
  operator()(const first_argument_type& x,
             const second_argument_type& y) const
  { return op1(op2(x), op3(y)); }
};

template < class Op1, class Op2, class Op3 >
inline Binary_compose_2< Op1, Op2, Op3 >
compose2_2(const Op1& op1, const Op2& op2, const Op3& op3)
{ return Binary_compose_2< Op1, Op2, Op3 >(op1, op2, op3); }


template < class Op >
class Compare_to_less
  : public Op
{
public:
  typedef Op      Type;
  typedef bool    result_type;

  Compare_to_less(const Op& op) : Op(op) {}

  template < typename A1 >
  bool
  operator()(const A1 &a1) const
  { return Op::operator()(a1) == SMALLER; }

  template < typename A1, typename A2 >
  bool
  operator()(const A1 &a1, const A2 &a2) const
  { return Op::operator()(a1, a2) == SMALLER; }

  template < typename A1, typename A2, typename A3 >
  bool
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
  { return Op::operator()(a1, a2, a3) == SMALLER; }

  template < typename A1, typename A2, typename A3, typename A4 >
  bool
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
  { return Op::operator()(a1, a2, a3, a4) == SMALLER; }

  // More can be added.
};

template < class Op >
inline Compare_to_less<Op>
compare_to_less(const Op& op)
{ return Compare_to_less<Op>(op); }


/*!\brief
 * Functor class to determine lexicographical order of pairs
 */
template < class T1, class T2,
           class Less1 = std::less<T1>, class Less2 = std::less<T2> >
struct Pair_lexicographical_less_than {
    typedef bool result_type;
    typedef std::pair<T1,T2> first_argument_type;
    typedef std::pair<T1,T2> second_argument_type;
    /*!\brief
     * returns \c true iff \c x is lexicograhically smaller than \c y
     * using \c Less1 and \c Less2 functors.
     */
    bool operator () (const std::pair<T1,T2>& x, const std::pair<T1,T2>& y) const {
        Less1 lt1;
        Less2 lt2;

        if (lt1(x.first, y.first)) {
            return true;
        } else if (lt1(y.first, x.first)) {
            return false;
        } else /* neither is less than the other */ {
            return lt2(x.second, y.second);
        }
    }
};


} //namespace CGAL

#endif // CGAL_FUNCTION_OBJECTS_H
