// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
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
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H
#define CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H

CGAL_BEGIN_NAMESPACE

namespace QP_solver_impl {
  
  template<class Rep>
  class Unbounded_direction_iterator {
  private:
    QP_solver<Rep> &solver;               // reference to solver
    int pos;                              // current position
    
  public:  // public types
    typedef Unbounded_direction_iterator<Rep> Self;
    typedef std::random_access_iterator_tag   iterator_category;
    typedef typename Rep::ET                  value_type;
    typedef int                               difference_type;
    
  private: // constructor
    friend class QP_solver<Rep>;
    
    // Note: The solver's routines unbounded_direction_begin/end()
    // construct instances of this class such that the range spanned
    // by the constructed iterators has length qp_n.
    Unbounded_direction_iterator(QP_solver<Rep>& solver,int pos) :
      solver(solver), pos(pos) {
      CGAL_qpe_assertion(pos <= solver.qp_n);
    }
    
  public:
    bool operator==(const Self& i) const
    { return pos == i.pos; }

    bool operator!=(const Self& i) const
    { return !(*this == i); }

    bool operator< (const Self& i) const
    { return pos < i.pos; }

    value_type operator*() const
    { 
      // Note: the vector we return here is described in documentation/
      // Test_suite.tex.
      CGAL_qpe_assertion(pos < solver.qp_n);

      const int i = solver.in_B[pos];
      if (i >= 0) {                  // basic variable?
	return solver.q_x_O[i];
      } else {                       // non-basic variable?
	if (pos == solver.j)         // most recent entering variable?
	  return -solver.d;
	return solver.et0;
      }
      return solver.et0;
    }
  
    Self& operator++(   ) { ++pos; return *this; }
    Self  operator++(int) { Self tmp = *this; ++(*this); return tmp; }
    Self& operator--(   ) { --pos; return *this; }
    Self  operator--(int) { Self tmp = *this; --(*this); return tmp; }
  
    value_type operator[](difference_type i) const
    {
      // todo: could be optimized slightly.
      return *(*this + i);
    }
  
    Self& operator+=(difference_type n) { pos += n; return *this; }
    Self& operator-=(difference_type n) { pos -= n; return *this; }
    Self  operator+ (difference_type n) const {
      Self tmp = *this;
      return tmp += n;
    }
    Self  operator- (difference_type n) const {
      Self tmp = *this;
      return tmp -= n;
    }
    difference_type operator-(const Self& i) const { return pos - i.pos; }
  };

} // namespace QP_solver_impl

template<class Rep>
typename QP_solver<Rep>::Unbounded_direction_iterator
QP_solver<Rep>::unbounded_direction_begin()
{
  return QP_solver_impl::Unbounded_direction_iterator<Rep>(*this,0);
}

template<class Rep>
typename QP_solver<Rep>::Unbounded_direction_iterator
QP_solver<Rep>::unbounded_direction_end()
{
  return QP_solver_impl::Unbounded_direction_iterator<Rep>(*this,qp_n);
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H

// ===== EOF ==================================================================
