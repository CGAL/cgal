// Copyright (c) 2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

// CGAL_CFG_CONVERSION_OPERATOR_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_CONVERSION_OPERATOR_BUG is set, if a compiler
//| crashes with some conversion operators.  G++ 3.3.0 is affected by
//| this bug (it hits Darwin severely since it is the system compiler).

template<class I1, class I2>
struct iterator_restrict_traits {
  typedef I1 iterator_category;
};

template<class T>
struct scalar_expression {
  typedef T value_type;
};

template<class E, class F>
class vector_scalar_unary:
  public scalar_expression<typename F::result_type> {
public:
  typedef typename F::result_type value_type;
  typedef typename E::const_iterator::iterator_category iterator_category;

  operator value_type () const {
    return evaluate (iterator_category ());
  }
};

template<class E1, class E2, class F>
class vector_scalar_binary:
  public scalar_expression<typename F::result_type> {
public:
  typedef typename F::result_type value_type;
  typedef typename iterator_restrict_traits<
	    typename E1::const_iterator::iterator_category,
	    typename E2::const_iterator::iterator_category>::iterator_category
iterator_category;

  operator value_type () const {
    return evaluate (iterator_category ());
  }
};

int main() {
    return 0;
}
