// Copyright (c) 2003  Utrecht University (The Netherlands),
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
// Author(s)     : various

// CGAL_CFG_USING_NAMESPACE_BUG.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_USING_NAMESPACE_BUG is set, if a compiler cannot
//| handle using directives (using namespace ...) across different
//| namespaces. Created to workaround a cl1300 bug

namespace CGAL {
  namespace CommonFunctors {
    template < class K >
    struct F {};
  }

  namespace CartesianFunctors {
    using namespace CommonFunctors;
  }

  template < class FT >
  struct Cartesian {
    typedef Cartesian<FT>  Self;
    typedef CartesianFunctors::F<Self> Func;
  };

} // end namespace CGAL

int main() {
  CGAL::Cartesian<double>::Func f;
  (void) f;
  return 0;
}
