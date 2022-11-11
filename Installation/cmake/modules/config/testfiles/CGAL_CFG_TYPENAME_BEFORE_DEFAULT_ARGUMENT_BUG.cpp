// Copyright (c) 2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>

struct Class {
  struct Nested_class {
    Nested_class(){}
  };
};


template <class Cl>
void
function(const Cl& cl,
         typename Cl::Nested_class start =  typename Cl::Nested_class())
{}


int
main()
{
  Class cl;
  function(cl);
  return 0;
}
