// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_SLIVERS_EXUDER_CELL_ATTRIBUTES_TRAITS_H
#define CGAL_MESH_3_SLIVERS_EXUDER_CELL_ATTRIBUTES_TRAITS_H

#include <CGAL/license/Mesh_3.h>


#include <boost/mpl/has_xxx.hpp>

namespace CGAL {
namespace Mesh_3 {

// The following macro defines a metafunction
// has_Slivers_exuder_attributes so that
// has_Slivers_exuder_attributes<T>::value is true iff
// a nested type T::Slivers_exuder_attributes exists.
BOOST_MPL_HAS_XXX_TRAIT_DEF(Slivers_exuder_attributes)

struct Empty_class {};

template <class Cell, bool>
struct Slivers_ex_att_t_aux
{
  typedef Empty_class Cell_attributes;

  Cell_attributes get_attributes(const Cell* ) const
  {
    return Cell_attributes();
  }

  void restore_attributes(const Cell*,
                                         const Cell_attributes&)
  {
  }
}; // end struct Slivers_ex_att_t_aux<Cell, bool>

template <class Cell>
struct Slivers_ex_att_t_aux<Cell, true>
{
  typedef typename Cell::Slivers_exuder_attributes Cell_attributes;

  Cell_attributes get_attributes(Cell* c) const
  {
    return c->slivers_exuder_get_attributes();
  }

  void restore_attributes(Cell* c, const Cell_attributes& attr)
  {
    return c->slivers_exuder_restore_attributes(attr);
  }
}; // end partial specialisation Slivers_ex_att_t_aux<Cell, true>

template <class Cell>
struct Slivers_exuder_cell_attributes_traits
  : public Slivers_ex_att_t_aux<Cell,
                                has_Slivers_exuder_attributes<Cell>::value >
{
};

} // end namespace Mesh_3
} // end namespace CGAL


#endif // CGAL_MESH_3_SLIVERS_EXUDER_CELL_ATTRIBUTES_TRAITS_H
