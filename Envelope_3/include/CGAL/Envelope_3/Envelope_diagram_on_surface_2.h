// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//             Baruch Zukerman        <baruchzu@post.tau.ac.il>
//             Ophir Setter           <ophirset@post.tau.ac.il>
//             Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIAGRAM_ON_SURFACE_2_H
#define CGAL_ENVELOPE_DIAGRAM_ON_SURFACE_2_H

#include <CGAL/license/Envelope_3.h>

#include <memory>
#include <type_traits>
#include <utility>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>
#include <CGAL/Arrangement_2/arrangement_type_traits.h>
#include <CGAL/Envelope_3/Envelope_pm_dcel.h>

namespace CGAL {

/*! \class Envelope_diagram_on_surface_2
 * Representation of an envelope diagram (a minimization diagram or a
 * maximization diagram).
 */
template <typename GeomTraits_,
          typename TopolTraits_ =
            typename Default_planar_topology<GeomTraits_,
                                             Envelope_3::Envelope_pm_dcel<GeomTraits_,
                                                                          typename GeomTraits_::Xy_monotone_surface_3>
                                            >::Traits>
class Envelope_diagram_on_surface_2 :
    public Arrangement_on_surface_2<GeomTraits_, TopolTraits_> {
public:
  using Traits_3 = GeomTraits_;
  using TopolTraits = TopolTraits_;
  using Xy_monotone_surface_3 = typename Traits_3::Xy_monotone_surface_3;

protected:
  using Self = Envelope_diagram_on_surface_2<Traits_3, TopolTraits>;

  friend class Arr_accessor<Self>;

public:
  using Base = Arrangement_on_surface_2<Traits_3, TopolTraits>;

  // The following is not needed anymore, but kept for backward compatibility
  using Arrangement = Base;

  using Face = typename Base::Face;
  using Surface_iterator = typename Face::Data_iterator;
  using Surface_const_iterator = typename Face::Data_const_iterator;

  /*! A shared pointer to the (immutable) geometry traits. */
  using Shared_geometry_traits = typename Base::Shared_geometry_traits;

  /*! constructs default. */
  Envelope_diagram_on_surface_2() = default;

  /*! constructs given a shared traits object. The diagram (co-)owns the traits. */
  explicit Envelope_diagram_on_surface_2(Shared_geometry_traits tr) : Base(std::move(tr)) {}

  /*! constructs given a traits object. The caller retains ownership of the traits and must keep it alive as long as
   * the diagram (or any copy of it) exists.
   */
  Envelope_diagram_on_surface_2(const Traits_3* tr) : Base(tr) {}
};

/*! \class Envelope_diagram_2
 * Representation of a planar envelope diagram (a minimization diagram or a
 * maximization diagram).
 */
template <typename GeomTraits,
          typename Dcel_ = Envelope_3::Envelope_pm_dcel<GeomTraits, typename GeomTraits::Xy_monotone_surface_3>>
class Envelope_diagram_2 :
  public Envelope_diagram_on_surface_2<GeomTraits, typename Default_planar_topology<GeomTraits, Dcel_>::Traits> {
public:
  using Traits_3 = GeomTraits;
  using Xy_monotone_surface_3 = typename Traits_3::Xy_monotone_surface_3;

protected:
  using Env_dcel = Dcel_;
  using Self = Envelope_diagram_2<Traits_3, Env_dcel>;

  friend class Arr_accessor<Self>;

public:
  using Topology_traits = typename Default_planar_topology<Traits_3, Env_dcel>::Traits;
  using Base = Envelope_diagram_on_surface_2<Traits_3, Topology_traits>;
  using Surface_iterator = typename Base::Surface_iterator;
  using Surface_const_iterator = typename Base::Surface_const_iterator;

  // The following is not needed anymore, but kept for backward compatibility
  using Arrangement = typename Base::Base;

  /*! A shared pointer to the (immutable) geometry traits. */
  using Shared_geometry_traits = typename Base::Shared_geometry_traits;

  /*! constructs default. */
  Envelope_diagram_2() = default;

  /*! constructs given a shared traits object. The diagram (co-)owns the traits. */
  explicit Envelope_diagram_2(Shared_geometry_traits tr) : Base(std::move(tr)) {}

  /*! constructs given a traits object. The caller retains ownership of the traits and must keep it alive as long as
   * the diagram (or any copy of it) exists.
   */
  Envelope_diagram_2(const Traits_3* tr) : Base(tr) {}
};

//-----------------------------------------------------------------------------
// Specializations of is_arrangement_2 for the envelope diagrams.
//
template <typename GeomTraits_, typename TopolTraits_>
class is_arrangement_2<Envelope_diagram_on_surface_2<GeomTraits_, TopolTraits_>> : public std::true_type {};

template <typename GeomTraits_, typename Dcel_>
class is_arrangement_2<Envelope_diagram_2<GeomTraits_, Dcel_>> : public std::true_type {};

} // namespace CGAL

#endif
