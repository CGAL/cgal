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


#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>
#include <CGAL/Arrangement_2/arrangement_type_traits.h>
#include <CGAL/Envelope_3/Envelope_pm_dcel.h>

namespace CGAL {

/*! \class
 * Representation of an envelope diagram (a minimization diagram or a
 * maximization diagram).
 */
template <typename GeomTraits_, typename TopTraits_ =
          typename Default_planar_topology
            <GeomTraits_,
             Envelope_3::Envelope_pm_dcel
               <GeomTraits_,
                typename GeomTraits_::Xy_monotone_surface_3>>::Traits>
class Envelope_diagram_on_surface_2 :
    public Arrangement_on_surface_2<GeomTraits_, TopTraits_> {
public:
  using Traits_3 = GeomTraits_;
  using TopTraits = TopTraits_;
  using Xy_monotone_surface_3 = typename Traits_3::Xy_monotone_surface_3;

protected:
  using Self = Envelope_diagram_on_surface_2<Traits_3, TopTraits>;

  friend class Arr_accessor<Self>;

public:
  using Base = Arrangement_on_surface_2<Traits_3, TopTraits>;

  // The following is not needed anymore, but kept for backward compatibility
  using Arrangement = Base;

  using Face = typename Base::Face;
  using Surface_iterator = typename Face::Data_iterator;
  using Surface_const_iterator = typename Face::Data_const_iterator;

  /*! Default constructor. */
  Envelope_diagram_on_surface_2() : Base() {}

  /*! Constructor with a traits-class instance. */
  Envelope_diagram_on_surface_2(const Traits_3* tr) : Base(tr) {}
};

/*! \class
 * Representation of an envelope diagram (a minimization diagram or a
 * maximization diagram).
 */
template <typename GeomTraits,
          typename Dcel_ = Envelope_3::Envelope_pm_dcel
            <GeomTraits, typename GeomTraits::Xy_monotone_surface_3>>
class Envelope_diagram_2 :
  public Envelope_diagram_on_surface_2<GeomTraits,
                                       typename Default_planar_topology
                                         <GeomTraits, Dcel_>::Traits>
{
public:
  using Traits_3 = GeomTraits;
  using Xy_monotone_surface_3 = typename Traits_3::Xy_monotone_surface_3;

protected:
  using Env_dcel = Dcel_;
  using Self = Envelope_diagram_2<Traits_3, Env_dcel>;

  friend class Arr_accessor<Self>;

public:
  using Topology_traits =
    typename Default_planar_topology< Traits_3, Env_dcel>::Traits;
  using Base = Envelope_diagram_on_surface_2<Traits_3, Topology_traits>;
  using Surface_iterator = typename Base::Surface_iterator;
  using Surface_const_iterator = typename Base::Surface_const_iterator;

  // The following is not needed anymore, but kept for backward compatibility
  using Arrangement = typename Base::Base;

  /*! Default constructor. */
  Envelope_diagram_2() : Base() {}

  /*! Constructor with a traits-class instance. */
  Envelope_diagram_2(const Traits_3* tr) : Base(tr) {}
};

//--------------------------------  Envelope_on_surface_3
// specialization
template <typename GeomTraits_, typename TopTraits_>
class is_arrangement_2<Envelope_diagram_on_surface_2<GeomTraits_, TopTraits_>> :
    public std::true_type
{};

// specialization
template <typename GeomTraits_, typename Dcel_>
class is_arrangement_2<Envelope_diagram_2<GeomTraits_, Dcel_>> :
    public std::true_type
{};

// /*! \class
//  * Representation of an envelope diagram (a minimization diagram or a
//  * maximization diagram).
//  */
// template <typename T_Traits,
// #ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
//           template <class T1, class T2>
// #endif
//           class T_Dcel = Envelope_3::Envelope_pm_dcel>
// class Envelope_diagram_2 :
//   public Arrangement_2<T_Traits,
// #ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
//                        T_Dcel<T_Traits,
//                               typename T_Traits::Xy_monotone_surface_3>
// #else
//                        typename T_Dcel::template Dcel<T_Traits,
//                                                       typename T_Traits::Xy_monotone_surface_3>
// #endif
//                        >
// {
// public:
//   typedef T_Traits                                      Traits_3;
//   typedef typename Traits_3::Xy_monotone_surface_3      Xy_monotone_surface_3;

// protected:
//   typedef T_Dcel<Traits_3, Xy_monotone_surface_3>       Env_dcel;
//   typedef Envelope_diagram_2<Traits_3, T_Dcel>          Self;
//   friend class Arr_accessor<Self>;

// public:
//   typedef Arrangement_2<Traits_3, Env_dcel>             Base;
//   typedef typename Env_dcel::Dcel_data_const_iterator   Surface_const_iterator;

//   /*! Default constructor. */
//   Envelope_diagram_2() :
//     Base()
//   {}

//   /*! Constructor with a traits-class instance. */
//   Envelope_diagram_2 (Traits_3* tr) :
//     Base (tr)
//   {}

// };

} //namespace CGAL

#endif
