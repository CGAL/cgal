// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>
//                 Ophir Setter           <ophirset@post.tau.ac.il>
//                 Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIAGRAM_ON_SURFACE_2_H
#define CGAL_ENVELOPE_DIAGRAM_ON_SURFACE_2_H

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
template <class GeomTraits_, class TopTraits_ = 
          typename Default_planar_topology< 
                               GeomTraits_,
                               Envelope_3::Envelope_pm_dcel<
                                 GeomTraits_, 
                                 typename GeomTraits_::Xy_monotone_surface_3
                               > >::Traits
          > 
class Envelope_diagram_on_surface_2 :
public Arrangement_on_surface_2<GeomTraits_, TopTraits_>
{
public:
  typedef GeomTraits_                                   Traits_3;
  typedef TopTraits_                                    TopTraits;
  typedef typename Traits_3::Xy_monotone_surface_3      Xy_monotone_surface_3;

protected:
  typedef Envelope_3::Envelope_pm_dcel<Traits_3,
                           Xy_monotone_surface_3>       Env_dcel;
  typedef Envelope_diagram_on_surface_2<Traits_3, TopTraits>       Self;
  friend class Arr_accessor<Self>;

public:
  typedef Arrangement_on_surface_2<Traits_3, 
    TopTraits>                                          Base;
  // This is yacky, but we have not choice because of observer stuff.
  typedef Base                                          Arrangement;

  typedef typename Env_dcel::Dcel_data_const_iterator   Surface_const_iterator;

  /*! Default constructor. */
  Envelope_diagram_on_surface_2() :
    Base()
  {}

  /*! Constructor with a traits-class instance. */
  Envelope_diagram_on_surface_2 (Traits_3* tr) :
    Base (tr)
  {}

};


/*! \class
 * Representation of an envelope diagram (a minimization diagram or a
 * maximization diagram).
 */

template <class GeomTraits_, 
          class Dcel_ = Envelope_3::Envelope_pm_dcel< 
            GeomTraits_, typename GeomTraits_::Xy_monotone_surface_3
            > 
          >
class Envelope_diagram_2 :
  public Envelope_diagram_on_surface_2< GeomTraits_,
     typename Default_planar_topology<GeomTraits_,
                                      Dcel_>::Traits
    >
{
public:
  typedef GeomTraits_                                   Traits_3;
  typedef typename Traits_3::Xy_monotone_surface_3      Xy_monotone_surface_3;
  
protected:
  typedef Dcel_                                         Env_dcel;
  typedef Envelope_diagram_2<Traits_3, Env_dcel>        Self;
  friend class Arr_accessor<Self>;

public:
  typedef typename Default_planar_topology< Traits_3,
                                   Env_dcel >::Traits   Topology_traits;
  typedef Envelope_diagram_on_surface_2<Traits_3, Topology_traits>
                                                        Base;
  typedef typename Env_dcel::Dcel_data_const_iterator   Surface_const_iterator;

  // This is yacky, but we have not choice because of observer stuff.
  typedef typename Base::Base                           Arrangement;


  /*! Default constructor. */
  Envelope_diagram_2() :
    Base()
  {}

  /*! Constructor with a traits-class instance. */
  Envelope_diagram_2 (Traits_3* tr) :
    Base (tr)
  {}

};


//--------------------------------  Envelope_on_surface_3
// specialization
template <class GeomTraits_, class TopTraits_>
class is_arrangement_2< 
  Envelope_diagram_on_surface_2<GeomTraits_, TopTraits_>
> : public boost::true_type
{};

// specialization
template <class GeomTraits_, class DCEL_>
class is_arrangement_2< 
  Envelope_diagram_2<GeomTraits_, DCEL_>
> : public boost::true_type
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
