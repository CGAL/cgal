// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CELL_CONST_ITERATORS_H
#define CGAL_CELL_CONST_ITERATORS_H 1

#include <CGAL/Cell_iterators.h>

// TODO do all the orbit iterator of any orbit ?

namespace CGAL {

  /** @file Cell_const_iterators.h
   * The cell const iterators. There are 3 classes:
   *  - CMap_cell_const_iterator<Map,Ite,i,dim>
   *  - CMap_one_dart_per_incident_cell_const_iterator<Map,Ite,i,dim>
   *  - CMap_one_dart_per_cell_const_iterator<Map,Ite,i,dim>
   */

  //****************************************************************************
  template <typename Map_,typename Ite,
            unsigned int i,unsigned int dim=Map_::dimension>
  class CMap_cell_const_iterator: public CMap_cell_iterator<Map_,Ite,i,dim,true>
  {
  public:
    typedef CMap_cell_iterator<Map_,Ite,i,dim,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /// Main constructor.
    CMap_cell_const_iterator(const Map_& amap,
                             Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_cell_const_iterator
    (const CMap_cell_iterator<Map_,Ite,i,dim,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template <typename Map_,unsigned int i,unsigned int j,
            unsigned int dim=Map_::dimension>
  class CMap_one_dart_per_incident_cell_const_iterator:
    public CMap_one_dart_per_incident_cell_iterator<Map_,i,j,dim,true>
  {
  public:
    typedef CMap_one_dart_per_incident_cell_iterator<Map_,i,j,dim,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /// Main constructor.
    CMap_one_dart_per_incident_cell_const_iterator(const Map_& amap,
                                                   Dart_const_handle adart):
      Base(amap, adart)
    {}
    /// Constructor from non const version.
    CMap_one_dart_per_incident_cell_const_iterator
    (const CMap_one_dart_per_incident_cell_iterator<Map_,i,j,dim,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template <typename Map_,unsigned int i,unsigned int dim=Map_::dimension>
  class CMap_one_dart_per_cell_const_iterator:
    public CMap_one_dart_per_cell_iterator<Map_,i,dim,true>
  {
  public:
    typedef CMap_one_dart_per_cell_iterator<Map_,i,dim,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /// Main constructor.
    CMap_one_dart_per_cell_const_iterator(const Map_& amap): Base(amap)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    CMap_one_dart_per_cell_const_iterator(const Map_& amap,
                                          Dart_const_handle adart):
      Base(amap, adart)
    {}
    /// Constructor from non const version.
    CMap_one_dart_per_cell_const_iterator
    (const CMap_one_dart_per_cell_iterator<Map_,i,dim,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_CELL_CONST_ITERATORS_H
//******************************************************************************
