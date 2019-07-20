// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_SLIVER_CRITERIA_H
#define CGAL_MESH_3_SLIVER_CRITERIA_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/min_dihedral_angle.h>
#include <CGAL/Mesh_3/radius_ratio.h>
#include <CGAL/FPU.h> // for CGAL::IA_force_to_double
#include <vector>

namespace CGAL {

namespace Mesh_3 {

template<typename Tr, //Triangulation
         bool update_sliver_cache = true,
         typename Cell_vector_ = std::vector<typename Tr::Cell_handle> >
class Sliver_criterion
{
public:
  typedef typename Tr::Geom_traits Gt;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Gt::Tetrahedron_3 Tetrahedron_3;
  typedef Cell_vector_ Cell_vector;

public:
  virtual double get_default_value() const = 0;
  virtual double get_max_value() const = 0;
  //Sliver_perturber performs perturbation "unit-per-unit"
  // so it needs to know how much is a unit for each criterion
  virtual double get_perturbation_unit() const = 0;

  // returns the value of the criterion for t
  virtual double operator()(Cell_handle cell) const
  {
    if(update_sliver_cache)
    {
      if( ! cell->is_cache_valid() )
      {
        // cell->sliver_value() is stored in a plain 64 bits floating point
        // number, and the value computed by operator() might be 
        // computed using the 80 bits floating point registers of the x87
        // unit. 
        // This could cause comparisons between sliver values to be 
        // inconsistent when comparing caches and registers values
        // (see also the comment below)
        // IA_force_to_double is available in CGAL/FPU.h
        double value 
          = CGAL::IA_force_to_double(operator()(tr_.tetrahedron(cell)));
        cell->set_sliver_value(value);
      }
      return cell->sliver_value();
    }  
    return operator()(tr_.tetrahedron(cell));
  }

  virtual double operator()(const Tetrahedron_3& t) const = 0;

  virtual void before_move(const Cell_vector& cells) const = 0;
  virtual bool valid_move(const Cell_vector& cells, 
                          const bool soft = false) const = 0;

  // Sliver bound
  void set_sliver_bound(double bound) { sliver_bound_ = bound; }
  double sliver_bound() const { return sliver_bound_; }

public:
  Sliver_criterion(const double& bound,
                   const Tr& tr)
    : tr_(tr)
    , sliver_bound_(bound)
  {}

  virtual ~Sliver_criterion(){}

protected:
  const Tr& tr_;
  double sliver_bound_;
};

template<typename SliverCriterion, typename Cell_vector>
class Min_value;

template <typename Tr,
          bool update_sliver_cache = true>
class Min_dihedral_angle_criterion 
  : public Sliver_criterion<Tr, update_sliver_cache>
{
protected:
  typedef Sliver_criterion<Tr, update_sliver_cache> Base;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;
  typedef typename Base::Gt             Gt;
  
public:
  typedef typename Base::Cell_handle    Cell_handle;

  virtual double get_default_value() const { return 12.; }
  virtual double get_max_value() const { return 90.; }
  virtual double get_perturbation_unit() const { return 1.; }

#if ( _MSC_VER != 1800 )
  using Base::operator();
#else
  virtual double operator()(Cell_handle cell) const
  {
    return Base::operator()(cell);
  }
#endif

  virtual double operator()(const Tetrahedron_3& t) const
  {
    return CGAL::to_double(minimum_dihedral_angle(t, this->tr_.geom_traits()));
  }

  virtual void before_move(const Cell_vector& cells) const
  {
    Min_value<Min_dihedral_angle_criterion, Cell_vector> min_value_op(*this);
    min_value_before_move_ = min_value_op(cells);
  }
  virtual bool valid_move(const Cell_vector& cells,
                          const bool soft = false) const
  {
    Min_value<Min_dihedral_angle_criterion, Cell_vector> min_value_op(*this);
    double min_val = min_value_op(cells);
    return (min_val > min_value_before_move_) 
        || (soft && min_val > this->sliver_bound_);
  }

public:
  Min_dihedral_angle_criterion(const double& sliver_bound,
                               const Tr& tr)
    : Base(sliver_bound, tr), min_value_before_move_(0.)
  {}
  
private:
  mutable double min_value_before_move_;
};


template <typename Tr,
          bool update_sliver_cache = true>
class Radius_ratio_criterion 
  : public Sliver_criterion<Tr, update_sliver_cache>
{
protected:
  typedef Sliver_criterion<Tr, update_sliver_cache> Base;
  typedef typename Base::Gt             Gt;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;
  typedef Radius_ratio_criterion<Tr, update_sliver_cache> RR_criterion;
  
public:

  virtual double get_default_value() const { return 0.25; }
  virtual double get_max_value() const { return 1.; }
  virtual double get_perturbation_unit() const { return 0.05; }

  using Base::operator();

  virtual double operator()(const Tetrahedron_3& t) const
  {
    return CGAL::to_double(radius_ratio(t, this->tr_.geom_traits()));
  }

  virtual void before_move(const Cell_vector& cells) const
  {
    Min_value<RR_criterion, Cell_vector> min_value_op(*this);
    min_value_before_move_ = min_value_op(cells);
  }
  virtual bool valid_move(const Cell_vector& cells,
                          const bool soft = false) const
  {
    Min_value<RR_criterion, Cell_vector> min_value_op(*this);
    double min_val = min_value_op(cells);
    return (min_val > min_value_before_move_) 
        || (soft && min_val > this->sliver_bound_);
  }

public:
  Radius_ratio_criterion(const double& sliver_bound,
               const Tr& tr)
    : Base(sliver_bound, tr)
  {}
  
private:
  mutable double min_value_before_move_;
};


template<typename SliverCriterion, typename Cell_vector>
class Min_value
{
public:
  double operator()(const Cell_vector& cells)
  {
    double minimum = criterion_.get_max_value();
    typename Cell_vector::const_iterator it;
    for(it = cells.begin(); it != cells.end(); ++it)
    {
      typename SliverCriterion::Cell_handle c = *it;
      minimum = (std::min)(minimum, criterion_(c));
    }
    return minimum;
  }

  Min_value(const SliverCriterion& criterion)
    : criterion_(criterion)
  {}

private:
  SliverCriterion criterion_;
};
  
} // end namespace Mesh_3
  
} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_SLIVER_CRITERIA_H
