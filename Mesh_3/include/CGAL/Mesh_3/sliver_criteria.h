// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <CGAL/Mesh_3/Sliver_value_cache.h>

#include <CGAL/FPU.h> // for CGAL::IA_force_to_double
#include <CGAL/Default.h>
#include <CGAL/Has_member.h>

#include <vector>

namespace CGAL {
namespace Mesh_3 {

enum class Sliver_caching_policy
{
  DEFAULT = 0,
  ALWAYS,
  NEVER,
  BELOW_BOUND // Only cache values below the sliver bound
};

template<typename Tr,
         typename Cache = typename CGAL::Mesh_3::Default_sliver_cache<Tr>::type,
         Sliver_caching_policy caching_policy_ = Sliver_caching_policy::DEFAULT,
         typename Cell_vector_ = std::vector<typename Tr::Cell_handle> >
class Sliver_criterion
{
public:
  typedef typename Tr::Geom_traits GT;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename GT::Tetrahedron_3 Tetrahedron_3;
  typedef Cell_vector_ Cell_vector;

public:
  virtual double get_default_value() const = 0;
  virtual double get_max_value() const = 0;
  // Sliver_perturber performs perturbation "unit-per-unit"
  // so it needs to know how much is a unit for each criterion
  virtual double get_perturbation_unit() const = 0;

  // Reset the entire cache (if any)
  void reset_cache() const {
    cache_.clear();
  }

  // Reset the cache for a single cell (if any)
  void reset_cache(const Cell_handle& cell) const {
    cache_.invalidate(cell);
  }

  void set_cache(const Cell_handle& cell, const double value) const {
    cache_.set(cell, value);
  }

  bool is_sliver(Cell_handle cell) const
  {
    if (caching_policy == Sliver_caching_policy::BELOW_BOUND) {
      return cache_.has_value(cell);
    } else {
      // if we cache everything, the operator() will return the cached value
      return operator()(cell) < this->sliver_bound_;
    }
  }

  // returns the value of the criterion for t
  virtual double operator()(Cell_handle cell) const
  {
    if (caching_policy == Sliver_caching_policy::NEVER) {
      return operator()(tr_.tetrahedron(cell));
    } else {
      if (cache_.has_value(cell)) {
        return cache_.get(cell);
      } else {
        // cell->sliver_value() is stored in a plain 64 bits floating point
        // number, and the value computed by operator() might be
        // computed using the 80 bits floating point registers of the x87
        // unit.
        // This could cause comparisons between sliver values to be
        // inconsistent when comparing caches and registers values
        // (see also the comment below)
        // IA_force_to_double is available in CGAL/FPU.h
        const double value = CGAL::IA_force_to_double(operator()(tr_.tetrahedron(cell)));
        if (caching_policy == Sliver_caching_policy::BELOW_BOUND && value > this->sliver_bound_) {
          return value;
        } else {
          cache_.set(cell, value);
          return value;
        }
      }
    }
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
                   const Tr& tr,
                   const Cache& cache = Cache())
    : tr_(tr)
    , sliver_bound_(bound)
    , cache_(cache)
  {
    caching_policy =
      (caching_policy_ == Sliver_caching_policy::DEFAULT)
        ? (has_set_sliver_value<typename Tr::Cell>::value ? Sliver_caching_policy::ALWAYS
                                                          : Sliver_caching_policy::BELOW_BOUND)
        : caching_policy_;
  }

  virtual ~Sliver_criterion(){}

protected:
  const Tr& tr_;
  double sliver_bound_;

  CGAL_GENERATE_MEMBER_DETECTOR(set_sliver_value);
  Sliver_caching_policy caching_policy;
  mutable Cache cache_;
};

template<typename SliverCriterion, typename Cell_vector>
class Min_value;

template <typename Tr,
          typename Cache = typename CGAL::Mesh_3::Default_sliver_cache<Tr>::type,
          Sliver_caching_policy caching_policy = Sliver_caching_policy::BELOW_BOUND>
class Min_dihedral_angle_criterion
  : public Sliver_criterion<Tr, Cache, caching_policy>
{
protected:
  typedef Sliver_criterion<Tr, Cache, caching_policy> Base;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;
  typedef typename Base::GT             GT;

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
                               const Tr& tr,
                               const Cache& cache = Cache())
    : Base(sliver_bound, tr, cache), min_value_before_move_(0.)
  {}

private:
  mutable double min_value_before_move_;
};


template <typename Tr,
          typename Cache = typename CGAL::Mesh_3::Default_sliver_cache<Tr>::type,
          Sliver_caching_policy caching_policy = Sliver_caching_policy::BELOW_BOUND>
class Radius_ratio_criterion
  : public Sliver_criterion<Tr, Cache, caching_policy>
{
  typedef Radius_ratio_criterion<Tr, Cache, caching_policy> Self;

protected:
  typedef Sliver_criterion<Tr, Cache, caching_policy> Base;
  typedef typename Base::GT             GT;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;

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
    Min_value<Self, Cell_vector> min_value_op(*this);
    min_value_before_move_ = min_value_op(cells);
  }
  virtual bool valid_move(const Cell_vector& cells,
                          const bool soft = false) const
  {
    Min_value<Self, Cell_vector> min_value_op(*this);
    double min_val = min_value_op(cells);
    return (min_val > min_value_before_move_)
        || (soft && min_val > this->sliver_bound_);
  }

public:
  Radius_ratio_criterion(const double& sliver_bound,
                         const Tr& tr,
                         const Cache& cache = Cache())
    : Base(sliver_bound, tr, cache)
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
