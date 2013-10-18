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
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_SLIVER_CRITERIA_H
#define CGAL_MESH_3_SLIVER_CRITERIA_H

#include <CGAL/Mesh_3/min_dihedral_angle.h>
#include <CGAL/Mesh_3/radius_ratio.h>
#include <vector>

namespace CGAL {

namespace Mesh_3 {

template<typename Tr, //Triangulation
         typename Cell_vector = std::vector<typename Tr::Cell_handle> >
class Sliver_criterion
{
public:
  typedef typename Tr::Geom_traits K;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename K::Tetrahedron_3 Tetrahedron_3;
  typedef Cell_vector Cell_vector;

public:
  virtual const double get_max_value() const = 0;
  //Sliver_perturber performs perturbation "unit-per-unit"
  // so it needs to know how much is a unit for each criterion
  virtual const double get_perturbation_unit() const = 0;

  // returns the value of the criterion for t
  virtual double operator()(Cell_handle cell) const
  {
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
    : sliver_bound_(bound),
      tr_(tr)
  {}

protected:
  const Tr& tr_;
  double sliver_bound_;
};
  
template <typename Tr>
class Min_dihedral_angle_criterion : public Sliver_criterion<Tr>
{
protected:
  typedef Sliver_criterion<Tr> Base;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;
  typedef Min_dihedral_angle_criterion<Tr> Dihedral_angle_criterion;
  
public:
  static double default_value;
  static double max_value;
  static double min_value;

  virtual const double get_max_value() const { return 90.; }
  virtual const double get_perturbation_unit() const { return 1.; }

  using Base::operator();

  virtual double operator()(const Tetrahedron_3& t) const
  {
    return CGAL::to_double(minimum_dihedral_angle(t, K()));
  }

  virtual void before_move(const Cell_vector& cells) const
  {
    Min_value<Dihedral_angle_criterion, Cell_vector> min_value_op(*this);
    min_value_before_move_ = min_value_op(cells);
  }
  virtual bool valid_move(const Cell_vector& cells,
                          const bool soft = false) const
  {
    Min_value<Dihedral_angle_criterion, Cell_vector> min_value_op(*this);
    double min_val = min_value_op(cells);
    return (min_val > min_value_before_move_) 
        || (soft && min_val > sliver_bound_);
  }

public:
  Dihedral_angle_criterion(const double& sliver_bound,
                           const Tr& tr)
    : Base(sliver_bound, tr)
  {}
  
private:
  mutable double min_value_before_move_;
};

template<typename Tr> double Min_dihedral_angle_criterion<Tr>::default_value = 12.;
template<typename Tr> double Min_dihedral_angle_criterion<Tr>::max_value = 90.; 
template<typename Tr> double Min_dihedral_angle_criterion<Tr>::min_value = 0.; 
  
template <typename Tr>
class Radius_ratio_criterion : public Sliver_criterion<Tr>
{
protected:
  typedef Sliver_criterion<Tr> Base;
  typedef typename Base::K              K;
  typedef typename Base::Tetrahedron_3  Tetrahedron_3;
  typedef typename Base::Cell_vector    Cell_vector;
  typedef Radius_ratio_criterion<Tr> RR_criterion;
  
public:
  static double default_value;
  static double max_value;
  static double min_value;

  virtual const double get_max_value() const { return 1.; }
  virtual const double get_perturbation_unit() const { return 0.05; }

  using Base::operator();

  virtual double operator()(const Tetrahedron_3& t) const
  {
    return CGAL::to_double(radius_ratio(t, K()));
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
        || (soft && min_val > sliver_bound_);
  }

public:
  RR_criterion(const double& sliver_bound,
               const Tr& tr)
    : Base(sliver_bound, tr)
  {}
  
private:
  mutable double min_value_before_move_;
};

template<typename Tr> double Radius_ratio_criterion<Tr>::default_value = 0.25; 
template<typename Tr> double Radius_ratio_criterion<Tr>::max_value = 1.;
template<typename Tr> double Radius_ratio_criterion<Tr>::min_value = 0.; 


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




#endif // CGAL_MESH_3_SLIVER_CRITERIA_H
