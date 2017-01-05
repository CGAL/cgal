// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_BASE_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_BASE_H

#include <vector>

namespace CGAL {

namespace Classification {
  
/*!
  \ingroup PkgClassification

  \brief Abstract class describing a classification attribute that
  associates a scalar value to each item of the classification input.
*/


class Attribute_base
{
  double m_weight; 
public:

  /// \cond SKIP_IN_MANUAL
  double mean;
  double max;

  virtual ~Attribute_base() { }
  /// \endcond

  /*!
    \brief Returns the weight of the attribute between 0 and +&infin;.
  */
  double& weight() { return m_weight; }
  
  /*!
    \brief Returns `abstract_attribute` and should be overloaded
    by a inherited class with a specific name.
  */
  virtual std::string name() { return "abstract_attribute"; }

  /*!
    \brief Returns the value taken by the attribute at the given
    index. This method must be implemented by inherited classes.
  */
  virtual double value (std::size_t index) = 0;

  /// \cond SKIP_IN_MANUAL
  virtual double normalized (std::size_t index)
  {
    return (std::max) (0., (std::min) (1., value(index) / m_weight));
  }
  virtual double favored (std::size_t index) { return (1. - normalized (index)); }
  virtual double penalized (std::size_t index) { return normalized (index); }
  //  virtual double ignored (std::size_t index) { return (std::min) (favored(index), penalized(index)); }
  virtual double ignored (std::size_t) { return 0.5; }

  void compute_mean_max (std::vector<double>& vect, double& mean, double& max)
  {
    mean = 0.;
    max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    
    for (std::size_t i = 0; i < vect.size(); ++ i)
      {
        mean += vect[i];
        if (vect[i] > max)
          max = vect[i];
        if (vect[i] < min)
          min = vect[i];
      }
    //    std::cerr << name() << " Min/max = " << min << " / " << max << std::endl;
    mean /= vect.size();

  }
  /// \endcond
};


#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassification

  \brief Handle to a classification `Attribute_base`.

  \cgalModels Handle
*/
class Attribute_handle { };
#else
typedef boost::shared_ptr<Attribute_base> Attribute_handle;
#endif


} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_BASE_H
