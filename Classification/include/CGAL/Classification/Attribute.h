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

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_H

#include <vector>

namespace CGAL {

namespace Classification {
  
/*!
  \ingroup PkgClassification

  \brief Abstract class describing a classification attribute.

  A classification attribute associates a scalar value to each item
  of the classification input.
*/


class Attribute
{
public:
  double weight; ///< Weight of the attribute between 0 and +&infin;.

  
  /// \cond SKIP_IN_MANUAL
  double mean;
  double max;

  virtual ~Attribute() { }
  /// \endcond

  /*!
    \brief Returns the ID of the attribute.

    This method returns _abstract\_attribute_ and should be overloaded
    by a child class with a specific ID.
  */
  virtual std::string id() { return "abstract_attribute"; }

  /*!
    \brief Returns the value taken by the attribute on the given item.

    This method must be implemented by inherited classes.

    \param index Index of the input item.

    \return the value of the attribute on the item of index `index`
  */
  virtual double value (std::size_t index) = 0;

  /// \cond SKIP_IN_MANUAL
  virtual double normalized (std::size_t index)
  {
    return (std::max) (0., (std::min) (1., value(index) / weight));
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
    //    std::cerr << id() << " Min/max = " << min << " / " << max << std::endl;
    mean /= vect.size();

  }
  /// \endcond
};


#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassification

  \brief Handle to a classification `Attribute`.

  \cgalModels Handle
*/
class Attribute_handle { };
#else
typedef boost::shared_ptr<Attribute> Attribute_handle;
#endif
  
} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_ATTRIBUTE_H
