// Copyright (c) 2017 GeometryFactory Sarl (France).
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

#ifndef CGAL_CLASSIFICATION_LABEL_H
#define CGAL_CLASSIFICATION_LABEL_H

#include <CGAL/Classification/Feature/Effect.h>

#include <boost/shared_ptr.hpp>

namespace CGAL {

namespace Classification {

/*!
\ingroup PkgClassification

\brief %Classification label (for example: vegetation, ground, etc.)
defined as a set of relationships with classification features.

*/
class Label
{
public:
  
private:
  /// \cond SKIP_IN_MANUAL
  std::string m_name;
  std::map<Feature_handle, Feature::Effect> m_feature_effects;
  /// \endcond

public:

  /*! 
    \param name name of the classification label (e.g. vegetation).
  */ 
  Label (std::string name) : m_name (name) { }

  /*! 
    \brief Sets the effect of feature `att` on the classification label.
  */ 
  void set_feature_effect (Feature_handle att, Feature::Effect effect)
  {
    m_feature_effects[att] = effect;
  }

  /*!
    \brief Returns the effect of feature `att` on the classification label.
   */
  Feature::Effect feature_effect (Feature_handle att) 
  {
    std::map<Feature_handle, Feature::Effect>::iterator
      search = m_feature_effects.find (att);
    return (search == m_feature_effects.end () ? Feature::NEUTRAL : search->second);
  }

  const std::string& name() const { return m_name; }
  
  /// \cond SKIP_IN_MANUAL
  void info()
  {
    std::cerr << "Feature " << m_name << ": ";
    for (std::map<Feature_handle, Feature::Effect>::iterator it = m_feature_effects.begin();
         it != m_feature_effects.end(); ++ it)
      {
        if (it->second == Feature::NEUTRAL)
          continue;
        
        std::cerr << it->first;
        if (it->second == Feature::FAVORING) std::cerr << " (favored), ";
        else if (it->second == Feature::PENALIZING) std::cerr << " (penalized), ";
      }
    std::cerr << std::endl;
  }
  /// \endcond

};

#ifdef DOXYGEN_RUNNING
/*!
  \ingroup PkgClassification

  \brief %Handle to a classification `Label`.

  \cgalModels Handle
*/
class Label_handle { };
#else
typedef boost::shared_ptr<Label> Label_handle;
#endif

} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_LABEL_H
