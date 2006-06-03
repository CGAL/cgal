// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_CGM_IOSTREAM_H
#define CGAL_CGM_IOSTREAM_H

/*! \file Implmentation of importer and exporter IO operations of a
 * Cubical_gaussian_map_3 object
 */

#include <CGAL/basic.h>
#include <CGAL/Cubical_gaussian_map_3.h>
#include <CGAL/IO/Arr_iostream.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE

/*! Exporter operator of a Projected_normal object
 * \param os the output stream
 * \param proj_normal the Projected_normal object
 */
template <class Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline std::ostream &
operator<<(std::ostream & os,
           const typename Cubical_gaussian_map_3<Kernel,T_Dcel>::Projected_normal &
           proj_normal)
{
  return os << proj_normal.get_projected_normal() << ", "
            << std::hex << proj_normal.get_faces_mask() << ", "
            << proj_normal.get_num_faces();
}

/*! Exporter operator of a Cubical_gaussian_map_3 object
 * \param os the output stream
 * \param cgm the Cubical_gaussian_map_3 object
 */
template <class Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline
std::ostream & operator<<(std::ostream & os,
                          const Cubical_gaussian_map_3<Kernel,T_Dcel> & cgm) 
{
  typedef Cubical_gaussian_map_3<Kernel,T_Dcel> Cgm;
  for (unsigned int i = 0; i < Cgm::NUM_FACES; ++i) {
    const typename Cgm::Arrangement & arr = cgm.arrangement(i);
    os << arr;
  }
  return os;
}

/*! Importer operator of a Cubical_gaussian_map_3 object
 * \param is the input stream
 * \param cgm the Cubical_gaussian_map_3 object
 */
template <class Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline std::istream & operator>>(std::istream & is,
                                 Cubical_gaussian_map_3<Kernel,T_Dcel> & cgm)
{
  typedef Cubical_gaussian_map_3<Kernel,T_Dcel> Cgm;
  for (unsigned int i = 0; i < Cubical_gaussian_map_3::NUM_FACES; ++i) {
    const typename Cgm::Arrangement & arr = cgm.get_arrangement(i);
    arr >> is;
  }
  return is;
}

CGAL_END_NAMESPACE

#endif

