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

#ifndef CGAL_POLYHEDRAL_CGM_IOSTREAM_H
#define CGAL_POLYHEDRAL_CGM_IOSTREAM_H

/*! \file Implmentation of importer and exporter IO operations of a
 * Polyhedral_cgm object
 */

#include <CGAL/basic.h>
#include <CGAL/Polyhedral_cgm.h>
#include <CGAL/IO/Cgm_iostream.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE

/*! Exporter operator of a Polyhedral_cgm object
 * \param os the output stream
 * \param cgm the Polyhedral_cgm object
 */
template <class Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline std::ostream & operator<<(std::ostream & os,
                                 const Polyhedral_cgm<Kernel,T_Dcel> & cgm) 
{
  const Cubical_gaussian_map_3<Kernel,T_Dcel> * tmp = &cgm;
  os << *tmp;
  return os;
}

/*! Importer operator of a Polyhedral_cgm object
 * \param is the input stream
 * \param cgm the Polyhedral_cgm object
 */
template <class Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline std::istream & operator>>(std::istream & is,
                                 Polyhedral_cgm<Kernel,T_Dcel> & cgm)
{
  static_cast<Cubical_gaussian_map_3<Kernel,T_Dcel>(cgm) >> is;
  return is;
}

CGAL_END_NAMESPACE

#endif
