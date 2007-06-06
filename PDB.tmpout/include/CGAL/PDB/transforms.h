/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#ifndef CGAL_DSR_PDB_TRANSFORMS_H
#define CGAL_DSR_PDB_TRANSFORMS_H
#include <CGAL/PDB/basic.h>
CGAL_PDB_BEGIN_NAMESPACE


  class Protein;
  class Transform;
  class PDB;

  //! Apply a Transform to all atoms of a Protein file
  void transform_protein(const Transform &t, Protein &p);
 

  //! Apply a Transform to all atoms of a PDB file
  void transform_pdb(const Transform &t, PDB &pdb);

CGAL_PDB_END_NAMESPACE

#endif
