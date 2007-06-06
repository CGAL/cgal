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


#ifndef CGAL_DSRPDB_RMS_H
#define CGAL_DSRPDB_RMS_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/Matrix.h>
#include <cmath>
#include <vector>
CGAL_PDB_BEGIN_NAMESPACE
//! Compute the all-atom cRMS with no alignment. 
/*!  This function assumes that a, and b have already been rigidly
  aligned externally. 
    
  \pre {a and b must have the same number of atoms. This is checked.}

  See pdb_distance.cc for an example.
*/
double no_align_cRMS(const Protein &a, const Protein &b);


//! Compute the calpha cRMS without aligning the proteins.
/*!
  \pre {a and b must have the same number of atoms. This is checked.}
    
  See pdb_distance.cc for an example.
*/
double no_align_ca_cRMS(const Protein &a, const Protein &b);


//! Compute the all-atom cRMS between two proteins
/*!
  \pre {a and b must have the same number of atoms. This is checked.}
*/
double cRMS(const Protein &a, const Protein &b);


//! Compute the C_alpha cRMS between two proteins
/*!
  \pre {a and b must have the same number of atoms. This is checked.}
*/
double ca_cRMS(const Protein &a, const Protein &b);


//! Compute the all-atom dRMS between two proteins
/*!
  \pre {a and b must have the same number of atoms. This is checked.}

  See pdb_distance.cc for an example.
*/
double dRMS(const Protein &a, const Protein &b);


//! Compute the C_alpha cRMS between two proteins
/*!
  \pre {a and b must have the same number of atoms. This is checked.}

  See pdb_distance.cc for an example.
*/
double ca_dRMS(const Protein &a, const Protein &b);


//! Return the distance matrix of a protein
/*!
  See pdb_distance_matrix.cc for an example.
*/
Matrix distance_matrix(const Protein &a);


//! Return the C_alpha distance matrix of a protein
Matrix ca_distance_matrix(const Protein &a);

//! Return the backbone atom distance matrix of a protein
Matrix backbone_distance_matrix(const Protein &a);


//! Return the distance matrix of a protein
Matrix distance_matrix(const Model &a);


//! Return the C_alpha distance matrix of a protein
Matrix ca_distance_matrix(const Model &a);

//! Return the backbone atom distance matrix of a protein
Matrix backbone_distance_matrix(const Model &a);

CGAL_PDB_END_NAMESPACE
#endif
