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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_PDB_UTILS_H
#define DSR_PDB_UTILS_H

namespace dsrpdb_internal {
  typedef enum Line_type {HEADER, DBREF, SEQRES, ATOM, HETATM, MASTER, ENDMDL, OTHER, TER, MODEL, END} Line_type;
  
  Line_type line_type (const char* line);
  extern const char atom_line_iformat_[];
}

#endif
