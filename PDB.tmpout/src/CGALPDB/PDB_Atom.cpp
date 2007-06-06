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

#include <CGAL/PDB/Atom.h>
#include <CGAL/PDB/internal/Error_logger.h>
#include <cstring>
#include <sstream>
#include <limits>
CGAL_PDB_BEGIN_NAMESPACE

static double nan__=std::numeric_limits<double>::signaling_NaN();

//enum Type {INVALID=0, C,N,H,O, S, FE, PT};
// Taken from the wikipedia, so don't trust them too much.
double Atom::radii_[]={0, 1.70, 1.55, 1.20, 1.52, 1.85, 1.10, 1.75, nan__, nan__, nan__, nan__};
  
Atom::Type Atom::string_to_type(const char *cp) {
  if (std::strchr(cp, 'C') != NULL) return C;
  else if (std::strchr(cp,'N') != NULL) return N;
  else if (std::strchr(cp,'S') != NULL) return S;
  else if (std::strchr(cp,'O') != NULL) return O;
  else if (std::strchr(cp,'H') != NULL) return H;
  else if (std::strstr(cp,"FE") != NULL) return FE; 
  else if (std::strstr(cp,"PT") != NULL) return PT;  
  else {
    std::ostringstream em;
    em << "Couldn't parse atom type of " << cp;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(em.str().c_str());
    return INVALID;
  }
}
CGAL_PDB_END_NAMESPACE
