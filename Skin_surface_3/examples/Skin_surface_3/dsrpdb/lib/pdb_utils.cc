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

#include "pdb_utils.h"

#include <string>
#include <iostream>
#include <set>
#include <dsrpdb_internal/Error_logger.h>

namespace dsrpdb_internal {


  //ATOM    812  OG ASER   106       -.072  22.447  10.384   .50 11.73      1ECD 918
  //ATOM     60  O   SER L   9     -26.231  10.210  -4.537  1.00 26.25      7FAB 160

  const char atom_line_iformat_[]=
  "ATOM  %5d%*1c%4c%1c%3c%*c%1c%4d%1c%*1c%*1c%*1c%8f%8f%8f%6f%6f%*1c%*1c%*1c%*1c%*1c%*1c%4c%2c%2c";

  Line_type line_type (const char* line) {
    if (line[0]=='\0') return OTHER;
    //char s[7];
    //strncpy (s, line, 6); s[6] = '\0';
    std::string word(line,0,6);
    if (word== "DBREF ") return DBREF;
    else if (word=="SEQRES") return SEQRES;
    else if (word== "ATOM  ") return ATOM;
    else if (word== "HETATM") return HETATM;
    else if (word== "MASTER") return MASTER;
    else if (word== "ENDMDL") return ENDMDL;
    else if (word== "END   " || word== "END  " || word== "END " || word== "END") return END;
    else if (   word == "HEADER" || word == "TITLE " || word == "COMPND" || word == "SOURCE"
	      || word == "KEYWDS" || word == "EXPDTA" || word == "AUTHOR" || word == "REVDAT"
		|| word == "JRNL  " || word == "REMARK" || word == "HELIX " || word == "SHEET "
		|| word == "SITE  " || word == "CRYST1" || word == "ORIGX1" || word == "ORIGX2"
		|| word == "ORIGX3" || word == "SCALE1" || word == "SCALE2" || word == "SCALE3"
		|| word == "SEQADV" || word == "TURN  " || word == "FORMUL" || word == "HETNAM"
		|| word == "SSBOND" || word == "MODRES" || word == "CONECT" || word == "HET   "
		|| word == "CISPEP") return HEADER;
    else if (word== "TER   " || word == "TER") return TER;
    else if (word== "MODEL ") return MODEL;
    else {
      error_logger.new_warning(std::string("\"" + word + "\" is not a known line type.").c_str());
      return OTHER;
    }
  }
}
