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

#include <dsrpdb/Protein.h>
#include <dsrpdb/Residue.h>
#include "pdb_utils.h"
#include <dsrpdb_internal/Error_logger.h>
#include <sstream>

namespace dsrpdb {


 

  /*static unsigned int getSequenceBeginning (const char* line)
  {

    char s[5];

    strncpy (s, line + 14, 4); s[4] = '\0';
    return atoi (s);
    }*/


  /*static  unsigned int getLineSerialNumber (const char* line) {
    char s[3];
    strncpy (s, line + 8, 2); s[2] = '\0';
    return atoi (s);
    }*/
  
  /*static unsigned int getNumberOfResidues (const char* line)
  {

    char s[5];

    strncpy (s, line + 13, 4); s[4] = '\0';
    return atoi (s);

    }*/
  
 unsigned int  getSequenceNumber (const char* line)
  {

    char s[5];

    strncpy (s, line + 22, 4); s[4] = '\0';
    return atoi (s);

  }

  void Protein::process_line(const char *line) {
    dsrpdb_internal::Line_type lt= dsrpdb_internal::line_type(line);
    if (lt== dsrpdb_internal::TER) {

    } else if (lt == dsrpdb_internal::ATOM) {
      // the read values are not zero padded so we must fill the buffers for strings with 0s.
      int snum=-1;
      char name[5]={'\0','\0','\0','\0','\0',};
      char alt='\0';
      char resname[4]={'\0','\0','\0','\0'};
      char chain;
      int resnum=-1;
      char insertion_residue_code;
      float x,y,z;
      float occupancy, tempFactor;
      char segID[5]={'\0','\0','\0','\0','\0'};
      char element[3]={'\0','\0','\0'};
      char charge[3]={'\0','\0','\0'};
      int numscan= sscanf(line, dsrpdb_internal::atom_line_iformat_,
			  //"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
			  &snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
			  &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);
      if (chain_==' ') chain_=chain;
      if (!(chain_==' ' || chain== chain_)){
	std::ostringstream oss;
	oss << "Confusion over chain numbers. Expected " << chain_ << " got " << chain 
	    << " on line:\n" << line;
	dsrpdb_internal::error_logger.new_warning(oss.str().c_str());
	return;
      }

      if (resnum < 0) {
	std::ostringstream oss;
	oss << "Got negative residue index on line " << line;
	dsrpdb_internal::error_logger.new_warning(oss.str().c_str());
	return;
      }

      Residue::Index resindex(resnum);

      if (insertion_residue_code != ' ' 
	  && (!residues_.empty() && resindex == residues_.back().index())){
	static bool skip_alt=false;
	if (!skip_alt){
	  dsrpdb_internal::error_logger.new_warning("Skipping alternate residue.");
	}
	skip_alt=true;
	return;
      }
      
      
      Residue::Atom_label al= Residue::atom_label(name);
      
      if (al != Residue::AL_OTHER) {
	if (residues_.empty() || residues_.back().index() != resindex) {
	  std::string nm(line, 17, 3);
	  Residue::Type rl= Residue::type(resname);
	  residues_.push_back(Residue(rl));
	  residues_.back().set_index(resindex);
	}

	Atom::Index sindex(snum);
	

	Atom a;
	//a.set_label(Residue::atom_label(al));
	a.set_cartesian_coords(Point(x,y,z));
	a.set_index(sindex);
	if (numscan >10) {
	  a.set_occupancy(occupancy);
	}
	if (numscan >11) {
	  a.set_temperature_factor(tempFactor);
	}
	a.set_segment_id(segID);
	a.set_element(element);
	a.set_charge(charge);
      
	if (static_cast<unsigned int>(residues_.back().index()) 
	    != static_cast<unsigned int>(resnum)){
	  std::ostringstream oss;
	  oss << "Confusion over residue numbers. Expected" << residues_.back().index()
	      << " got " << resnum 
	      << " on line:\n" << line << std::endl;
	  dsrpdb_internal::error_logger.new_fatal_error(oss.str().c_str());
	  return;
	}
	if (residues_.back().type() != Residue::type(resname)){
	  std::ostringstream oss;
	  oss << "Confusion over residue types. Expected" 
	      << Residue::type_string(residues_.back().type())
	      << " got " << Residue::type_string(Residue::type(resname))
	      << " on line:\n" << line << std::endl;
	  dsrpdb_internal::error_logger.new_fatal_error(oss.str().c_str());
	}
	//assert(residues_.back().type() == Residue::type(resname));
	
	residues_.back().set_atom(al, a);
	
	//residue(resnum-1)->set_coords (al, Point(x,y,z));
	//residue(resnum-1)->set_index(al, snum);
	//++cur_atom_index;
      } else {
	/*std::ostringstream out;
	out << "Unhandled atom type: " << name;
	dsrpdb_internal::error_logger.new_warning(out.str().c_str());*/
      }
    } else {
      assert(0);
    }
  }


  Protein::Protein(std::istream &in, bool print_errors) {
    char line[1000];
    chain_=' ';
  
    dsrpdb_internal::error_logger.set_is_output(print_errors);

    do {
      in.getline (line, 1000);
      if (!in) break;

      dsrpdb_internal::Line_type lt= dsrpdb_internal::line_type(line);
      if (lt== dsrpdb_internal::HEADER) {
	header_.push_back(std::string(line));
      } else if (lt== dsrpdb_internal::DBREF){
	header_.push_back(std::string(line));
	//seqBegin = getSequenceBeginning (line);
      } else if (lt== dsrpdb_internal::SEQRES){
	header_.push_back(std::string(line));
      } else if (lt== dsrpdb_internal::ATOM){
	process_line(line);
      } else if (lt== dsrpdb_internal::TER || lt == dsrpdb_internal::END){
	break;
      }
    } while (true);
    set_has_bonds(true);
    dsrpdb_internal::error_logger.dump();
  }
 
  void Protein::write_pdb(std::ostream &out) const {
    assert(!residues_.empty());
    for (unsigned int i=0; i< header_.size(); ++i){
      out << header_[i] << std::endl;
    }

    char line[81];
  
    sprintf(line, "MODEL %8d         ", 1);
    out << line << std::endl;
  
    //    int anum=1;
    write(out);
  
    out << "ENDMDL                       " << std::endl;
  }

  void Protein::write(std::ostream &out) const {
    char line[81];
    //    int anum=1;
 
    for (unsigned int i = 0; i < residues_.size(); i++) {
      const Residue &res= residues_[i];
      //Residue::Label rl =  res.label();
      //residues_[i]->atoms();
      res.write(chain_, out);
    }
    const char *terformat="TER   %5d      %3s %c%3d%c";
    if (!residues_.empty()) {
      sprintf(line, terformat, static_cast<unsigned int>(residues_.back().last_atom_index())+1, 
	      Residue::type_string(residues_.back().type()).c_str(), chain(), static_cast<unsigned int>(residues_.back().index()),' ');
      out << line << std::endl;
    }
  }
};
