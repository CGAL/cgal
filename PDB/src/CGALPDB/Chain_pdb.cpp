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

#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/internal/pdb_utils.h>
#include <CGAL/PDB/internal/Error_logger.h>
#include <sstream>
CGAL_PDB_BEGIN_NAMESPACE

 

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

void Chain::process_line(const char *line) {
  CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
  if (lt== CGAL_PDB_INTERNAL_NS::TER) {

  } else if (lt == CGAL_PDB_INTERNAL_NS::ATOM) {
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
    int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::atom_line_iformat_,
			//"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
			&snum, name, &alt, resname, &chain, &resnum,
			&insertion_residue_code,
			&x,&y,&z, &occupancy, &tempFactor, segID, element, charge);
    /*if (chain_==' ') chain_=chain;
      if (!(chain_==' ' || chain== chain_)){
      std::ostringstream oss;
      oss << "Confusion over chain numbers. Expected " << chain_ << " got " << chain 
      << " on line:\n" << line;
      CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
      return;
      }*/

    if (resnum < 0) {
      std::ostringstream oss;
      oss << "Got negative residue index on line " << line;
      CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
      return;
    }

    Monomer_key resindex(resnum);
      
    Monomer::Atom_key al= Monomer::atom_key(name);
      
    if (al != Monomer::AL_OTHER) {
      Monomer *cur_residue=NULL;
  
      if (insertion_residue_code != ' '){
	if (insert_residues_.find(resindex) == insert_residues_.end()
	    || insert_residues_[resindex].find(IR_key(insertion_residue_code)) 
	    == insert_residues_[resindex].end()) {
	  Monomer::Type rl= Monomer::type(resname);
	  insert_residues_[resindex][IR_key(insertion_residue_code)]=Monomer(rl); 
	} 
	cur_residue= &insert_residues_[resindex][IR_key(insertion_residue_code)];
      } else {
	if (residues_.find(resindex) == residues_.end()) {
	  std::string nm(line, 17, 3);
	  Monomer::Type rl= Monomer::type(resname);
	  Monomer m(rl);
	  Monomer_iterator mit= residues_.insert(Container::value_type(resindex, m));
	  cur_residue= &mit->monomer();
	  CGAL_assertion(mit->key()== resindex);
	  //residues_[resindex]=m;
	  CGAL_postcondition(m.type()== cur_residue->type());
	} else {
	  CGAL_assertion(residues_.find(resindex) < residues_.end());
	  CGAL_assertion(residues_.find(resindex) >= residues_.begin());
	  cur_residue =&residues_.find(resindex)->monomer();
	  CGAL_assertion(residues_.find(resindex)->key()== resindex);
	}
      }    
      
      Atom a;
      //a.set_label(Residue::atom_label(al));
      a.set_point(Point(x,y,z));

      if (numscan >10) {
	a.set_occupancy(occupancy);
      }
      if (numscan >11) {
	a.set_temperature_factor(tempFactor);
      }
      a.set_segment_id(segID);
      a.set_element(element);
      a.set_charge(charge);
      
      /*if (cur_residue->index().to_index() 
	!= static_cast<unsigned int>(resnum)){
	std::ostringstream oss;
	oss << "Confusion over residue numbers. Expected" << cur_residue->index()
	<< " got " << resnum 
	<< " on line:\n" << line << std::endl;
	CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error(oss.str().c_str());
	return;
	}*/
      if (cur_residue->type() != Monomer::type(resname)){
	std::ostringstream oss;
	oss << "Confusion over residue types. Expected " 
	    << Monomer::type_string(cur_residue->type())
	    << " got " << Monomer::type_string(Monomer::type(resname))
	    << " on line:\n" << line << std::endl;
	CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error(oss.str().c_str());
      }
      //assert(cur_residue->type() == Residue::type(resname));
      Monomer::Atom_iterator rit= cur_residue->insert(al, a);
      if (rit == cur_residue->atoms_end()) {
	std::ostringstream oss;
	oss << "Error adding atom to residue " <<  resindex << std::endl;
	CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
      }
	
      //residue(resnum-1)->set_coords (al, Point(x,y,z));
      //residue(resnum-1)->set_index(al, snum);
      //++cur_atom_index;
    } else {
      /*std::ostringstream out;
	out << "Unhandled atom type: " << name;
	CGAL_PDB_INTERNAL_NS::error_logger.new_warning(out.str().c_str());*/
    }
  } else {
    assert(0);
  }
}


Chain::Chain(std::istream &in, bool print_errors) {
  char line[1000];
  
  CGAL_PDB_INTERNAL_NS::error_logger.set_is_output(print_errors);

  do {
    in.getline (line, 1000);
    if (!in) break;

    CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
    if (lt== CGAL_PDB_INTERNAL_NS::HEADER
	|| lt== CGAL_PDB_INTERNAL_NS::DBREF
	|| lt== CGAL_PDB_INTERNAL_NS::SEQRES
	|| lt == CGAL_PDB_INTERNAL_NS::COMPND){
      header_.push_back(std::string(line));
    } else if (lt== CGAL_PDB_INTERNAL_NS::ATOM){
      process_line(line);
    } else if (lt== CGAL_PDB_INTERNAL_NS::TER 
	       || lt == CGAL_PDB_INTERNAL_NS::END){
      break;
    } else {
      std::cout << "Skipping line: " << line << std::endl;
    }
  } while (true);
  set_has_bonds(true);
  CGAL_PDB_INTERNAL_NS::error_logger.dump();
}
 
void Chain::write_pdb(std::ostream &out) const {
  assert(!residues_.empty());
  for (unsigned int i=0; i< header_.size(); ++i){
    out << header_[i] << std::endl;
  }

  char line[81];
  
  sprintf(line, "MODEL %8d         ", 1);
  out << line << std::endl;
  
  //    int anum=1;
  write(' ' , 1, out);
  
  out << "ENDMDL                       " << std::endl;
}

int Chain::write(char chain, int start_index, std::ostream &out) const {
  char line[81];
  //    int anum=1;
  Monomer_key last_resindex;
  Monomer::Type last_type= Monomer::INV;
  for (Monomer_const_iterator it = monomers_begin(); it != monomers_end(); ++it) {
    const Monomer &res= it->monomer();
    //Residue::Label rl =  res.label();
    //residues_[i]->atoms();
    start_index= res.write(chain, it->key().to_index(), ' ', start_index, out);
    
    IR_Map::const_iterator irit= insert_residues_.find(it->key());
    if (irit!= insert_residues_.end()) {
      for (unsigned int i=0; i< irit->data().size(); ++i){
	start_index= 
	  irit->data().find(IR_key(i))->data().write(chain, it->key().to_index(),
						     irit->data().find(IR_key(i))->key().to_index(), start_index, out);
      }
    }
    last_resindex= it->key();
    last_type= it->data().type();
  }
  const char *terformat="TER   %5d      %3s %c%3d%c";
  if (!residues_.empty()) {
    sprintf(line, terformat, start_index, 
	    Monomer::type_string(last_type).c_str(), chain, 
	    last_resindex.to_index(),' ');
    out << line << std::endl;
  }
  return start_index+1;
}

CGAL_PDB_END_NAMESPACE
