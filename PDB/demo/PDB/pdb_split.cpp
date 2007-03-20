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

#include <CGAL/PDB/PDB.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cctype>
#include <boost/program_options.hpp>

/*!
  \example pdb_split.cc

  This example shows how to split a pdb file various ways. A pdb can
  be split into different models, chains or event cut a particular
  chain.
*/

int main(int argc, char *argv[]){
  bool split_chains=false;
  std::string input_file, output_template;
  int split_domain=-1;
  bool print_help=false;
  bool verbose=false;
  char split_chain='\0';

  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out verbose messages about reading and writing pdb files")
    ("split-chains,c", boost::program_options::bool_switch(&split_chains),
     "Split all chains into separate files.")
    ("select-chain,C", boost::program_options::value<char>(&split_chain),
      "Select this chain only.")
    ("domain-split,s", boost::program_options::value<int>(&split_domain),
     "Split a chain into domains at this residue.");
  po.add_options()
    ("input-pdb", boost::program_options::value< std::string>(&input_file),
     "input file")
    ("output-pdb-template", boost::program_options::value< std::string>(&output_template),
     "A sprintf style string that will be used to generate the names for the output files.");

  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input-pdb", 1);
  p.add("output-pdb-template", 2);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm);



  //boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  //boost::program_options::notify(vm);

  if (input_file.empty()
      || output_template.empty()
      || (split_chains && split_chain != '\0')
      || (split_chains && split_domain != -1)
      || (split_chain != '\0' && split_domain != -1)
      || print_help) {
    std::cout << "This program splits a pdb file with multiple models or multiple domains into multiple files each with one model .\n";
    std::cout << "useage: " << argv[0]
	      << " [-c] input-pdb template%d[%c].pdb\n" << std::endl;
    std::cout << "The second argument is an sprintf style string that will be used to generate the names for the output files.\n\n";
    std::cout << o << "\n";
    return EXIT_SUCCESS;
  }



  if (split_chain != '\0'){
    std::cout << "Splitting on chain " << split_chain << std::endl;
  } else if (split_domain != -1) {
    std::cout << "Splitting on residue " << split_domain << std::endl;
  } else {
    std::cout << "Splitting into chains " << std::endl;
  }


  // std::cout << input_file << " " << output_template << " " << split_domain << " " << split_chains << std::endl;

  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  //= new char[strlen(argv[2]+1000)];
  CGAL_PDB_NS::PDB pdb(in, verbose);

  std::cout << "Input PDB has " << pdb.number_of_models() << " models." << std::endl;

  if (split_domain != -1 && pdb.number_of_models()!= 1){
    std::cerr << "Splitting a domain can only work if there is only one model and chain.\n";
    return EXIT_FAILURE;
  }

  for (unsigned int i=0; i< pdb.number_of_models(); ++i){
    const CGAL_PDB_NS::Model &m= pdb.model(i);
    if (split_chain != '\0' || split_chains || split_domain!= -1) {
      std::cout << "Model " << i << " has " << m.number_of_chains() << " chains."<< std::endl;
      if (split_domain != -1 && m.number_of_chains()!= 1){
	std::cerr << "Splitting a domain can only work if there is only one model and chain.\n";
	return EXIT_FAILURE;
      }
      for (unsigned int j=0; j< m.number_of_chains(); ++j){
	const CGAL_PDB_NS::Protein &p= m.chain(j);

	if (split_chain != '\0') {
	  if (p.chain() == split_chain || p.chain() == std::toupper(split_chain) ) {
	    std::cout << "Writing chain " << p.chain() << std::endl;
	    assert(split_domain ==-1);

	    std::ofstream out(output_template.c_str());
	    CGAL_PDB_NS::PDB npdb;
	    CGAL_PDB_NS::Model nm;
	    nm.new_chain(p);
	    npdb.new_model(nm);
	    npdb.set_header(pdb.header_begin(), pdb.header_end());
	    npdb.write(out);
	  } else {
	    std::cout << "Skipping chain " << p.chain() << std::endl;
	  }
	} else if (split_chains) {
	  std::cout << "Writing chain " << p.chain() << std::endl;
	  assert(split_domain ==-1);
	  char buf[100000];
	  if (pdb.number_of_models()==1) {
	    sprintf(buf, output_template.c_str(),p.chain());
	  } else {
	    sprintf(buf, output_template.c_str(), m.index(), p.chain());
	  }
	  std::ofstream out(buf);
	  CGAL_PDB_NS::PDB npdb;
	  CGAL_PDB_NS::Model nm;
	  nm.new_chain(p);
	  npdb.new_model(nm);
	  npdb.set_header(pdb.header_begin(), pdb.header_end());
	  npdb.write(out);
	} else {
	  CGAL_PDB_NS::Protein ps[2]={CGAL_PDB_NS::Protein(), CGAL_PDB_NS::Protein()};
	  for (CGAL_PDB_NS::Protein::Const_residues_iterator rit = p.residues_begin();
	       rit != p.residues_end(); ++rit){
	    //rit->write('c', std::cout);
	    if (rit->index().to_index()
		< static_cast<unsigned int>(split_domain)){
	      ps[0].new_residue(*rit);
	    } else {
	      ps[1].new_residue(*rit);
	    }
	  }

	  assert(ps[0].number_of_residues() + ps[1].number_of_residues()
		 == p.number_of_residues());

	  for (unsigned int j=0; j< 2; ++j){
	    std::cout << "Writing pdb with " << ps[j].number_of_residues() << " residues.\n";
	    char buf[100000];
	    sprintf(buf, output_template.c_str(), j);
	    std::ofstream out(buf);
	      CGAL_PDB_NS::PDB npdb;
	      CGAL_PDB_NS::Model nm(0);
	      nm.new_chain(ps[j]);
	      npdb.new_model(nm);
	      npdb.set_header(pdb.header_begin(), pdb.header_end());
	      npdb.write(out);
	    }
	  }
      }
    } else {
      char buf[100000];
      sprintf(buf, argv[2], m.index());
      std::ofstream out(buf);
      CGAL_PDB_NS::PDB npdb;
      npdb.new_model(m);
      npdb.set_header(pdb.header_begin(), pdb.header_end());
      npdb.write(out);
    }
  }
  //delete[] buf;
  return EXIT_SUCCESS;
}
