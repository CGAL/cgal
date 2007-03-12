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

#include <CGAL/PDB/PDB.h>
#include <fstream>
#include <cassert>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/cgal.h>

#include "check_equal.h"


int main(int , char *[]){
	//dsr::Residueres= dsr::Residue(dsr::Residue::VAL);
	//res.write(std::cout);
	std::string	argv1="data/check_pdb.pdb";
	std::ifstream in(argv1.c_str());
	CGAL_PDB_NS::PDB p(in);
	//p.write(std::cout);
	std::ostringstream of;
	
	std::cout << "There are " << p.number_of_models() << " models." << std::endl;
	
	for (unsigned int i=0; i< p.number_of_models(); ++i){
		const CGAL_PDB_NS::Model &m= p.model(i);
		std::cout << "Model " << i << " has " << m.number_of_chains() << " chains" << std::endl;
		for (unsigned int j=0; j< m.number_of_chains(); ++j){
			std::cout << "Chain " << j << " has " << m.chain(j).number_of_residues() << " residues" << std::endl;
		}
	}
	
	p.write(of);
	//of << std::fflush;
	
#ifdef PDB_USE_CGAL
	if (false){
		typedef CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> Tr;
		
		std::vector<Tr::Weighted_point> wps;
		CGAL_PDB_NS::all_weighted_points(p.model(0), Tr(), 
										 std::back_inserter(wps));
		std::copy(wps.begin(), wps.end(),
				  std::ostream_iterator<Tr::Weighted_point>(std::cout, "\n"));
	}
#endif
	std::ifstream nif(argv1.c_str());
	std::istringstream ofis(of.str().c_str());
	check_equal(nif, ofis);
	
	return return_code__;
}
