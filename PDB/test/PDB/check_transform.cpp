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
#include <CGAL/PDB/range.h>
#include <fstream>
#include <cassert>
#include <CGAL/PDB/align.h>
#include <CGAL/PDB/distance.h>
#include <CGAL/PDB/Transform.h>

int main(int argc, char *argv[]){
  using namespace CGAL::PDB;
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  std::ifstream in("data/check_refine_alignment_0.pdb");
  PDB p(in);
  
  std::cout << "There are " << p.models().size() << " models." << std::endl;
  
  CGAL_PDB_FOREACH(PDB::Model_pair &m, p.models()) {
    std::cout << "Model " << m.key() << " has " << m.model().chains().size() << " chains" << std::endl;
    CGAL_PDB_FOREACH(Chain &c, make_chain_range(m.model().chains())) {
      /*Transform tr(0.309976, 0.851651, -0.422618,4,
                     -0.741526, 0.494764, 0.453154,5,
                     0.595025, 0.172916, 0.78488,6);*/
      Transform tr(1,0,0,1,
                   0,1,0,0,
                   0,0,1,0);
      write(tr);
      Point pt(1,2,3);
      std::cout << pt << std::endl;
      std::cout << tr(pt) << std::endl;
      
      Chain mc= c;
    
      //CGAL_PDB_FOREACH(Monomer &r, make_monomer_range(mc.monomers())) {
      CGAL_PDB_FOREACH(Atom &a, make_atom_range(mc.atoms())) {
        a.set_point(tr(a.point()));
      }
      //}
      double di= CGAL::squared_distance(mc.atoms().begin()->atom().point(),
                                        c.atoms().begin()->atom().point());
      std::cout<< "moved distance is " << di << std::endl;
      double di2= CGAL::squared_distance(*make_point_range(make_atom_range(mc.atoms())).begin(),
                                         *make_point_range(make_atom_range(c.atoms())).begin());
      std::cout<< "moved distance2 is " << di2 << std::endl;
      

        double err= cRMS(make_point_range(make_atom_range(c.atoms())),
                         make_point_range(make_atom_range(mc.atoms())));
        std::cout << "Err: " << err << std::endl;

        Transform trp= transform_taking_first_to_second( make_point_range(make_atom_range(c.atoms())), make_point_range(make_atom_range(mc.atoms())));

        double mdiff=0;
        for (unsigned int i=0; i< 4; ++i) {
          for (unsigned int j=0; j< 3; ++j) {
            mdiff+= std::abs(tr.cartesian(i,j)-trp.cartesian(i,j));
          }
        }
        CGAL_assertion(mdiff<.1);
        write(trp);
 
        err= cRMS(make_point_range(make_atom_range(c.atoms())),
                  make_point_range(make_atom_range(mc.atoms())));
        std::cout << "Err: " << err << std::endl;

        for_each(make_atom_range(c.atoms()), Transform_atom(trp));

        err= cRMS(make_point_range(make_atom_range(c.atoms())),
                  make_point_range(make_atom_range(mc.atoms())));
        std::cout << "Final err: " << err << std::endl;
        CGAL_assertion(err < 1e-5);
    }
  }

  
  return EXIT_SUCCESS;
}
