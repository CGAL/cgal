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

#include "dsrpdb/Residue.h"
#include "Residue_data.h"
#include <cassert>

namespace dsrpdb {
  namespace Residue_data {
    bool amino_acid_initialized_=false;
    typedef std::vector<std::vector<Residue::Atom_label> > Clean_atom_fallbacks;
    std::vector<Amino_acid_data > amino_acid_data_;



  
    //! The per atom data associating string names with labels and atom types
    /*!  To add a new atom or a different string refering to the same
      atom add a line to this array.
    */
    Atom_data atom_name_data_[]= {
      {" N  ",Residue::AL_N, Atom::N},
      {" H  ",Residue::AL_H, Atom::H},
      {"1H  ",Residue::AL_1H, Atom::H},{" H1 ",Residue::AL_1H, Atom::H},
      {"2H  ",Residue::AL_2H, Atom::H},{" H2 ",Residue::AL_2H, Atom::H},
      {"3H  ",Residue::AL_3H, Atom::H},{" H3 ",Residue::AL_3H, Atom::H},
      {" C  ",Residue::AL_C, Atom::C},
      {" O  ",Residue::AL_O, Atom::O},
      {" OXT",Residue::AL_OXT, Atom::O},
      {" CH3",Residue::AL_CH3, Atom::C},

      {" CA ",Residue::AL_CA, Atom::C},
      {" HA ",Residue::AL_HA, Atom::H},
      {"1HA ",Residue::AL_1HA, Atom::H},{" HA1",Residue::AL_1HA, Atom::H},
      {"2HA ",Residue::AL_2HA, Atom::H},{" HA2",Residue::AL_2HA, Atom::H},

      {" CB ",Residue::AL_CB, Atom::C},
      {" HB ",Residue::AL_HB, Atom::H},
      {"1HB ",Residue::AL_1HB, Atom::H}, {" HB1",Residue::AL_1HB, Atom::H}, 
      {"2HB ", Residue::AL_2HB, Atom::H}, {" HB2", Residue::AL_2HB, Atom::H},
      {"3HB ", Residue::AL_3HB, Atom::H}, {" HB3", Residue::AL_3HB, Atom::H},

      {" CG ",Residue::AL_CG, Atom::C},
      {" CG1",Residue::AL_CG1, Atom::C},
      {" CG2",Residue::AL_CG2, Atom::C},
      {" HG ",Residue::AL_HG, Atom::H},
      {"1HG ", Residue::AL_1HG, Atom::H}, {" HG1", Residue::AL_1HG, Atom::H},
      {"2HG ", Residue::AL_2HG, Atom::H}, {" HG2", Residue::AL_2HG, Atom::H},
      //{"HG1",Residue::AL_HG1},
      {"1HG1",Residue::AL_1HG1, Atom::H},{"HG11",Residue::AL_1HG1, Atom::H},
      {"2HG1",Residue::AL_2HG1, Atom::H},{"HG12",Residue::AL_2HG1, Atom::H},
      {"3HG1",Residue::AL_3HG1, Atom::H},{"HG13",Residue::AL_3HG1, Atom::H},
      {"1HG2",Residue::AL_1HG2, Atom::H},{"HG21",Residue::AL_1HG2, Atom::H},
      {"2HG2",Residue::AL_2HG2, Atom::H},{"HG22",Residue::AL_2HG2, Atom::H},
      {"3HG2",Residue::AL_3HG2, Atom::H},{"HG23",Residue::AL_3HG2, Atom::H},
      {" OG ",Residue::AL_OG, Atom::O},
      {" OG1",Residue::AL_OG1, Atom::O},
      {" SG ",Residue::AL_SG, Atom::S},

      {" CD ",Residue::AL_CD, Atom::C},
      {" CD1",Residue::AL_CD1, Atom::C},
      {" CD2",Residue::AL_CD2, Atom::C},
      //{"HD1",Residue::AL_HD1},
      //{"HD2",Residue::AL_HD2},
      {" HD ",Residue::AL_HD, Atom::H},
      {"1HD ",Residue::AL_1HD, Atom::H},{" HD1",Residue::AL_1HD, Atom::H},
      {"2HD ",Residue::AL_2HD, Atom::H},{" HD2",Residue::AL_2HD, Atom::H},
      {"3HD ",Residue::AL_3HD, Atom::H},{" HD3",Residue::AL_3HD, Atom::H},
      {"1HD1",Residue::AL_1HD1, Atom::H}, {"HD11",Residue::AL_1HD1, Atom::H}, 
      {"2HD1",Residue::AL_2HD1, Atom::H},{"HD12",Residue::AL_2HD1, Atom::H},
      {"3HD1",Residue::AL_3HD1, Atom::H},{"HD13",Residue::AL_3HD1, Atom::H},
      {"1HD2",Residue::AL_1HD2, Atom::H},{"HD21",Residue::AL_1HD2, Atom::H},
      {"2HD2",Residue::AL_2HD2, Atom::H},{"HD22",Residue::AL_2HD2, Atom::H},
      {"3HD2",Residue::AL_3HD2, Atom::H},{"HD23",Residue::AL_3HD2, Atom::H},
      {" SD ",Residue::AL_SD, Atom::S},
      {" OD1",Residue::AL_OD1, Atom::O},
      {" OD2",Residue::AL_OD2, Atom::O},
      {" ND1",Residue::AL_ND1, Atom::N},
      {" ND2",Residue::AL_ND2, Atom::N},

      {" CE ",Residue::AL_CE, Atom::C},
      {" CE1",Residue::AL_CE1, Atom::C},
      {" CE2",Residue::AL_CE2, Atom::C},
      {" CE3",Residue::AL_CE3, Atom::C},
      {" HE ",Residue::AL_HE, Atom::H},
      {"1HE ",Residue::AL_1HE, Atom::H},{" HE1",Residue::AL_1HE, Atom::H},
      {"2HE ",Residue::AL_2HE, Atom::H},{" HE2",Residue::AL_2HE, Atom::H},
      {"3HE ",Residue::AL_3HE, Atom::H},{" HE3",Residue::AL_3HE, Atom::H},
      //{"HE1",Residue::AL_HE1},
      //{"HE2",Residue::AL_HE2},
      //{"HE3",Residue::AL_HE3},
      {"1HE2",Residue::AL_1HE2, Atom::H},{"HE21",Residue::AL_1HE2, Atom::H},
      {"2HE2",Residue::AL_2HE2, Atom::H},{"HE22",Residue::AL_2HE2, Atom::H},
      {" OE1",Residue::AL_OE1, Atom::O},
      {" OE2",Residue::AL_OE2, Atom::O},
      {" NE ",Residue::AL_NE, Atom::N},
      {" NE1",Residue::AL_NE1, Atom::N},
      {" NE2",Residue::AL_NE2, Atom::N},

      {" CZ ",Residue::AL_CZ, Atom::C},
      {" CZ2",Residue::AL_CZ2, Atom::C},
      {" CZ3",Residue::AL_CZ3, Atom::C},
      {" NZ ",Residue::AL_NZ, Atom::N},
      {" HZ ",Residue::AL_HZ, Atom::H},
      {"1HZ ",Residue::AL_1HZ, Atom::H},{" HZ1",Residue::AL_1HZ, Atom::H},
      {"2HZ ",Residue::AL_2HZ, Atom::H},{" HZ2",Residue::AL_2HZ, Atom::H},
      {"3HZ ",Residue::AL_3HZ, Atom::H},{" HZ3",Residue::AL_3HZ, Atom::H},
      //{"HZ1",Residue::AL_HZ2},
      //{"HZ2",Residue::AL_HZ2},
      //{"HZ3",Residue::AL_HZ3},

      {" CH2",Residue::AL_CH2, Atom::C},
      {" NH1",Residue::AL_NH1, Atom::N},
      {" NH2",Residue::AL_NH2, Atom::N},
      {" OH",Residue::AL_OH, Atom::O},
      {" HH",Residue::AL_HH, Atom::H},
      
      {"1HH1", Residue::AL_1HH1, Atom::H}, {"HH11", Residue::AL_1HH1, Atom::H},
      {"2HH1", Residue::AL_2HH1, Atom::H}, {"HH12", Residue::AL_2HH1, Atom::H},
      {" HH2",Residue::AL_HH2, Atom::H},
      {"1HH2", Residue::AL_1HH2, Atom::H}, {"HH21", Residue::AL_1HH2, Atom::H},
      {"2HH2", Residue::AL_2HH2, Atom::H}, {"HH22", Residue::AL_2HH2, Atom::H},
      {" HH ", Residue::AL_1HH2, Atom::H},

      {"HH31", Residue::AL_1HH3, Atom::H}, 
      {"HH32", Residue::AL_2HH3, Atom::H}, 
      {"HH33", Residue::AL_3HH3, Atom::H},
      {"UNKN", Residue::AL_INVALID, Atom::INVALID}};
    

    Atom_data clean_atom_name_data_[sizeof(atom_name_data_)
				    /sizeof(Atom_data)+1];

    Atom_fallback_data atom_fallback_data_[]=
      {{Residue::AL_CD1, Residue::AL_CD}, 
       {Residue::AL_1HA, Residue::AL_HA},
       {Residue::AL_1HB, Residue::AL_HB},
       {Residue::AL_1HD, Residue::AL_HD},
       {Residue::AL_1HE, Residue::AL_HE},
       {Residue::AL_1HZ, Residue::AL_HZ},
       //{Residue::AL_1HH, Residue::AL_HH},
       {Residue::AL_1HD1, Residue::AL_1HD},
       {Residue::AL_2HD1, Residue::AL_2HD},
       {Residue::AL_3HD1, Residue::AL_3HD},
       {Residue::AL_INVALID, Residue::AL_INVALID}};
    Clean_atom_fallbacks clean_atom_fallbacks_;

#define BEGIN_RES(name) case name:{
#define DEFINE_ATOMS(...) static Residue::Atom_label cl[]={__VA_ARGS__}; dat.atms_=cl;
#define DEFINE_BONDS(...) static Residue::Atom_label cbl[]={__VA_ARGS__}; for (int i=0; cbl[i] != Residue::AL_INVALID; ++i) { bool found=false; for (int j=0; cl[j] != Residue::AL_INVALID; ++j) if (cl[j]==cbl[i]) found=true; assert(found || cbl[i]== Residue::AL_CA || cbl[i]== Residue::AL_N);}; dat.bnds_=cbl; 
#define DEFINE_EXTREMES(...) static Residue::Atom_label ext[]={__VA_ARGS__}; for (int i=0; ext[i] != Residue::AL_INVALID; ++i) { bool found=false; for (int j=0; cl[j] != Residue::AL_INVALID; ++j) if (cl[j]==ext[i]) found=true; assert(found);} dat.extr_=ext; 
#define END_RES break;}


    //! This function returns the per-residue atom and bond information.
    /*!
      To add an atom or bond to a residue add them to the appropriate list below.
    */
    Residue_init_data 
    get_residue_initialization_data(Residue::Type rl){
      static Residue::Atom_label end=Residue::AL_INVALID;
      Residue_init_data dat;
      dat.atms_=&end;
      dat.bnds_=&end;
      dat.extr_=&end;
      switch(rl){

	BEGIN_RES(Residue::VAL);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_HB, Residue::AL_CG1,  Residue::AL_CG2, 
		     Residue::AL_1HG1, Residue::AL_2HG1, Residue::AL_3HG1, Residue::AL_1HG2,
		     Residue::AL_2HG2, Residue::AL_3HG2, Residue::AL_INVALID );
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_HB, Residue::AL_CB, 
		     Residue::AL_CG1, Residue::AL_CB,
		     Residue::AL_CG2, Residue::AL_CB, 
		     Residue::AL_1HG1, Residue::AL_CG1, 
		     Residue::AL_2HG1, Residue::AL_CG1,
		     Residue::AL_3HG1, Residue::AL_CG1, 
		     Residue::AL_1HG2, Residue::AL_CG2, 
		     Residue::AL_2HG2, Residue::AL_CG2,
		     Residue::AL_3HG2, Residue::AL_CG2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CG1, Residue::AL_CG2, Residue::AL_INVALID); 
	END_RES;

	BEGIN_RES(Residue::TYR);// was HD1, HD2, HE1, HE2
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_CD1, Residue::AL_CD2, Residue::AL_1HD, Residue::AL_2HD, Residue::AL_CE1,
		      Residue::AL_CE2, Residue::AL_1HE, Residue::AL_2HE, Residue::AL_CZ, Residue::AL_OH, Residue::AL_HH, Residue::AL_INVALID );
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_CD1, Residue::AL_CG, 
		     Residue::AL_CD2, Residue::AL_CG,
		     // was HD1, HD2
		     Residue::AL_1HD, Residue::AL_CD1, 
		     Residue::AL_2HD, Residue::AL_CD2, 
		     Residue::AL_CE1, Residue::AL_CD1,
		     Residue::AL_CE2, Residue::AL_CD2, 
		     // HE1,2
		     Residue::AL_1HE, Residue::AL_CE1, 
		     Residue::AL_2HE, Residue::AL_CE2,
		     Residue::AL_CZ, Residue::AL_CE1,
		     Residue::AL_CZ, Residue::AL_CE2,
		     Residue::AL_OH, Residue::AL_CZ, 
		     Residue::AL_HH, Residue::AL_OH, 
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CZ, Residue::AL_INVALID); // or Residue::AL_OH
	END_RES;

	BEGIN_RES(Residue::TRP);
	// was HD1, HE1, HE3, HZ2, HZ3
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_CD1, Residue::AL_CD2, Residue::AL_HD, Residue::AL_NE1, Residue::AL_CE2,
		      Residue::AL_CE3, Residue::AL_1HE, Residue::AL_3HE, Residue::AL_CZ2, Residue::AL_CZ3, Residue::AL_2HZ, Residue::AL_3HZ, Residue::AL_CH2,
		      Residue::AL_HH2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_CD1, Residue::AL_CG, 
		     Residue::AL_CD2, Residue::AL_CG,
		     // was 1HD
		     Residue::AL_HD, Residue::AL_CD1, 
		     Residue::AL_NE1, Residue::AL_CD1, 
		     Residue::AL_CE2, Residue::AL_CG,
		     Residue::AL_CD2, Residue::AL_CE2, 
		     Residue::AL_CE3, Residue::AL_CD2, 
		     // was HE1,3
		     Residue::AL_1HE, Residue::AL_NE1,
		     Residue::AL_3HE, Residue::AL_CE3, 
		     Residue::AL_CZ2, Residue::AL_CE2,
		     Residue::AL_CZ3, Residue::AL_CE3,
		     Residue::AL_2HZ, Residue::AL_CZ2,
		     Residue::AL_3HZ, Residue::AL_CZ3, 
		     Residue::AL_CH2, Residue::AL_CZ2,
		     Residue::AL_CH2, Residue::AL_CZ3, 
		     Residue::AL_HH2, Residue::AL_CH2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CH2, Residue::AL_INVALID);
	END_RES;


	BEGIN_RES(Residue::THR);
	// was HG1
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_HB, Residue::AL_CG2, Residue::AL_OG1, Residue::AL_1HG2, Residue::AL_2HG2, Residue::AL_3HG2,
		     Residue::AL_1HG,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_HB, Residue::AL_CB, 
		     Residue::AL_CG2, Residue::AL_CB,
		     Residue::AL_OG1, Residue::AL_CB,
		     Residue::AL_1HG2, Residue::AL_CG2,
		     Residue::AL_2HG2, Residue::AL_CG2,
		     Residue::AL_3HG2, Residue::AL_CG2,
		     // was HG1
		     Residue::AL_1HG, Residue::AL_OG1,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CG2, Residue::AL_OG1, Residue::AL_INVALID); 
	END_RES;

	BEGIN_RES(Residue::SER);
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_OG, Residue::AL_HG,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_OG, Residue::AL_CB,
		     Residue::AL_HG, Residue::AL_OG,
		     Residue::AL_INVALID);
	END_RES;

	BEGIN_RES(Residue::PRO);
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, Residue::AL_CD, Residue::AL_1HD, Residue::AL_2HD,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB,
		     Residue::AL_1HG, Residue::AL_CG,
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_CD, Residue::AL_CG,
		     Residue::AL_1HD, Residue::AL_CD,
		     Residue::AL_2HD, Residue::AL_CD,
		     Residue::AL_CD, Residue::AL_N,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CG, Residue::AL_INVALID); // maybe Residue::AL_CG
	END_RES;
	
	BEGIN_RES(Residue::PHE);
	// was HD1, HD2, HE1,2
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_CD1, Residue::AL_CD2, Residue::AL_1HD, Residue::AL_2HD, Residue::AL_CE1,
		      Residue::AL_CE2, Residue::AL_1HE, Residue::AL_2HE, Residue::AL_CZ, Residue::AL_HZ,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_CD1, Residue::AL_CG,
		     Residue::AL_CD2, Residue::AL_CG,
		     // was HD1, HD2
		     Residue::AL_1HD, Residue::AL_CD1,
		     Residue::AL_2HD, Residue::AL_CD2, 
		     Residue::AL_CE1, Residue::AL_CD1,
		     Residue::AL_CE2, Residue::AL_CD2,
		     // HE1,2
		     Residue::AL_1HE, Residue::AL_CE1, 
		     Residue::AL_2HE, Residue::AL_CE2,
		     Residue::AL_CZ, Residue::AL_CE1, 
		     Residue::AL_CZ, Residue::AL_CE2, Residue::AL_HZ,
		     Residue::AL_CZ,Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CZ, Residue::AL_INVALID); //
	END_RES;
	
	BEGIN_RES(Residue::NH2);
	DEFINE_ATOMS(Residue::AL_N, Residue::AL_1H, Residue::AL_2H,Residue::AL_INVALID);
	DEFINE_BONDS( Residue::AL_1H, Residue::AL_N, 
		      Residue::AL_2H, Residue::AL_N,
		      Residue::AL_INVALID);
	END_RES;
	
	BEGIN_RES(Residue::MET);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, 
		     Residue::AL_SD, Residue::AL_CE, Residue::AL_1HE,
		     Residue::AL_2HE, Residue::AL_3HE,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_1HG, Residue::AL_CG,
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_SD, Residue::AL_CG,
		     Residue::AL_CE, Residue::AL_SD, 
		     Residue::AL_1HE, Residue::AL_CE,
		     Residue::AL_2HE, Residue::AL_CE, 
		     Residue::AL_3HE, Residue::AL_CE,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CE, Residue::AL_INVALID); 
	END_RES;

	BEGIN_RES(Residue::LYS);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, Residue::AL_CD, Residue::AL_1HD, Residue::AL_2HD,
		     Residue::AL_CE, Residue::AL_1HE, Residue::AL_2HE, Residue::AL_NZ, Residue::AL_1HZ, Residue::AL_2HZ, Residue::AL_3HZ,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB,
		     Residue::AL_1HG, Residue::AL_CG, 
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_CD, Residue::AL_CG,
		     Residue::AL_1HD,Residue::AL_CD, 
		     Residue::AL_2HD, Residue::AL_CD,
		     Residue::AL_CE, Residue::AL_CD, 
		     Residue::AL_1HE, Residue::AL_CE, 
		     Residue::AL_2HE, Residue::AL_CE,
		     Residue::AL_NZ, Residue::AL_CE,
		     Residue::AL_1HZ, Residue::AL_NZ, 
		     Residue::AL_2HZ, Residue::AL_NZ,
		     Residue::AL_3HZ, Residue::AL_NZ,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_NZ, Residue::AL_INVALID);
	END_RES;
	
	BEGIN_RES(Residue::LEU);
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_HG, Residue::AL_CD1, Residue::AL_CD2, Residue::AL_1HD1, Residue::AL_2HD1,
		      Residue::AL_3HD1, Residue::AL_1HD2, Residue::AL_2HD2, Residue::AL_3HD2,Residue::AL_INVALID);
	DEFINE_BONDS( Residue::AL_CB, Residue::AL_CA,
		      Residue::AL_1HB, Residue::AL_CB,
		      Residue::AL_2HB, Residue::AL_CB,
		      Residue::AL_CG, Residue::AL_CB, 
		      Residue::AL_HG, Residue::AL_CG,
		      Residue::AL_CD1, Residue::AL_CG,
		      Residue::AL_CD2, Residue::AL_CG, 
		      Residue::AL_1HD1, Residue::AL_CD1,
		      Residue::AL_2HD1, Residue::AL_CD1,
		      Residue::AL_3HD1, Residue::AL_CD1,
		      Residue::AL_1HD2, Residue::AL_CD2,
		      Residue::AL_2HD2, Residue::AL_CD2,
		      Residue::AL_3HD2, Residue::AL_CD2,
		      Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CD1, Residue::AL_CD2, Residue::AL_INVALID);
	END_RES;

	BEGIN_RES(Residue::ILE);
	// was 1HD
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_HB, Residue::AL_CG1, Residue::AL_CG2, Residue::AL_1HG1, Residue::AL_2HG1, Residue::AL_CD, Residue::AL_1HD, Residue::AL_2HD,
		     Residue::AL_3HD, Residue::AL_1HG2, Residue::AL_2HG2, Residue::AL_3HG2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_HB, Residue::AL_CB,
		     Residue::AL_CG1, Residue::AL_CB,
		     Residue::AL_CG2, Residue::AL_CB,
		     Residue::AL_1HG1, Residue::AL_CG1,
		     Residue::AL_2HG1, Residue::AL_CG1,
		     Residue::AL_CD, Residue::AL_CG1, 
		     Residue::AL_1HD, Residue::AL_CD,
		     Residue::AL_2HD, Residue::AL_CD,
		     Residue::AL_3HD, Residue::AL_CD,
		     Residue::AL_1HG2, Residue::AL_CG2,
		     Residue::AL_2HG2, Residue::AL_CG2,
		     Residue::AL_3HG2, Residue::AL_CG2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CD, Residue::AL_INVALID); // maybe add CG2
	END_RES;
	
	BEGIN_RES(Residue::HIS);
	// was HD1, HD2 HE1
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_ND1, Residue::AL_CD2, Residue::AL_1HD, Residue::AL_2HD, Residue::AL_CE1,
		     Residue::AL_NE2, Residue::AL_1HE,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_ND1, Residue::AL_CG,
		     Residue::AL_CD2, Residue::AL_CG,
		     Residue::AL_1HD, Residue::AL_ND1, 
		     Residue::AL_2HD, Residue::AL_CD2,
		     Residue::AL_CE1, Residue::AL_ND1,
		     Residue::AL_CE1, Residue::AL_NE2, 
		     Residue::AL_1HE, Residue::AL_CE1, 
		     Residue::AL_NE2, Residue::AL_CD2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CE1, Residue::AL_INVALID);
	END_RES;
	
	
	BEGIN_RES(Residue::GLY);
	DEFINE_ATOMS(Residue::AL_2HA,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_2HA, Residue::AL_CA,Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_INVALID);
	END_RES;	
	
	BEGIN_RES(Residue::GLU);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, Residue::AL_CD, Residue::AL_OE1, Residue::AL_OE2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB,
		     Residue::AL_1HG, Residue::AL_CG, 
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_CD, Residue::AL_CG, 
		     Residue::AL_OE1, Residue::AL_CD, 
		     Residue::AL_OE2, Residue::AL_CD,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CD, Residue::AL_INVALID); // maybe should be Residue::AL_OE[12]
	END_RES;	
	
	BEGIN_RES(Residue::GLN);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, Residue::AL_CD, Residue::AL_OE1, Residue::AL_NE2,
		     Residue::AL_1HE2, Residue::AL_2HE2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_1HG, Residue::AL_CG,
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_CD, Residue::AL_CG, 
		     Residue::AL_OE1, Residue::AL_CD, 
		     Residue::AL_NE2, Residue::AL_CD,
		     Residue::AL_1HE2, Residue::AL_NE2, 
		     Residue::AL_2HE2, Residue::AL_NE2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CD, Residue::AL_INVALID); // or Residue::AL_NE2, Residue::AL_OE1
	END_RES;	
	
	BEGIN_RES(Residue::CYS);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_SG, Residue::AL_HG,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_SG, Residue::AL_CB, 
		     Residue::AL_HG, Residue::AL_SG,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_SG, Residue::AL_INVALID);
	END_RES;	
	
	BEGIN_RES(Residue::ASP);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_OD1, Residue::AL_OD2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA,
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_OD1, Residue::AL_CG, 
		     Residue::AL_OD2, Residue::AL_CG,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CG, Residue::AL_INVALID); // or Residue::AL_OD[12]
	END_RES;	
	
	BEGIN_RES(Residue::ASN);
	DEFINE_ATOMS( Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_OD1, Residue::AL_ND2, Residue::AL_1HD2, Residue::AL_2HD2,Residue::AL_INVALID);
	DEFINE_BONDS( Residue::AL_CB, Residue::AL_CA,
		      Residue::AL_1HB, Residue::AL_CB, 
		      Residue::AL_2HB, Residue::AL_CB,
		      Residue::AL_CG, Residue::AL_CB,
		      Residue::AL_OD1, Residue::AL_CG,
		      Residue::AL_ND2, Residue::AL_CG,
		      Residue::AL_1HD2, Residue::AL_ND2, 
		      Residue::AL_2HD2, Residue::AL_ND2 ,
		      Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CG, Residue::AL_INVALID); // or Residue::AL_OD1, Residue::AL_ND1
	END_RES;	
	
	BEGIN_RES(Residue::ARG);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_CG, Residue::AL_1HG, Residue::AL_2HG, Residue::AL_CD, Residue::AL_1HD, Residue::AL_2HD,
		     Residue::AL_NE, Residue::AL_HE, Residue::AL_CZ, Residue::AL_NH1, Residue::AL_NH2, Residue::AL_1HH1, Residue::AL_2HH1, Residue::AL_1HH2,
		     Residue::AL_2HH2,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CB, Residue::AL_CA, 
		     Residue::AL_1HB, Residue::AL_CB, 
		     Residue::AL_2HB, Residue::AL_CB,
		     Residue::AL_CG, Residue::AL_CB, 
		     Residue::AL_1HG, Residue::AL_CG, 
		     Residue::AL_2HG, Residue::AL_CG,
		     Residue::AL_CD, Residue::AL_CG,
		     Residue::AL_1HD, Residue::AL_CD,
		     Residue::AL_2HD, Residue::AL_CD,
		     Residue::AL_NE, Residue::AL_CD, 
		     Residue::AL_HE, Residue::AL_NE, 
		     Residue::AL_CZ, Residue::AL_NE,
		     Residue::AL_NH1, Residue::AL_CZ, 
		     Residue::AL_NH2, Residue::AL_CZ, 
		     Residue::AL_1HH1, Residue::AL_NH1,
		     Residue::AL_2HH1, Residue::AL_NH1, 
		     Residue::AL_1HH2, Residue::AL_NH2, 
		     Residue::AL_2HH2, Residue::AL_NH2,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CZ, Residue::AL_INVALID); // maybe Residue::AL_NH2m Residue::AL_NH1
	END_RES;	
	
	BEGIN_RES(Residue::ALA);
	DEFINE_ATOMS(Residue::AL_CB, Residue::AL_1HB, Residue::AL_2HB, Residue::AL_3HB,Residue::AL_INVALID);
	DEFINE_BONDS( Residue::AL_CB, Residue::AL_CA,
		      Residue::AL_1HB, Residue::AL_CB, 
		      Residue::AL_2HB, Residue::AL_CB,
		      Residue::AL_3HB, Residue::AL_CB,
		      Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_CB, Residue::AL_INVALID);
	END_RES;	
	
	BEGIN_RES(Residue::ACE);
	DEFINE_ATOMS(Residue::AL_CH3, Residue::AL_C, Residue::AL_O, Residue::AL_1HH3, Residue::AL_2HH3, Residue::AL_3HH3,Residue::AL_INVALID);
	DEFINE_BONDS(Residue::AL_CH3, Residue::AL_C, 
		     Residue::AL_O, Residue::AL_C,
		     Residue::AL_C, Residue::AL_1HH3,
		     Residue::AL_C, Residue::AL_2HH3, 
		     Residue::AL_C, Residue::AL_3HH3,
		     Residue::AL_INVALID);
	DEFINE_EXTREMES(Residue::AL_INVALID);
	END_RES;	
      default:
	;
      }
    
      return dat;
      //return std::pair<Residue::Atom_label*, Residue::Atom_label*>(als, bls);
    }

    void do_initialize() {
      //if (amino_acid_initialized_) return;
      assert(amino_acid_initialized_==false);
      amino_acid_initialized_=true;
      unsigned int i=0;
      for (; atom_name_data_[i].l != Residue::AL_INVALID; ++i){
	clean_atom_name_data_[i].l= atom_name_data_[i].l;
	clean_atom_name_data_[i].t= atom_name_data_[i].t;
	sscanf(atom_name_data_[i].s, "%4s", clean_atom_name_data_[i].s);
      }
      clean_atom_name_data_[i].l= Residue::AL_INVALID;

      unsigned int max_atoms=0;
      for (unsigned int i=0;atom_fallback_data_[i].l != Residue::AL_INVALID; ++i){
	max_atoms= std::max(max_atoms, static_cast<unsigned int>(atom_fallback_data_[i].l));
      }
      clean_atom_fallbacks_.resize(max_atoms+1);

      for (unsigned int i=0;atom_fallback_data_[i].l != Residue::AL_INVALID; ++i){
	clean_atom_fallbacks_[atom_fallback_data_[i].l].push_back(atom_fallback_data_[i].lr);
      }


      Residue::Type all_res[]={Residue::GLY, Residue::ALA, 
			       Residue::VAL, Residue::LEU, 
			       Residue::ILE, Residue::SER, 
			       Residue::THR, Residue::CYS, 
			       Residue::MET, Residue::PRO,
			       Residue::ASP, Residue::ASN, 
			       Residue::GLU, Residue::GLN, 
			       Residue::LYS, Residue::ARG,
			       Residue::HIS, Residue::PHE, 
			       Residue::TYR, Residue::TRP,
			       Residue::ACE, Residue::NH2, 
			       Residue::INV};
      unsigned int num_res= sizeof(all_res)/sizeof(Residue::Type);


      amino_acid_data_.resize(num_res);


      for (unsigned int i=0; i< num_res; ++i){
	Residue::Type cur_res=all_res[i];
	amino_acid_data_[cur_res]= Amino_acid_data();
	if (cur_res != Residue::ACE && cur_res != Residue::NH2) {
	  // set up shared atoms
	  Residue::Atom_label bl[]={ Residue::AL_N, Residue::AL_H,
				     Residue::AL_1H, Residue::AL_2H, 
				     Residue::AL_3H, Residue::AL_CA, 
				     Residue::AL_HA, Residue::AL_1HA,
				     Residue::AL_C, Residue::AL_O,
				     Residue::AL_OXT};
	  for (unsigned int j=0; j< sizeof(bl)/sizeof(Residue::Atom_label); ++j){
	    amino_acid_data_[cur_res].atoms.push_back(bl[j]);
	  }
	  // set up basic bonds
	  Residue::Atom_label bd[][2]={ {Residue::AL_N, Residue::AL_H}, {Residue::AL_N, Residue::AL_1H},
					{Residue::AL_N, Residue::AL_2H}, {Residue::AL_N, Residue::AL_3H},
					{Residue::AL_N, Residue::AL_CA}, {Residue::AL_CA, Residue::AL_HA},
					{Residue::AL_CA, Residue::AL_1HA}, {Residue::AL_CA, Residue::AL_C},
					{Residue::AL_C, Residue::AL_O}, {Residue::AL_O, Residue::AL_OXT} };
	  for (unsigned int j=0; j< sizeof(bd)/sizeof(Possible_bond); ++j){
	    amino_acid_data_[cur_res].bonds.push_back(Possible_bond(bd[j][0], bd[j][1]));
	  }
	}

	/*
	  (Residue::Residue::AL_CB, Residue::Residue::AL_HB, Residue::Residue::AL_CG1,  Residue::Residue::AL_CG2, 
	  Residue::Residue::AL_1HG1, Residue::Residue::AL_2HG1, Residue::Residue::AL_3HG1, Residue::Residue::AL_1HG2,
	  Residue::Residue::AL_2HG2, Residue::Residue::AL_3HG2, Residue::Residue::AL_INVALID) ,
	  ()
	*/

	Residue_init_data data= get_residue_initialization_data(cur_res);
      
      
	for (unsigned int j=0; data.atms_[j] != Residue::AL_INVALID; ++j){
	  amino_acid_data_[cur_res].atoms.push_back(data.atms_[j]);
	}
	for (unsigned int j=0; data.bnds_[j] != Residue::AL_INVALID; j+=2){
	  amino_acid_data_[cur_res].bonds.push_back(Possible_bond(data.bnds_[j], data.bnds_[j+1]));
	}
	for (unsigned int j=0; data.extr_[j] != Residue::AL_INVALID; ++j){
	  amino_acid_data_[cur_res].extremes.push_back(data.extr_[j]);
	}

	if (false){
	  std::cout << "For residue " << Residue::type_string(cur_res)
		    << " there are " << amino_acid_data_[cur_res].atoms.size() << " atoms, "
		    << amino_acid_data_[cur_res].bonds.size() << " bonds, and "
		    << amino_acid_data_[cur_res].extremes.size() << " extremes." << std::endl;
	  assert(amino_acid_data_[cur_res].extremes.size() < 4);
	}
      
      }
    }


    Residue::Atom_label fix_atom_label(Residue::Type label, Residue::Atom_label al) {
      for (unsigned int i=0; i< amino_acid_data_[label].atoms.size(); ++i){
	if (Residue_data::amino_acid_data_[label].atoms[i] == al) return al;
      }
      //HERE
      /*if (clean_atom_fallbacks_.find(al) == clean_atom_fallbacks_.end()){
	dsrpdb_internal::error_logger.new_warning((std::string("No atom named ")+atom_label_string(al) 
	+" in residue "+ type_string(label_)).c_str());
	return AL_INVALID;
	}*/
      for (unsigned int j=0; j< clean_atom_fallbacks_[al].size(); ++j){
	for (unsigned int i=0; i< amino_acid_data_[label].atoms.size(); ++i){
	  if (amino_acid_data_[label].atoms[i] == clean_atom_fallbacks_[al][j]) 
	    return clean_atom_fallbacks_[al][j];
	}
      }
      return Residue::AL_INVALID;
    }


  }

}
