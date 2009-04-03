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

#include <CGAL/PDB/Monomer.h>
#include <CGAL/PDB/internal/Monomer_data.h>
#include <cassert>
#include <cstdio>
namespace CGAL { namespace PDB {

namespace Monomer_data {
  bool amino_acid_initialized_=false;
  typedef std::vector<std::vector<Monomer::Atom_key> > Clean_atom_fallbacks;
  std::vector<Amino_acid_data > amino_acid_data_;



  
  //! The per atom data associating string names with labels and atom types
  /*!  To add a new atom or a different string refering to the same
    atom add a line to this array.
  */
  Atom_data atom_name_data_[]= {
    {" N  ",Monomer::AL_N, Atom::N},
    {" H  ",Monomer::AL_H, Atom::H},
    {"1H  ",Monomer::AL_1H, Atom::H},{" H1 ",Monomer::AL_1H, Atom::H},
    {"2H  ",Monomer::AL_2H, Atom::H},{" H2 ",Monomer::AL_2H, Atom::H},
    {"3H  ",Monomer::AL_3H, Atom::H},{" H3 ",Monomer::AL_3H, Atom::H},
    {" C  ",Monomer::AL_C, Atom::C},
    {" O  ",Monomer::AL_O, Atom::O},
    {" OXT",Monomer::AL_OXT, Atom::O},
    {" CH3",Monomer::AL_CH3, Atom::C},

    {" CA ",Monomer::AL_CA, Atom::C},
    {" HA ",Monomer::AL_HA, Atom::H},
    {"1HA ",Monomer::AL_1HA, Atom::H},{" HA1",Monomer::AL_1HA, Atom::H},
    {"2HA ",Monomer::AL_2HA, Atom::H},{" HA2",Monomer::AL_2HA, Atom::H},

    {" CB ",Monomer::AL_CB, Atom::C},
    {" HB ",Monomer::AL_HB, Atom::H},
    {"1HB ",Monomer::AL_1HB, Atom::H}, {" HB1",Monomer::AL_1HB, Atom::H}, 
    {"2HB ", Monomer::AL_2HB, Atom::H}, {" HB2", Monomer::AL_2HB, Atom::H},
    {"3HB ", Monomer::AL_3HB, Atom::H}, {" HB3", Monomer::AL_3HB, Atom::H},

    {" CG ",Monomer::AL_CG, Atom::C},
    {" CG1",Monomer::AL_CG1, Atom::C},
    {" CG2",Monomer::AL_CG2, Atom::C},
    {" HG ",Monomer::AL_HG, Atom::H},
    {"1HG ", Monomer::AL_1HG, Atom::H}, {" HG1", Monomer::AL_1HG, Atom::H},
    {"2HG ", Monomer::AL_2HG, Atom::H}, {" HG2", Monomer::AL_2HG, Atom::H},
    //{"HG1",Monomer::AL_HG1},
    {"1HG1",Monomer::AL_1HG1, Atom::H},{"HG11",Monomer::AL_1HG1, Atom::H},
    {"2HG1",Monomer::AL_2HG1, Atom::H},{"HG12",Monomer::AL_2HG1, Atom::H},
    {"3HG1",Monomer::AL_3HG1, Atom::H},{"HG13",Monomer::AL_3HG1, Atom::H},
    {"1HG2",Monomer::AL_1HG2, Atom::H},{"HG21",Monomer::AL_1HG2, Atom::H},
    {"2HG2",Monomer::AL_2HG2, Atom::H},{"HG22",Monomer::AL_2HG2, Atom::H},
    {"3HG2",Monomer::AL_3HG2, Atom::H},{"HG23",Monomer::AL_3HG2, Atom::H},
    {" OG ",Monomer::AL_OG, Atom::O},
    {" OG1",Monomer::AL_OG1, Atom::O},
    {" SG ",Monomer::AL_SG, Atom::S},

    {" CD ",Monomer::AL_CD, Atom::C},
    {" CD1",Monomer::AL_CD1, Atom::C},
    {" CD2",Monomer::AL_CD2, Atom::C},
    //{"HD1",Monomer::AL_HD1},
    //{"HD2",Monomer::AL_HD2},
    {" HD ",Monomer::AL_HD, Atom::H},
    {"1HD ",Monomer::AL_1HD, Atom::H},{" HD1",Monomer::AL_1HD, Atom::H},
    {"2HD ",Monomer::AL_2HD, Atom::H},{" HD2",Monomer::AL_2HD, Atom::H},
    {"3HD ",Monomer::AL_3HD, Atom::H},{" HD3",Monomer::AL_3HD, Atom::H},
    {"1HD1",Monomer::AL_1HD1, Atom::H}, {"HD11",Monomer::AL_1HD1, Atom::H}, 
    {"2HD1",Monomer::AL_2HD1, Atom::H},{"HD12",Monomer::AL_2HD1, Atom::H},
    {"3HD1",Monomer::AL_3HD1, Atom::H},{"HD13",Monomer::AL_3HD1, Atom::H},
    {"1HD2",Monomer::AL_1HD2, Atom::H},{"HD21",Monomer::AL_1HD2, Atom::H},
    {"2HD2",Monomer::AL_2HD2, Atom::H},{"HD22",Monomer::AL_2HD2, Atom::H},
    {"3HD2",Monomer::AL_3HD2, Atom::H},{"HD23",Monomer::AL_3HD2, Atom::H},
    {" SD ",Monomer::AL_SD, Atom::S},
    {" OD1",Monomer::AL_OD1, Atom::O},
    {" OD2",Monomer::AL_OD2, Atom::O},
    {" ND1",Monomer::AL_ND1, Atom::N},
    {" ND2",Monomer::AL_ND2, Atom::N},

    {" CE ",Monomer::AL_CE, Atom::C},
    {" CE1",Monomer::AL_CE1, Atom::C},
    {" CE2",Monomer::AL_CE2, Atom::C},
    {" CE3",Monomer::AL_CE3, Atom::C},
    {" HE ",Monomer::AL_HE, Atom::H},
    {"1HE ",Monomer::AL_1HE, Atom::H},{" HE1",Monomer::AL_1HE, Atom::H},
    {"2HE ",Monomer::AL_2HE, Atom::H},{" HE2",Monomer::AL_2HE, Atom::H},
    {"3HE ",Monomer::AL_3HE, Atom::H},{" HE3",Monomer::AL_3HE, Atom::H},
    //{"HE1",Monomer::AL_HE1},
    //{"HE2",Monomer::AL_HE2},
    //{"HE3",Monomer::AL_HE3},
    {"1HE2",Monomer::AL_1HE2, Atom::H},{"HE21",Monomer::AL_1HE2, Atom::H},
    {"2HE2",Monomer::AL_2HE2, Atom::H},{"HE22",Monomer::AL_2HE2, Atom::H},
    {" OE1",Monomer::AL_OE1, Atom::O},
    {" OE2",Monomer::AL_OE2, Atom::O},
    {" NE ",Monomer::AL_NE, Atom::N},
    {" NE1",Monomer::AL_NE1, Atom::N},
    {" NE2",Monomer::AL_NE2, Atom::N},

    {" CZ ",Monomer::AL_CZ, Atom::C},
    {" CZ2",Monomer::AL_CZ2, Atom::C},
    {" CZ3",Monomer::AL_CZ3, Atom::C},
    {" NZ ",Monomer::AL_NZ, Atom::N},
    {" HZ ",Monomer::AL_HZ, Atom::H},
    {"1HZ ",Monomer::AL_1HZ, Atom::H},{" HZ1",Monomer::AL_1HZ, Atom::H},
    {"2HZ ",Monomer::AL_2HZ, Atom::H},{" HZ2",Monomer::AL_2HZ, Atom::H},
    {"3HZ ",Monomer::AL_3HZ, Atom::H},{" HZ3",Monomer::AL_3HZ, Atom::H},
    //{"HZ1",Monomer::AL_HZ2},
    //{"HZ2",Monomer::AL_HZ2},
    //{"HZ3",Monomer::AL_HZ3},

    {" CH2",Monomer::AL_CH2, Atom::C},
    {" NH1",Monomer::AL_NH1, Atom::N},
    {" NH2",Monomer::AL_NH2, Atom::N},
    {" OH",Monomer::AL_OH, Atom::O},
    {" HH",Monomer::AL_HH, Atom::H},
      
    {"1HH1", Monomer::AL_1HH1, Atom::H}, {"HH11", Monomer::AL_1HH1, Atom::H},
    {"2HH1", Monomer::AL_2HH1, Atom::H}, {"HH12", Monomer::AL_2HH1, Atom::H},
    {" HH2",Monomer::AL_HH2, Atom::H},
    {"1HH2", Monomer::AL_1HH2, Atom::H}, {"HH21", Monomer::AL_1HH2, Atom::H},
    {"2HH2", Monomer::AL_2HH2, Atom::H}, {"HH22", Monomer::AL_2HH2, Atom::H},
    {" HH ", Monomer::AL_1HH2, Atom::H},

    {"HH31", Monomer::AL_1HH3, Atom::H}, 
    {"HH32", Monomer::AL_2HH3, Atom::H}, 
    {"HH33", Monomer::AL_3HH3, Atom::H},

    {"P   ", Monomer::AL_P, Atom::P},
    {"O1P ", Monomer::AL_OP1, Atom::O},
    {"O2P ", Monomer::AL_OP2, Atom::O},
    {"O5* ", Monomer::AL_O5p, Atom::O},
    {"C5* ", Monomer::AL_H5p, Atom::C}, 
    {"H5**", Monomer::AL_H5pp, Atom::H},
    {"C4* ", Monomer::AL_C4p, Atom::C},
    {"H4* ", Monomer::AL_H4p, Atom::H},
    {"O4* ", Monomer::AL_O4p, Atom::O},
    {"C1* ", Monomer::AL_C1p, Atom::C},
    {"H1* ", Monomer::AL_H1p, Atom::H},
    {"C3* ", Monomer::AL_C3p, Atom::C},
    {"H3* ", Monomer::AL_H3p, Atom::H},
    {"O3* ", Monomer::AL_O3p, Atom::O},
    {"C2* ", Monomer::AL_C2p, Atom::C},
    {"H2* ", Monomer::AL_H2p, Atom::H},
    {"H2**", Monomer::AL_H2pp, Atom::H},
    {"O2* ", Monomer::AL_O2p, Atom::O},
    {"HO2*", Monomer::AL_HO2p, Atom::O},
    {"N9 ", Monomer::AL_N9, Atom::N},
    {"C8 ", Monomer::AL_C8, Atom::C},
    {"H8 ", Monomer::AL_H8, Atom::H},
    {"N7 ", Monomer::AL_N7, Atom::N},
    {"C5 ", Monomer::AL_C5, Atom::C},
    {"C4 ", Monomer::AL_C4, Atom::C},
    {"N3 ", Monomer::AL_N3, Atom::N},
    {"C2 ", Monomer::AL_C2, Atom::C},
    {"H2 ", Monomer::AL_H2, Atom::H},
    {"N1 ", Monomer::AL_N1, Atom::N},
    {"C6 ", Monomer::AL_C6, Atom::C},
    {"N6 ", Monomer::AL_N6, Atom::N},
    {"H61 ", Monomer::AL_H61, Atom::H},
    {"H62 ", Monomer::AL_H62, Atom::H},
    {"O6  ", Monomer::AL_O6, Atom::O},
    {"H1  ", Monomer::AL_H1, Atom::H},
    {"N2  ", Monomer::AL_N2, Atom::N},
    {"H21 ", Monomer::AL_H21, Atom::H},
    {"H22 ", Monomer::AL_H22, Atom::H},
    
    {"H6  ", Monomer::AL_H6, Atom::H},
    {"H5  ", Monomer::AL_H5, Atom::H},
    {"O2  ", Monomer::AL_O2, Atom::O},
    {"N4  ", Monomer::AL_N4, Atom::N},
    {"H41 ", Monomer::AL_H41, Atom::H},
    {"H42 ", Monomer::AL_H42,Atom::H},
    {"H3  ", Monomer::AL_H3, Atom::H},
    {"O4  ", Monomer::AL_O4, Atom::O},
    {"C7  ", Monomer::AL_C7, Atom::C},
    {"H71 ", Monomer::AL_H71, Atom::H},
    {"H72 ", Monomer::AL_H72,Atom::H},
    {"H73 ", Monomer::AL_H73,Atom::H},

    {"UNKN", Monomer::AL_INVALID, Atom::INVALID}};
    

  Atom_data clean_atom_name_data_[sizeof(atom_name_data_)
				  /sizeof(Atom_data)+1];

  Atom_fallback_data atom_fallback_data_[]=
    {{Monomer::AL_CD1, Monomer::AL_CD}, 
     {Monomer::AL_1HA, Monomer::AL_HA},
     {Monomer::AL_1HB, Monomer::AL_HB},
     {Monomer::AL_1HD, Monomer::AL_HD},
     {Monomer::AL_1HE, Monomer::AL_HE},
     {Monomer::AL_1HZ, Monomer::AL_HZ},
     //{Monomer::AL_1HH, Monomer::AL_HH},
     {Monomer::AL_1HD1, Monomer::AL_1HD},
     {Monomer::AL_2HD1, Monomer::AL_2HD},
     {Monomer::AL_3HD1, Monomer::AL_3HD},
     {Monomer::AL_INVALID, Monomer::AL_INVALID}};
  Clean_atom_fallbacks clean_atom_fallbacks_;

#define BEGIN_RES(name) case name:{
#define BEGIN_ATOMS static Monomer::Atom_key cl[]=
#define END_ATOMS ; dat.atms_=cl;
#define BEGIN_BONDS static Monomer::Atom_key cbl[]=
#define END_BONDS ; for (int i=0; cbl[i] != Monomer::AL_INVALID; ++i) {\
    bool found=false;							\
    for (int j=0; cl[j] != Monomer::AL_INVALID; ++j) {			\
      if (cl[j]==cbl[i]) found=true;					\
    }\
    assert(found || cbl[i]== Monomer::AL_CA || cbl[i]== Monomer::AL_N ||cbl[i]== Monomer::AL_C1p); \
  };									\
    dat.bnds_=cbl;
#define BEGIN_EXTREMES static Monomer::Atom_key ext[]=
#define END_EXTREMES ; for (int i=0; ext[i] != Monomer::AL_INVALID; ++i) { bool found=false; for (int j=0; cl[j] != Monomer::AL_INVALID; ++j) if (cl[j]==ext[i]) found=true; assert(found);} dat.extr_=ext; 
 
  /*#define DEFINE_ATOMS(...) static Monomer::Atom_key cl[]={__VA_ARGS__}; dat.atms_=cl;
#define DEFINE_BONDS(...) static Monomer::Atom_key cbl[]={__VA_ARGS__}; for (int i=0; cbl[i] != Monomer::AL_INVALID; ++i) { bool found=false; for (int j=0; cl[j] != Monomer::AL_INVALID; ++j) if (cl[j]==cbl[i]) found=true; assert(found || cbl[i]== Monomer::AL_CA || cbl[i]== Monomer::AL_N);}; dat.bnds_=cbl; 
#define DEFINE_EXTREMES(...) static Monomer::Atom_key ext[]={__VA_ARGS__}; for (int i=0; ext[i] != Monomer::AL_INVALID; ++i) { bool found=false; for (int j=0; cl[j] != Monomer::AL_INVALID; ++j) if (cl[j]==ext[i]) found=true; assert(found);} dat.extr_=ext; */
#define END_RES break;}


  //! This function returns the per-residue atom and bond information.
  /*!
    To add an atom or bond to a residue add them to the appropriate list below.
  */
  Monomer_init_data 
  get_residue_initialization_data(Monomer::Type rl){
    //using class Monomer;
    static Monomer::Atom_key end=Monomer::AL_INVALID;
    Monomer_init_data dat;
    dat.atms_=&end;
    dat.bnds_=&end;
    dat.extr_=&end;
    switch(rl){
      
      BEGIN_RES(Monomer::VAL);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_HB, Monomer::AL_CG1,  Monomer::AL_CG2, 
	  Monomer::AL_1HG1, Monomer::AL_2HG1, Monomer::AL_3HG1, Monomer::AL_1HG2,
	  Monomer::AL_2HG2, Monomer::AL_3HG2, Monomer::AL_INVALID }END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
	  Monomer::AL_HB, Monomer::AL_CB, 
		   Monomer::AL_CG1, Monomer::AL_CB,
		   Monomer::AL_CG2, Monomer::AL_CB, 
		   Monomer::AL_1HG1, Monomer::AL_CG1, 
		   Monomer::AL_2HG1, Monomer::AL_CG1,
		   Monomer::AL_3HG1, Monomer::AL_CG1, 
		   Monomer::AL_1HG2, Monomer::AL_CG2, 
		   Monomer::AL_2HG2, Monomer::AL_CG2,
		   Monomer::AL_3HG2, Monomer::AL_CG2,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CG1, Monomer::AL_CG2, Monomer::AL_INVALID}END_EXTREMES; 
      END_RES;

      BEGIN_RES(Monomer::TYR);// was HD1, HD2, HE1, HE2
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_CD1, Monomer::AL_CD2, Monomer::AL_1HD, Monomer::AL_2HD, Monomer::AL_CE1,
	  Monomer::AL_CE2, Monomer::AL_1HE, Monomer::AL_2HE, Monomer::AL_CZ, Monomer::AL_OH, Monomer::AL_HH, Monomer::AL_INVALID } END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_CD1, Monomer::AL_CG, 
		   Monomer::AL_CD2, Monomer::AL_CG,
		   // was HD1, HD2
		   Monomer::AL_1HD, Monomer::AL_CD1, 
		   Monomer::AL_2HD, Monomer::AL_CD2, 
		   Monomer::AL_CE1, Monomer::AL_CD1,
		   Monomer::AL_CE2, Monomer::AL_CD2, 
		   // HE1,2
		   Monomer::AL_1HE, Monomer::AL_CE1, 
		   Monomer::AL_2HE, Monomer::AL_CE2,
		   Monomer::AL_CZ, Monomer::AL_CE1,
		   Monomer::AL_CZ, Monomer::AL_CE2,
		   Monomer::AL_OH, Monomer::AL_CZ, 
		   Monomer::AL_HH, Monomer::AL_OH, 
		   Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CZ, Monomer::AL_INVALID}END_EXTREMES; // or Monomer::AL_OH
      END_RES;

      BEGIN_RES(Monomer::TRP);
      // was HD1, HE1, HE3, HZ2, HZ3
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_CD1, Monomer::AL_CD2, Monomer::AL_HD, Monomer::AL_NE1, Monomer::AL_CE2,
		    Monomer::AL_CE3, Monomer::AL_1HE, Monomer::AL_3HE, Monomer::AL_CZ2, Monomer::AL_CZ3, Monomer::AL_2HZ, Monomer::AL_3HZ, Monomer::AL_CH2,
	  Monomer::AL_HH2,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_CD1, Monomer::AL_CG, 
		   Monomer::AL_CD2, Monomer::AL_CG,
		   // was 1HD
		   Monomer::AL_HD, Monomer::AL_CD1, 
		   Monomer::AL_NE1, Monomer::AL_CD1, 
		   Monomer::AL_CE2, Monomer::AL_CG,
		   Monomer::AL_CD2, Monomer::AL_CE2, 
		   Monomer::AL_CE3, Monomer::AL_CD2, 
		   // was HE1,3
		   Monomer::AL_1HE, Monomer::AL_NE1,
		   Monomer::AL_3HE, Monomer::AL_CE3, 
		   Monomer::AL_CZ2, Monomer::AL_CE2,
		   Monomer::AL_CZ3, Monomer::AL_CE3,
		   Monomer::AL_2HZ, Monomer::AL_CZ2,
		   Monomer::AL_3HZ, Monomer::AL_CZ3, 
		   Monomer::AL_CH2, Monomer::AL_CZ2,
		   Monomer::AL_CH2, Monomer::AL_CZ3, 
		   Monomer::AL_HH2, Monomer::AL_CH2,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CH2, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;


      BEGIN_RES(Monomer::THR);
      // was HG1
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_HB, Monomer::AL_CG2, Monomer::AL_OG1, Monomer::AL_1HG2, Monomer::AL_2HG2, Monomer::AL_3HG2,
	  Monomer::AL_1HG,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_HB, Monomer::AL_CB, 
		   Monomer::AL_CG2, Monomer::AL_CB,
		   Monomer::AL_OG1, Monomer::AL_CB,
		   Monomer::AL_1HG2, Monomer::AL_CG2,
		   Monomer::AL_2HG2, Monomer::AL_CG2,
		   Monomer::AL_3HG2, Monomer::AL_CG2,
		   // was HG1
		   Monomer::AL_1HG, Monomer::AL_OG1,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CG2, Monomer::AL_OG1, Monomer::AL_INVALID}END_EXTREMES; 
      END_RES;

      BEGIN_RES(Monomer::SER);
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_OG, Monomer::AL_HG,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_OG, Monomer::AL_CB,
		   Monomer::AL_HG, Monomer::AL_OG,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_OG, Monomer::AL_INVALID}END_EXTREMES; 
      END_RES;

      BEGIN_RES(Monomer::PRO);
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, Monomer::AL_CD, Monomer::AL_1HD, Monomer::AL_2HD,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB,
		   Monomer::AL_1HG, Monomer::AL_CG,
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_CD, Monomer::AL_CG,
		   Monomer::AL_1HD, Monomer::AL_CD,
		   Monomer::AL_2HD, Monomer::AL_CD,
		   Monomer::AL_CD, Monomer::AL_N,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CG, Monomer::AL_INVALID}END_EXTREMES; // maybe Monomer::AL_CG
      END_RES;
	
      BEGIN_RES(Monomer::PHE);
      // was HD1, HD2, HE1,2
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_CD1, Monomer::AL_CD2, Monomer::AL_1HD, Monomer::AL_2HD, Monomer::AL_CE1,
	  Monomer::AL_CE2, Monomer::AL_1HE, Monomer::AL_2HE, Monomer::AL_CZ, Monomer::AL_HZ,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_CD1, Monomer::AL_CG,
		   Monomer::AL_CD2, Monomer::AL_CG,
		   // was HD1, HD2
		   Monomer::AL_1HD, Monomer::AL_CD1,
		   Monomer::AL_2HD, Monomer::AL_CD2, 
		   Monomer::AL_CE1, Monomer::AL_CD1,
		   Monomer::AL_CE2, Monomer::AL_CD2,
		   // HE1,2
		   Monomer::AL_1HE, Monomer::AL_CE1, 
		   Monomer::AL_2HE, Monomer::AL_CE2,
		   Monomer::AL_CZ, Monomer::AL_CE1, 
		   Monomer::AL_CZ, Monomer::AL_CE2, Monomer::AL_HZ,
	  Monomer::AL_CZ,Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CZ, Monomer::AL_INVALID}END_EXTREMES; //
      END_RES;
	
      BEGIN_RES(Monomer::NH2);
      BEGIN_ATOMS{Monomer::AL_N, Monomer::AL_1H, Monomer::AL_2H,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{ Monomer::AL_1H, Monomer::AL_N, 
		    Monomer::AL_2H, Monomer::AL_N,
	  Monomer::AL_INVALID}END_BONDS;
      END_RES;
	
      BEGIN_RES(Monomer::MET);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, 
		   Monomer::AL_SD, Monomer::AL_CE, Monomer::AL_1HE,
		   Monomer::AL_2HE, Monomer::AL_3HE,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_1HG, Monomer::AL_CG,
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_SD, Monomer::AL_CG,
		   Monomer::AL_CE, Monomer::AL_SD, 
		   Monomer::AL_1HE, Monomer::AL_CE,
		   Monomer::AL_2HE, Monomer::AL_CE, 
		   Monomer::AL_3HE, Monomer::AL_CE,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CE, Monomer::AL_INVALID}END_EXTREMES; 
      END_RES;

      BEGIN_RES(Monomer::LYS);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, Monomer::AL_CD, Monomer::AL_1HD, Monomer::AL_2HD,
		   Monomer::AL_CE, Monomer::AL_1HE, Monomer::AL_2HE, Monomer::AL_NZ, Monomer::AL_1HZ, Monomer::AL_2HZ, Monomer::AL_3HZ,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB,
		   Monomer::AL_1HG, Monomer::AL_CG, 
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_CD, Monomer::AL_CG,
		   Monomer::AL_1HD,Monomer::AL_CD, 
		   Monomer::AL_2HD, Monomer::AL_CD,
		   Monomer::AL_CE, Monomer::AL_CD, 
		   Monomer::AL_1HE, Monomer::AL_CE, 
		   Monomer::AL_2HE, Monomer::AL_CE,
		   Monomer::AL_NZ, Monomer::AL_CE,
		   Monomer::AL_1HZ, Monomer::AL_NZ, 
		   Monomer::AL_2HZ, Monomer::AL_NZ,
		   Monomer::AL_3HZ, Monomer::AL_NZ,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_NZ, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;
	
      BEGIN_RES(Monomer::LEU);
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_HG, Monomer::AL_CD1, Monomer::AL_CD2, Monomer::AL_1HD1, Monomer::AL_2HD1,
		    Monomer::AL_3HD1, Monomer::AL_1HD2, Monomer::AL_2HD2, Monomer::AL_3HD2,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{ Monomer::AL_CB, Monomer::AL_CA,
		    Monomer::AL_1HB, Monomer::AL_CB,
		    Monomer::AL_2HB, Monomer::AL_CB,
		    Monomer::AL_CG, Monomer::AL_CB, 
		    Monomer::AL_HG, Monomer::AL_CG,
		    Monomer::AL_CD1, Monomer::AL_CG,
		    Monomer::AL_CD2, Monomer::AL_CG, 
		    Monomer::AL_1HD1, Monomer::AL_CD1,
		    Monomer::AL_2HD1, Monomer::AL_CD1,
		    Monomer::AL_3HD1, Monomer::AL_CD1,
		    Monomer::AL_1HD2, Monomer::AL_CD2,
		    Monomer::AL_2HD2, Monomer::AL_CD2,
		    Monomer::AL_3HD2, Monomer::AL_CD2,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CD1, Monomer::AL_CD2, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;

      BEGIN_RES(Monomer::ILE);
      // was 1HD
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_HB, Monomer::AL_CG1, Monomer::AL_CG2, Monomer::AL_1HG1, Monomer::AL_2HG1, Monomer::AL_CD, Monomer::AL_1HD, Monomer::AL_2HD,
		   Monomer::AL_3HD, Monomer::AL_1HG2, Monomer::AL_2HG2, Monomer::AL_3HG2,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_HB, Monomer::AL_CB,
		   Monomer::AL_CG1, Monomer::AL_CB,
		   Monomer::AL_CG2, Monomer::AL_CB,
		   Monomer::AL_1HG1, Monomer::AL_CG1,
		   Monomer::AL_2HG1, Monomer::AL_CG1,
		   Monomer::AL_CD, Monomer::AL_CG1, 
		   Monomer::AL_1HD, Monomer::AL_CD,
		   Monomer::AL_2HD, Monomer::AL_CD,
		   Monomer::AL_3HD, Monomer::AL_CD,
		   Monomer::AL_1HG2, Monomer::AL_CG2,
		   Monomer::AL_2HG2, Monomer::AL_CG2,
		   Monomer::AL_3HG2, Monomer::AL_CG2,
	  Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CD, Monomer::AL_INVALID}END_EXTREMES; // maybe add CG2
      END_RES;
	
      BEGIN_RES(Monomer::HIS);
      // was HD1, HD2 HE1
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_ND1, Monomer::AL_CD2, Monomer::AL_1HD, Monomer::AL_2HD, Monomer::AL_CE1,
		   Monomer::AL_NE2, Monomer::AL_1HE,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_ND1, Monomer::AL_CG,
		   Monomer::AL_CD2, Monomer::AL_CG,
		   Monomer::AL_1HD, Monomer::AL_ND1, 
		   Monomer::AL_2HD, Monomer::AL_CD2,
		   Monomer::AL_CE1, Monomer::AL_ND1,
		   Monomer::AL_CE1, Monomer::AL_NE2, 
		   Monomer::AL_1HE, Monomer::AL_CE1, 
		   Monomer::AL_NE2, Monomer::AL_CD2,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CE1, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;
	
	
      BEGIN_RES(Monomer::GLY);
      BEGIN_ATOMS{Monomer::AL_2HA,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_2HA, Monomer::AL_CA,Monomer::AL_INVALID}END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;	
	
      BEGIN_RES(Monomer::GLU);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, Monomer::AL_CD, Monomer::AL_OE1, Monomer::AL_OE2,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB,
		   Monomer::AL_1HG, Monomer::AL_CG, 
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_CD, Monomer::AL_CG, 
		   Monomer::AL_OE1, Monomer::AL_CD, 
		   Monomer::AL_OE2, Monomer::AL_CD,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CD, Monomer::AL_INVALID}END_EXTREMES; // maybe should be Monomer::AL_OE[12]
      END_RES;	
	
      BEGIN_RES(Monomer::GLN);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, Monomer::AL_CD, Monomer::AL_OE1, Monomer::AL_NE2,
	  Monomer::AL_1HE2, Monomer::AL_2HE2,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_1HG, Monomer::AL_CG,
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_CD, Monomer::AL_CG, 
		   Monomer::AL_OE1, Monomer::AL_CD, 
		   Monomer::AL_NE2, Monomer::AL_CD,
		   Monomer::AL_1HE2, Monomer::AL_NE2, 
		   Monomer::AL_2HE2, Monomer::AL_NE2,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CD, Monomer::AL_INVALID}END_EXTREMES; // or Monomer::AL_NE2, Monomer::AL_OE1
      END_RES;	
	
      BEGIN_RES(Monomer::CYS);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_SG, Monomer::AL_HG,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_SG, Monomer::AL_CB, 
		   Monomer::AL_HG, Monomer::AL_SG,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_SG, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;	
	
      BEGIN_RES(Monomer::ASP);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_OD1, Monomer::AL_OD2,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA,
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_OD1, Monomer::AL_CG, 
		   Monomer::AL_OD2, Monomer::AL_CG,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CG, Monomer::AL_INVALID}END_EXTREMES; // or Monomer::AL_OD[12]
      END_RES;	
	
      BEGIN_RES(Monomer::ASN);
      BEGIN_ATOMS{ Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_OD1, Monomer::AL_ND2, Monomer::AL_1HD2, Monomer::AL_2HD2,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{ Monomer::AL_CB, Monomer::AL_CA,
		    Monomer::AL_1HB, Monomer::AL_CB, 
		    Monomer::AL_2HB, Monomer::AL_CB,
		    Monomer::AL_CG, Monomer::AL_CB,
		    Monomer::AL_OD1, Monomer::AL_CG,
		    Monomer::AL_ND2, Monomer::AL_CG,
		    Monomer::AL_1HD2, Monomer::AL_ND2, 
		    Monomer::AL_2HD2, Monomer::AL_ND2 ,
		    Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CG, Monomer::AL_INVALID}END_EXTREMES; // or Monomer::AL_OD1, Monomer::AL_ND1
      END_RES;	
	
      BEGIN_RES(Monomer::ARG);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_CG, Monomer::AL_1HG, Monomer::AL_2HG, Monomer::AL_CD, Monomer::AL_1HD, Monomer::AL_2HD,
		   Monomer::AL_NE, Monomer::AL_HE, Monomer::AL_CZ, Monomer::AL_NH1, Monomer::AL_NH2, Monomer::AL_1HH1, Monomer::AL_2HH1, Monomer::AL_1HH2,
	  Monomer::AL_2HH2,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CB, Monomer::AL_CA, 
		   Monomer::AL_1HB, Monomer::AL_CB, 
		   Monomer::AL_2HB, Monomer::AL_CB,
		   Monomer::AL_CG, Monomer::AL_CB, 
		   Monomer::AL_1HG, Monomer::AL_CG, 
		   Monomer::AL_2HG, Monomer::AL_CG,
		   Monomer::AL_CD, Monomer::AL_CG,
		   Monomer::AL_1HD, Monomer::AL_CD,
		   Monomer::AL_2HD, Monomer::AL_CD,
		   Monomer::AL_NE, Monomer::AL_CD, 
		   Monomer::AL_HE, Monomer::AL_NE, 
		   Monomer::AL_CZ, Monomer::AL_NE,
		   Monomer::AL_NH1, Monomer::AL_CZ, 
		   Monomer::AL_NH2, Monomer::AL_CZ, 
		   Monomer::AL_1HH1, Monomer::AL_NH1,
		   Monomer::AL_2HH1, Monomer::AL_NH1, 
		   Monomer::AL_1HH2, Monomer::AL_NH2, 
		   Monomer::AL_2HH2, Monomer::AL_NH2,
		   Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CZ, Monomer::AL_INVALID}END_EXTREMES; // maybe Monomer::AL_NH2m Monomer::AL_NH1
      END_RES;	
	
      BEGIN_RES(Monomer::ALA);
      BEGIN_ATOMS{Monomer::AL_CB, Monomer::AL_1HB, Monomer::AL_2HB, Monomer::AL_3HB,Monomer::AL_INVALID} END_ATOMS;
      BEGIN_BONDS{ Monomer::AL_CB, Monomer::AL_CA,
		    Monomer::AL_1HB, Monomer::AL_CB, 
		    Monomer::AL_2HB, Monomer::AL_CB,
		    Monomer::AL_3HB, Monomer::AL_CB,
		    Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_CB, Monomer::AL_INVALID}END_EXTREMES;
      END_RES;	
	
      BEGIN_RES(Monomer::ACE);
      BEGIN_ATOMS{Monomer::AL_CH3, Monomer::AL_C, Monomer::AL_O, Monomer::AL_1HH3, Monomer::AL_2HH3, Monomer::AL_3HH3,Monomer::AL_INVALID}END_ATOMS;
      BEGIN_BONDS{Monomer::AL_CH3, Monomer::AL_C, 
	  Monomer::AL_O, Monomer::AL_C,
	  Monomer::AL_C, Monomer::AL_1HH3,
	  Monomer::AL_C, Monomer::AL_2HH3, 
	  Monomer::AL_C, Monomer::AL_3HH3,
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;

      END_RES;	
      BEGIN_RES(Monomer::ADE);
      BEGIN_ATOMS{Monomer::AL_N9, Monomer::AL_C8, Monomer::AL_H8,
	  Monomer::AL_N7, Monomer::AL_C5, Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_C2, Monomer::AL_H2, Monomer::AL_N1, Monomer::AL_C6, Monomer::AL_N6,
	  Monomer::AL_H61, Monomer::AL_H62} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_N9, Monomer::AL_C1p,
	  Monomer::AL_N9, Monomer::AL_C8,
	  Monomer::AL_C8, Monomer::AL_H8,
	  Monomer::AL_C8, Monomer::AL_N7,
	  Monomer::AL_N7, Monomer::AL_C5,
	  Monomer::AL_C4, Monomer::AL_C5,
	  Monomer::AL_C4, Monomer::AL_N9,
	  Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_N3, Monomer::AL_C2,
	  Monomer::AL_H2, Monomer::AL_C2,
	  Monomer::AL_C2, Monomer::AL_N1,
	  Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_C6, Monomer::AL_N6,
	  Monomer::AL_N6, Monomer::AL_H61,
	  Monomer::AL_N6, Monomer::AL_H62,
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;

      BEGIN_RES(Monomer::GUA);
      BEGIN_ATOMS{Monomer::AL_N9, Monomer::AL_C8, Monomer::AL_H8,
	  Monomer::AL_N7, Monomer::AL_C5, Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_C2, Monomer::AL_H2, Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_O6, Monomer::AL_H1, Monomer::AL_N2,
	  Monomer::AL_H21, Monomer::AL_H22} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_N9, Monomer::AL_C1p,
	  Monomer::AL_N9, Monomer::AL_C8,
	  Monomer::AL_C8, Monomer::AL_H8,
	  Monomer::AL_C8, Monomer::AL_N7,
	  Monomer::AL_N7, Monomer::AL_C5,
	  Monomer::AL_C4, Monomer::AL_C5,
	  Monomer::AL_C4, Monomer::AL_N9,
	  Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_N3, Monomer::AL_C2,
	  Monomer::AL_C2, Monomer::AL_N1,
	  Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_C6, Monomer::AL_O6,
	  Monomer::AL_N1, Monomer::AL_H1,
	  Monomer::AL_N2, Monomer::AL_C2,
	  Monomer::AL_N2, Monomer::AL_H21,  
	  Monomer::AL_N2, Monomer::AL_H22,  
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;


      BEGIN_RES(Monomer::CYT);
      BEGIN_ATOMS{Monomer::AL_N1, Monomer::AL_C6, Monomer::AL_H6, Monomer::AL_C5, Monomer::AL_H5, 
	  Monomer::AL_C4, Monomer::AL_N3, Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_N4, Monomer::AL_H41, Monomer::AL_H42} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_C1p, Monomer::AL_N1,
	  Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_C6, Monomer::AL_H6,
	  Monomer::AL_C6, Monomer::AL_C5,
	  Monomer::AL_C5, Monomer::AL_H5,
	  Monomer::AL_C5, Monomer::AL_C4,
	  Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_N3, Monomer::AL_C2,
	  Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_C2, Monomer::AL_N1,
	  Monomer::AL_C4, Monomer::AL_N4,
	  Monomer::AL_H41, Monomer::AL_N4,
	  Monomer::AL_N4, Monomer::AL_H42,
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;

      BEGIN_RES(Monomer::URA);
      BEGIN_ATOMS{Monomer::AL_N1, Monomer::AL_C6, Monomer::AL_H6, Monomer::AL_C5, Monomer::AL_H5, 
	  Monomer::AL_C4, Monomer::AL_N3, Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_H3, Monomer::AL_O4} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_C1p, Monomer::AL_N1,
	  Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_C6, Monomer::AL_H6,
	  Monomer::AL_C6, Monomer::AL_C5,
	  Monomer::AL_C5, Monomer::AL_H5,
	  Monomer::AL_C5, Monomer::AL_C4,
	  Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_N3, Monomer::AL_C2,
	  Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_C2, Monomer::AL_N1,
	  Monomer::AL_C4, Monomer::AL_O4,
	  Monomer::AL_N3, Monomer::AL_H3,
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;

      BEGIN_RES(Monomer::THY);
      BEGIN_ATOMS{Monomer::AL_N1, Monomer::AL_C6, Monomer::AL_H6, Monomer::AL_C5, Monomer::AL_C7,
	  Monomer::AL_H71, Monomer::AL_H72, Monomer::AL_H73, 
	  Monomer::AL_C4, Monomer::AL_N3, Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_H3, Monomer::AL_O4} END_ATOMS;
      BEGIN_BONDS{Monomer::AL_C1p, Monomer::AL_N1,
	  Monomer::AL_N1, Monomer::AL_C6,
	  Monomer::AL_C6, Monomer::AL_H6,
	  Monomer::AL_C6, Monomer::AL_C5,
	  Monomer::AL_C5, Monomer::AL_C7,
	  Monomer::AL_C7, Monomer::AL_H71,
	  Monomer::AL_C7, Monomer::AL_H72,
	  Monomer::AL_C7, Monomer::AL_H73,
	  Monomer::AL_C5, Monomer::AL_C4,
	  Monomer::AL_C4, Monomer::AL_N3,
	  Monomer::AL_N3, Monomer::AL_C2,
	  Monomer::AL_C2, Monomer::AL_O2,
	  Monomer::AL_C2, Monomer::AL_N1,
	  Monomer::AL_C4, Monomer::AL_O4,
	  Monomer::AL_N3, Monomer::AL_H3,
	  Monomer::AL_INVALID} END_BONDS;
      BEGIN_EXTREMES{Monomer::AL_INVALID}END_EXTREMES;
      END_RES;

    default:
      ;
    }
    
    return dat;
    //return std::pair<Monomer::Atom_key*, Monomer::Atom_key*>(als, bls);
  }

  void do_initialize() {
    //if (amino_acid_initialized_) return;
    assert(amino_acid_initialized_==false);
    amino_acid_initialized_=true;
    unsigned int i=0;
    for (; atom_name_data_[i].l != Monomer::AL_INVALID; ++i){
      clean_atom_name_data_[i].l= atom_name_data_[i].l;
      clean_atom_name_data_[i].t= atom_name_data_[i].t;
      std::sscanf(atom_name_data_[i].s, "%4s", clean_atom_name_data_[i].s);
    }
    clean_atom_name_data_[i].l= Monomer::AL_INVALID;

    /*unsigned int max_atoms=0;
    for (unsigned int i=0;atom_fallback_data_[i].l != Monomer::AL_INVALID; ++i){
      max_atoms= (std::max)(max_atoms, static_cast<unsigned int>(atom_fallback_data_[i].l));
      }*/
    clean_atom_fallbacks_.resize(Monomer::AL_LAST_LABEL+1);

    for (unsigned int i=0;atom_fallback_data_[i].l != Monomer::AL_INVALID; ++i){
      clean_atom_fallbacks_[atom_fallback_data_[i].l].push_back(atom_fallback_data_[i].lr);
    }


    const unsigned int num_res= static_cast<int>(Monomer::INV)+1;
    //Monomer::Type all_res[num_res];
    /*for(unsigned int i=0; i< num_res; ++i) {
      all_res[i]=Monomer::Type(i);
      }*/
    //sizeof(all_res)/sizeof(Monomer::Type);


    amino_acid_data_.resize(num_res);
    //http://www.bmrb.wisc.edu/referenc/nomenclature/

    for (unsigned int i=0; i< num_res; ++i){
      Monomer::Type cur_res=Monomer::Type(i);
      amino_acid_data_[cur_res]= Amino_acid_data();
      if (cur_res == Monomer::ADE || cur_res == Monomer::CYT
	  || cur_res == Monomer::GUA || cur_res == Monomer::URA || cur_res == Monomer::THY){
	Monomer::Atom_key bl[]={Monomer::AL_P, Monomer::AL_OP1, Monomer::AL_OP2,
				Monomer::AL_O5p, Monomer::AL_C5p, Monomer::AL_H5p,
				Monomer::AL_H5pp,
				Monomer::AL_C4p, Monomer::AL_H4p, Monomer::AL_O4p,
				Monomer::AL_C1p, Monomer::AL_H1p,
				Monomer::AL_C3p, Monomer::AL_H3p, Monomer::AL_O3p,
				Monomer::AL_C2p, Monomer::AL_H2p,
				Monomer::AL_H2pp, Monomer::AL_O2p, Monomer::AL_HO2p};
	for (unsigned int j=0; j< sizeof(bl)/sizeof(Monomer::Atom_key); ++j){
	  amino_acid_data_[cur_res].atoms.push_back(bl[j]);
	}
	Monomer::Atom_key bd[][2]={ {Monomer::AL_P, Monomer::AL_OP2}, {Monomer::AL_P, Monomer::AL_OP1},
				    {Monomer::AL_P, Monomer::AL_O5p}, {Monomer::AL_O5p, Monomer::AL_C5p},
				    {Monomer::AL_C5p, Monomer::AL_H5p}, {Monomer::AL_C5p, Monomer::AL_H5pp},
				    {Monomer::AL_C5p, Monomer::AL_C4p}, {Monomer::AL_C4p, Monomer::AL_H4p},
				    {Monomer::AL_C4p, Monomer::AL_C3p}, {Monomer::AL_C3p, Monomer::AL_H3p}, 
				    {Monomer::AL_C3p, Monomer::AL_O3p}, {Monomer::AL_C3p, Monomer::AL_C2p},
				    {Monomer::AL_C2p, Monomer::AL_H2pp}, {Monomer::AL_C2p, Monomer::AL_HO2p},
				    {Monomer::AL_C2p, Monomer::AL_C1p}, {Monomer::AL_C1p, Monomer::AL_H1p},
				    {Monomer::AL_C1p, Monomer::AL_O4p}, {Monomer::AL_O4p, Monomer::AL_C4p}};
	for (unsigned int j=0; j< sizeof(bd)/sizeof(Possible_bond); ++j){
	  amino_acid_data_[cur_res].bonds.push_back(Possible_bond(bd[j][0], bd[j][1]));
	}
      } else if (cur_res != Monomer::ACE && cur_res != Monomer::NH2) {
	// set up shared atoms
	Monomer::Atom_key bl[]={ Monomer::AL_N, Monomer::AL_H,
				   Monomer::AL_1H, Monomer::AL_2H, 
				   Monomer::AL_3H, Monomer::AL_CA, 
				   Monomer::AL_HA, Monomer::AL_1HA,
				   Monomer::AL_C, Monomer::AL_O,
				   Monomer::AL_OXT};
	for (unsigned int j=0; j< sizeof(bl)/sizeof(Monomer::Atom_key); ++j){
	  amino_acid_data_[cur_res].atoms.push_back(bl[j]);
	}
	// set up basic bonds
	Monomer::Atom_key bd[][2]={ {Monomer::AL_N, Monomer::AL_H}, {Monomer::AL_N, Monomer::AL_1H},
				      {Monomer::AL_N, Monomer::AL_2H}, {Monomer::AL_N, Monomer::AL_3H},
				      {Monomer::AL_N, Monomer::AL_CA}, {Monomer::AL_CA, Monomer::AL_HA},
				      {Monomer::AL_CA, Monomer::AL_1HA}, {Monomer::AL_CA, Monomer::AL_C},
				      {Monomer::AL_C, Monomer::AL_O}, {Monomer::AL_O, Monomer::AL_OXT} };
	for (unsigned int j=0; j< sizeof(bd)/sizeof(Possible_bond); ++j){
	  amino_acid_data_[cur_res].bonds.push_back(Possible_bond(bd[j][0], bd[j][1]));
	}
      }

      /*
	(Monomer::Monomer::AL_CB, Monomer::Monomer::AL_HB, Monomer::Monomer::AL_CG1,  Monomer::Monomer::AL_CG2, 
	Monomer::Monomer::AL_1HG1, Monomer::Monomer::AL_2HG1, Monomer::Monomer::AL_3HG1, Monomer::Monomer::AL_1HG2,
	Monomer::Monomer::AL_2HG2, Monomer::Monomer::AL_3HG2, Monomer::Monomer::AL_INVALID) ,
	()
      */

      Monomer_init_data data= get_residue_initialization_data(cur_res);
      
      
      for (unsigned int j=0; data.atms_[j] != Monomer::AL_INVALID; ++j){
	amino_acid_data_[cur_res].atoms.push_back(data.atms_[j]);
      }
      for (unsigned int j=0; data.bnds_[j] != Monomer::AL_INVALID; j+=2){
	amino_acid_data_[cur_res].bonds.push_back(Possible_bond(data.bnds_[j], data.bnds_[j+1]));
      }
      for (unsigned int j=0; data.extr_[j] != Monomer::AL_INVALID; ++j){
	amino_acid_data_[cur_res].extremes.push_back(data.extr_[j]);
      }

      if (false){
	std::cout << "For residue " << Monomer::type_string(cur_res)
		  << " there are " << amino_acid_data_[cur_res].atoms.size() << " atoms, "
		  << amino_acid_data_[cur_res].bonds.size() << " bonds, and "
		  << amino_acid_data_[cur_res].extremes.size() << " extremes." << std::endl;
	assert(amino_acid_data_[cur_res].extremes.size() < 4);
      }
      
    }
  }


  Monomer::Atom_key fix_atom_label(Monomer::Type label, Monomer::Atom_key al) {
    for (unsigned int i=0; i< amino_acid_data_[label].atoms.size(); ++i){
      if (Monomer_data::amino_acid_data_[label].atoms[i] == al) return al;
    }
    //HERE
    /*if (clean_atom_fallbacks_.find(al) == clean_atom_fallbacks_.end()){
      dsrpdb_internal::error_logger.new_warning((std::string("No atom named ")+atom_key_string(al) 
      +" in residue "+ type_string(label_)).c_str());
      return AL_INVALID;
      }*/
    for (unsigned int j=0; j< clean_atom_fallbacks_[al].size(); ++j){
      for (unsigned int i=0; i< amino_acid_data_[label].atoms.size(); ++i){
	if (amino_acid_data_[label].atoms[i] == clean_atom_fallbacks_[al][j]) 
	  return clean_atom_fallbacks_[al][j];
      }
    }
    return Monomer::AL_INVALID;
  }


}

}}
