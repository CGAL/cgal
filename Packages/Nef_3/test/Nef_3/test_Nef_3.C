// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : demo/Nef_3/nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Demo program maintaining a stack of Nef polyhedra in the space and
// a manipulation language for stack ops, file loading and saving, etc.
// ============================================================================

// set this macro if you have the OpenGL and glut based visualization
#define CGAL_NEF3_VISUALIZOR

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>

template<typename Kernel>
class test {

  typedef typename Kernel::Point_3                       Point_3;
  typedef typename Kernel::Plane_3                       Plane_3;
  typedef typename Kernel::Vector_3                      Vector_3;
  typedef typename Kernel::Aff_transformation_3          Aff_transformation_3;
  typedef CGAL::Polyhedron_3<Kernel>                     Polyhedron;
  typedef CGAL::SNC_items<Kernel, bool>                  SNC_items;
  typedef CGAL::SNC_structure<SNC_items>                 SNC_structure;
  typedef typename SNC_structure::Vertex_const_iterator  Vertex_const_iterator;
  typedef CGAL::Nef_polyhedron_3<SNC_items>              Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_explorer          SNC_explorer;
  typedef typename Nef_polyhedron::SM_explorer           SM_explorer;

public:
  test() {};

private:
  bool are_files_equal(char* name1, char* name2) {
    std::ifstream in1(name1);
    std::ifstream in2(name2);
    std::string s1;
    std::string s2;
    while(in1) {
      in1 >> s1;
      in2 >> s2;
      if(s1 != s2) {
	cerr << s1 << std::endl;
	cerr << s2 << std::endl;
	return false;
      }
    }
    if(in2)
      return false;
    return true;
  }

  bool does_nef3_equals_file(Nef_polyhedron& N, char* name, char* suffix) {
    char* fullname = new char[strlen(name)+strlen(suffix)+1];
    strcpy(fullname, name);
    strcat(fullname, suffix);
    std::ofstream out("temp.nef3");
    N.dump(out);
    return are_files_equal("temp.nef3",fullname);
    delete [] fullname;
  }

  Nef_polyhedron load_off( char* name) {
    Polyhedron poly;
    std::ifstream in(name);
    in >> poly;
    Nef_polyhedron N(poly);
    return N;
  }

  Nef_polyhedron load_nef3(char* name, char* suffix) {
    char* fullname = new char[strlen(name)+strlen(suffix)+1];
    strcpy(fullname, name);
    strcat(fullname, suffix);
    return Nef_polyhedron(fullname);
  }

public:
  void make_test(char* suffix) {
    
    Nef_polyhedron C;
    Nef_polyhedron N;
    Nef_polyhedron N1;
    Nef_polyhedron N2;
    Nef_polyhedron N3;

    C = load_off("cube.off");
    N = C;
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3",suffix));

    // ********************* Construction *****************************

    if(suffix[1] == 'E') {
   
      N = Nef_polyhedron(Plane_3(3,4,5,0)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-0.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(3,4,5,31));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-31.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-0.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,29)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-29.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(1,2,3,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(1,2,3,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-23.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(1,2,4,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-0.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(1,2,4,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-23.nef3",suffix)); 
    }      

    // ******************** marks_of_halfspheres *******************

    if(suffix[1] == 'E') {   
      N = load_nef3("marks_of_halfspheres.nef3",suffix);
      int i=0;
      SNC_explorer E(N.SNCexplorer());
      Vertex_const_iterator vi = E.vertices_begin();
    
      do {vi++;} while(++i < 8); //  -1 1 0
      CGAL_assertion(E.point(vi) == Point_3(-1,1,0));
      SM_explorer SME(N.SMexplorer(vi));
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 9); //  0 1 -1
      CGAL_assertion(E.point(vi) == Point_3(0,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 10); //  1 1 0
      CGAL_assertion(E.point(vi) == Point_3(1,1,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 11); //  0 1 1
      CGAL_assertion(E.point(vi) == Point_3(0,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 12); //  0 0 -1
      CGAL_assertion(E.point(vi) == Point_3(0,0,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 13); //  0 0 1
      CGAL_assertion(E.point(vi) == Point_3(0,0,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 14); //  -1 0 0
      CGAL_assertion(E.point(vi) == Point_3(-1,0,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 15); //  1 0 0
      CGAL_assertion(E.point(vi) == Point_3(1,0,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 17); //  -1 1 1
      CGAL_assertion(E.point(vi) == Point_3(-1,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 18); //  1 1 1
      CGAL_assertion(E.point(vi) == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 21); //  -1 1 -1
      CGAL_assertion(E.point(vi) == Point_3(-1,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 22); //  1 1 -1
      CGAL_assertion(E.point(vi) == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
    }

    // ******************** Simplification *************************

    N = load_off("cube+v.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3",suffix));

    N = load_off("cube+vee.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube1.nef3",suffix));

    N = load_off("cube+veeee.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube2.nef3",suffix)); 

    N = load_off("cube+vONe.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube3.nef3",suffix));

    N = load_off("cube.off");
    N1 = N;
    N2 = N;
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,0,0,1)));
    N = N.symmetric_difference(N1);
    CGAL_assertion(N.is_valid(0,0));
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(0,2,2,1)));
    N = N.join(N2);
    CGAL_assertion(N.is_valid(0,0));
    N = N.regularization();
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"simplifySface.nef3",suffix));    

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,1,0,1)));
    CGAL_assertion(N2.is_valid(0,0));
    N1=N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2,-1,4,2)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+plane.nef3",suffix));
 
    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,0,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,2,1)));
    N1 = N.join(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+line.nef3",suffix));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,2,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,0,1)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+vertex1.nef3",suffix));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,2,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,-1,1)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+vertex2.nef3",suffix));

    if(suffix[1] == 'E') {   
      N = C;
      CGAL_assertion(N.is_valid(0,0));
      N1 = N;
      N2 = N;
      N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,0,0,1)));
      CGAL_assertion(N2.is_valid(0,0));
      N1 = N1.intersection(N2);
      CGAL_assertion(N1.is_valid(0,0));
      N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,0,0,1)));
      CGAL_assertion(N1.is_valid(0,0));
      N2 = Nef_polyhedron("viereck.nef3");
      CGAL_assertion(N2.is_valid(0,0));
      N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,3,0,1)));
      CGAL_assertion(N2.is_valid(0,0));
      N2 = N1.difference(N2);
      CGAL_assertion(N2.is_valid(0,0));
      N = N.difference(N1);
      CGAL_assertion(N.is_valid(0,0));
      N = N.join(N2);
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"donotsimplify.nef3",suffix));
    }      

    // ******************* Union Find Test **********************************

    if(suffix[1] == 'E') {   
      N = load_nef3("unionfind.nef3", suffix);
      CGAL_assertion(N.is_valid(0,0));
      N1 =  N;
      N.transform(Aff_transformation_3(1, 0, 0, 
				       0, 0,-1,
				       0, 1, 0,
				       1));
      CGAL_assertion(N1.is_valid(0,0));
      N1.transform(Aff_transformation_3(1, 0, 0, 
					0, 0,-1,
					0, 1, 0,
					1));
      CGAL_assertion(N1.is_valid(0,0));

      N2 = N.symmetric_difference(N1);
      CGAL_assertion(N2.is_valid(0,0));
      N2 = N2.regularization();
      CGAL_assertion(N2.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(N2,"unionfindref.nef3",suffix));

      N2 = N.join(N1);
      CGAL_assertion(N2.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N2,"unionfindref.nef3",suffix));
    }      

    // ************ Pluecker coordinates + binary operation of halfspaces **********

    if(suffix[1] == 'E') {
      N = Nef_polyhedron(Nef_polyhedron::COMPLETE);  
      CGAL_assertion(N.is_valid(0,0));               
      N1 = Nef_polyhedron(Plane_3(1,0,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(-1,0,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,1,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,-1,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,0,1,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,0,-1,-1));
      N = N.intersection(N1);      
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube_created_from_halfspaces.nef3", suffix));
    }

    // ******************* Unary Operations ********************************
    
    if(suffix[1] == 'E') {
      Nef_polyhedron T("topologyE.nef3");
      N = T;
      N = N.boundary();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"boundary.nef3",suffix));
    
      N = T;
      N = N.interior();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"interior.nef3",suffix));
      
      N = T;
      N = N.closure();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"closure.nef3",suffix));
      
      N = T;
      N = N.regularization();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"regularization.nef3",suffix));
      
      N = T; 
      N = N.complement();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"complement.nef3",suffix));
    }

    // ******************* Binary Operations ********************************


  }

};

int main() {

  typedef CGAL::Gmpz                         NT;
  typedef CGAL::Simple_homogeneous<NT>       SH_Kernel;
  typedef CGAL::Extended_homogeneous_3<NT>   EH_Kernel;
  
  test<SH_Kernel> test_SH;
  test<EH_Kernel> test_EH;

  test_SH.make_test(".SH");
  test_EH.make_test(".EH");
}

