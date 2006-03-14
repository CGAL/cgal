#ifndef OGL_VISU
#define OGL_VISU

#include <math.h>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/basic.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "enriched_polyhedron.h" 

#define RIDGES_CALL_LIST 1

typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 CGAL_point;


//Mesh
//------------------
typedef float	                   FT;
typedef CGAL::Simple_cartesian<FT> K0;
typedef CGAL::Filtered_kernel <K0> kernel;

// CGAL polyhedron
typedef Enriched_polyhedron<kernel,Enriched_items> Mesh;

//Data structure for one line of the input file .4ogl.txt
struct data_line{
  int ridge_type;
  double strength, sharpness;
  std::list<CGAL_point> ridge_points;
  data_line(int ridge_type, double strength, double sharpness,
	    std::list<CGAL_point> ridge_points): 
  ridge_type(ridge_type), strength(strength), sharpness(sharpness),
    ridge_points(ridge_points)
    {};
};


typedef std::list<data_line* > DS;
typedef std::list<data_line* >::iterator DS_iterator;


//VARIABLES
//---------------

extern const char* file_off;
extern const char* file_res;
  
extern char presse;
extern int anglex,angley,x,y,xold,yold,rpresse, xrold,U_scale;
extern Mesh m_mesh;

//General ogl visu
void affichage();
void clavier(unsigned char touche,int x,int y);
void reshape(int x,int y);
void idle();
void mouse(int bouton,int etat,int x,int y);
void mousemotion(int x,int y);


//specific ridge visu
void clean_DS(DS& L);
void read_line(std::ifstream& stream_res, int& ridge_type, 
	       double& strength, double& sharpness, 
	       std::list<CGAL_point>& ridge_points);
void load_data_from_file(DS& l);
void draw_one_ridge(data_line* line);
void MakeCallList(DS& L);

void run_visu(int argc, char* argv[]);

#endif
