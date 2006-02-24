#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <list>

#include "../../include/CGAL/Ridges.h" 
#include "../../../Jet_fitting_3/include/CGAL/Monge_via_jet_fitting.h" 
#include "GSL.h" 
 
#include "PolyhedralSurf.h"
#include "PolyhedralSurf_rings.h"
#include "options.h"//parsing command line

//Kernel of the PolyhedralSurf
typedef double                DFT;
typedef CGAL::Cartesian<DFT>  Data_Kernel;
typedef Data_Kernel::Point_3  DPoint;
typedef Data_Kernel::Vector_3 DVector;
typedef PolyhedralSurf::Vertex Vertex;
typedef PolyhedralSurf::Vertex_iterator Vertex_iterator;
typedef PolyhedralSurf::Vertex Vertex;

typedef T_PolyhedralSurf_rings<PolyhedralSurf> Poly_rings;
//Kernel for local computations
typedef double                LFT;
typedef CGAL::Cartesian<LFT>  Local_Kernel;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel, Local_Kernel, GSL> My_Monge_via_jet_fitting;
typedef CGAL::Monge_rep<Data_Kernel> My_Monge_rep;
typedef CGAL::Monge_info<Local_Kernel> My_Monge_info;
      
//RIDGES
typedef CGAL::Ridge_line<PolyhedralSurf> Ridge_line;
// typedef CGAL::Ridge_approximation<PolyhedralSurf,
//  std::vector<Ridge_line*>::iterator > Ridge_approximation;
  
//Syntax requirred by Options
static const char *const optv[] = {
  "?|?",
  "f:fName <string>",	//name of the input off file
  "d:deg <int>",	//degree of the jet
  "m:mdegree <int>",	//degree of the Monge rep
  "a:nrings <int>",	//# rings
  "p:npoints <int>",	//# points
  "v|",//verbose?
  NULL
};
 
// default fct parameter values and global variables
unsigned int d_fitting = 3;
unsigned int d_monge = 3;
unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
unsigned int nb_points_to_use = 0;//
bool verbose = false;
unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

//gather points around the vertex v using rings on the polyhedralsurf
//for the fitting
void gather_fitting_points( Vertex* v, 
			    std::vector<DPoint> &in_points)
{
  //container to collect vertices of v on the PolyhedralSurf
  std::vector<Vertex*> current_ring, next_ring, gathered; 
  std::vector<Vertex*> *p_current_ring, *p_next_ring;

  //initialize
  unsigned int nbp = 0, //current nb of collected points
    ith = 0;	//i-th ring index
  current_ring.clear();
  next_ring.clear();
  p_current_ring = &current_ring;
  p_next_ring = &next_ring;
  gathered.clear();
  in_points.clear();  
  
  //DO NOT FORGET TO UNTAG AT THE END!
  v->setRingIndex(ith);
  //collect 0th ring : the vertex v!
  gathered.push_back(v);
  nbp = 1;
  //collect 1-ring
  ith = 1;
  nbp = Poly_rings::push_neighbours_of(v, ith, current_ring, gathered);
  //collect more neighbors depending on options...
 
    //OPTION -p nb, with nb != 0
    //for approximation/interpolation with a fixed nb of points, collect
    // enough rings and discard some points of the last collected ring 
    // to get the exact "nb_points_to_use"
  if ( nb_points_to_use != 0 ) {
    while( gathered.size() < nb_points_to_use ) {
      ith++;
      //using tags
      nbp += Poly_rings::
	collect_ith_ring_neighbours(ith, *p_current_ring,
				    *p_next_ring, gathered);
      //next round must be launched from p_nextRing...
      p_current_ring->clear();
      std::swap(p_current_ring, p_next_ring);
    }
    //clean up
    Poly_rings::reset_ring_indices(gathered);
    //discard non-required collected points of the last ring
    gathered.resize(nb_points_to_use, NULL);
    assert(gathered.size() == nb_points_to_use );
  }
  else{
    //OPTION -a nb, with nb = 0
    //  select the mini nb of rings needed to make approx possible
    if (nb_rings == 0) {
      while ( gathered.size() < min_nb_points ) {
	ith++;
	nbp += Poly_rings::
	  collect_ith_ring_neighbours(ith, *p_current_ring,
				      *p_next_ring, gathered);
	//next round must be launched from p_nextRing...
	p_current_ring->clear();
	std::swap(p_current_ring, p_next_ring);
      }
    }
    //OPTION -a nb, with nb = 1, nothing to do! we have already
    //      collected the 1 ring
    //OPTION -a nb, with nb > 1
    //for approximation with a fixed nb of rings, collect
    //      neighbors up to the "nb_rings"th ring
    if (nb_rings > 1)
      while (ith < nb_rings) {
	ith++;
	nbp += Poly_rings::
	  collect_ith_ring_neighbours(ith, *p_current_ring,
				      *p_next_ring, gathered);
	//next round must be launched from p_nextRing...
	p_current_ring->clear();
	std::swap(p_current_ring, p_next_ring);
      }
   
    //clean up
    Poly_rings::reset_ring_indices(gathered);
  } //END ELSE

  //store the gathered points
  std::vector<Vertex*>::iterator itb = gathered.begin(),
    ite = gathered.end();
  CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
  
}


void compute_differential_quantities(PolyhedralSurf& P)
{
  //container for approximation points
  std::vector<DPoint> in_points;
 
  //MAIN LOOP
  Vertex_iterator vitb = P.vertices_begin(), vite = P.vertices_end();
  for (; vitb != vite; vitb++) {
    //initialize
    Vertex* v = &(*vitb);
    in_points.clear();  
    My_Monge_rep monge_rep;
    My_Monge_info monge_info;
      
    //gather points arourd the vertex using rings
    gather_fitting_points( v, in_points);

    //exit if the nb of points is to small 
    if ( in_points.size() < min_nb_points ) exit(0);

    //For Ridges we need at least 3rd order info
    assert( d_monge >= 3);
    // run the main fct : perform the fitting
    My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				   d_fitting, d_monge, 
				   monge_rep, monge_info);
    
    //switch min-max ppal curv/dir wrt the mesh orientation
    const DVector normal_mesh = P.computeFacetsAverageUnitNormal(v);
    monge_rep.comply_wrt_given_normal(normal_mesh);
       
    //    Store monge data needed for ridge computations in v
    v->d1() = monge_rep.d1();
    v->d2() = monge_rep.d2();
    v->k1() = monge_rep.coefficients()[0];
    v->k2() = monge_rep.coefficients()[1];
    v->b0() = monge_rep.coefficients()[2];
    v->b3() = monge_rep.coefficients()[5];
    if ( d_monge >= 4) {
      //= 3*b1^2+(k1-k2)(c0-3k1^3)
      v->P1() =
	3*monge_rep.coefficients()[3]*monge_rep.coefficients()[3]
	+(v->k1()-v->k2())
	*(monge_rep.coefficients()[6]-3*v->k1()*v->k1()*v->k1()) ; 
      //= 3*b2^2+(k2-k1)(c4-3k2^3)
      v->P2() = 
	3*monge_rep.coefficients()[4]*monge_rep.coefficients()[4]
	+(-v->k1()+v->k2())
	*(monge_rep.coefficients()[10]-3*v->k2()*v->k2()*v->k2()) ; 
    }
  } //END FOR LOOP
}


int main(int argc, char *argv[])
{  
  char *if_name = NULL, //input file name
    *w_if_name = NULL;  //as above, but / replaced by _
  char* res4openGL_fname;
  char* verbose_fname;
  std::ofstream *out_4ogl = NULL, *out_verbose = NULL;

  int optchar;
  char *optarg;
  Options opts(*argv, optv);
  OptArgvIter iter(--argc, ++argv);
 
  while ((optchar = opts(iter, (const char *&) optarg))){
    switch (optchar){
    case 'f': if_name = optarg; break;
    case 'd': d_fitting = atoi(optarg); break;
    case 'm': d_monge = atoi(optarg); break;
    case 'a': nb_rings = atoi(optarg); break;
    case 'p': nb_points_to_use = atoi(optarg); break;
    case 'v': verbose=true; break;
    default:
      cerr << "Unknown command line option " << optarg;
      exit(0);
    }
  } 
  //modify global variables
  min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

  //prepare output file names
  assert(if_name != NULL);
  w_if_name = new char[strlen(if_name)+1];
  strcpy(w_if_name, if_name);
  for(unsigned int i=0; i<strlen(w_if_name); i++) 
    if (w_if_name[i] == '/') w_if_name[i]='_';
  cerr << if_name << '\n';  
  cerr << w_if_name << '\n';  

  res4openGL_fname = new char[strlen(w_if_name) + 10];// append .4ogl.txt
  sprintf(res4openGL_fname, "%s.4ogl.txt", w_if_name);
  out_4ogl = new std::ofstream(res4openGL_fname, std::ios::out);
  assert(out_4ogl!=NULL);
  //if verbose only...
  if(verbose){
    verbose_fname  = new char[strlen(w_if_name) + 10];// append .verb.txt
    sprintf(verbose_fname, "%s.verb.txt", w_if_name);
    out_verbose = new std::ofstream( verbose_fname, std::ios::out);
    assert(out_verbose != NULL);
    CGAL::set_pretty_mode(*out_verbose);
  }

  //load the model from <mesh.off>
  PolyhedralSurf P;
  std::ifstream stream(if_name);
  stream >> P;
  fprintf(stderr, "loadMesh %d Ves %d Facets\n",
	  P.size_of_vertices(), P.size_of_facets());
  if(verbose) 
    (*out_verbose) << "Polysurf with " << P.size_of_vertices()
		   << " vertices and " << P.size_of_facets()
		   << " facets. " << std::endl;
  
  //exit if not enough points in the model
  if (min_nb_points > P.size_of_vertices())    exit(0);
  //initialize Polyhedral data : length of edges, normal of facets and
  //monge data
  P.compute_edges_length();
  P.compute_facets_normals();
  compute_differential_quantities( P);
    
  ///////////////////////////////////////////
  //RIDGES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////////////////
  //recall:
  typedef CGAL::Ridge_approximation < PolyhedralSurf,
   back_insert_iterator< std::vector<Ridge_line*> > > Ridge_approximation;

  Ridge_approximation ridge_approximation;
  std::vector<Ridge_line*> ridge_lines;


 //  //1, ne compile pas
//   std::vector<Ridge_line*>::iterator it1 ;
//   it1 = ridge_approximation.compute_all_ridges(P, 
// 					      std::back_inserter(ridge_lines),
// 					      Ridge_approximation::Tag_3);  

  //2, plante
//   std::vector<Ridge_line*>::iterator iterb = ridge_lines.begin(), it2;
//    it2 = ridge_approximation.compute_all_ridges(P,
// 					      iterb,
// 					      Ridge_approximation::Tag_3);  
  

  back_insert_iterator<std::vector<Ridge_line*> > ii(ridge_lines);

  ridge_approximation.compute_all_ridges(P,
					 ii,
					 Ridge_approximation::Tag_3);  
 

     

  std::vector<Ridge_line*>::iterator iter_lines = ridge_lines.begin(), 
    iter_end = ridge_lines.end();

  // std::cout << std::endl << iter_end - iter_lines << std::endl;
  //OpenGL output
  
  for (;iter_lines!=iter_end;iter_lines++)
    {
      (*iter_lines)->dump_4ogl(*out_4ogl);
      //  (*iter_lines)->dump_4ogl(std::cout);
    }
    
  //  //verbose txt output 
  //       if (verbose){     
  // 	(*out_verbose) 	 << std::endl ;
  //       } //END IF 
  
  //cleanup filenames
  delete res4openGL_fname; 
  out_4ogl->close(); 
  delete out_4ogl;
  if(verbose) {
    delete verbose_fname;
    out_verbose->close(); 
    delete out_verbose;
  }
  return 1;
}
 
