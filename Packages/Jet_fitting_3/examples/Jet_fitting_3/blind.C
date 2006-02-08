#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

#include "../../include/CGAL/Monge_via_jet_fitting.h" 
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
typedef T_PolyhedralSurf_rings<PolyhedralSurf> Poly_rings;
//Kernel for local computations
typedef double                LFT;
typedef CGAL::Cartesian<LFT>  Local_Kernel;
typedef Monge_via_jet_fitting<Data_Kernel, Local_Kernel, GSL> My_Monge_via_jet_fitting;
typedef Monge_rep<Data_Kernel> My_Monge_rep;
typedef Monge_info<Local_Kernel> My_Monge_info;
         
//Syntax requirred by Options
static const char *const optv[] = {
  "?|?",
  "f:fName <string>",		//<name>???
  "d:deg <int>",	//degree of the jet
  "m:mdegree <int>",	//degree of the Monge rep
  "a:nrings <int>",	//# rings
  "p:npoints <int>",	//# points
  "v|",//verbose?
  NULL
};
 

int main(int argc, char *argv[])
{
  // fct parameters
  unsigned int d_fitting = 2;
  unsigned int d_monge = 2;
  unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
  unsigned int nb_points_to_use = 0;//
  bool verbose = false;
  
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

  //prepare output file names
  assert(if_name != NULL);
  w_if_name = new char[strlen(if_name)];
  strcpy(w_if_name, if_name);
  for(unsigned int i=0; i<strlen(w_if_name); i++) 
    if (w_if_name[i] == '/') w_if_name[i]='_';
  cerr << if_name << '\n';  
  cerr << w_if_name << '\n';  

  res4openGL_fname = new char[strlen(w_if_name) + 9];// append .4ogl.txt
  sprintf(res4openGL_fname, "%s.4ogl.txt", w_if_name);
  out_4ogl = new std::ofstream(res4openGL_fname, std::ios::out);
  assert(out_4ogl!=NULL);
  //if verbose only...
  if(verbose){
    verbose_fname  = new char[strlen(w_if_name) + 9];// append .verb.txt
    sprintf(verbose_fname, "%s.verb.txt", w_if_name);
    out_verbose = new std::ofstream( verbose_fname, std::ios::out);
    assert(out_verbose != NULL);
    CGAL::set_pretty_mode(*out_verbose);
  }
  unsigned int nb_vertices_considered = 0;//count vertices for verbose 
  //  output

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
  unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;
  if (min_nb_points > P.size_of_vertices())    exit(0);
  //initialize Polyhedral data : length of edges, normal of facets
  P.compute_edges_length();
  P.compute_facets_normals();
  
  //container for approximation points
  std::vector<DPoint> in_points;
  //container to collect vertices of the PolyhedralSurf
  std::vector<Vertex*> current_ring, next_ring, gathered; 
  std::vector<Vertex*> *p_current_ring, *p_next_ring;

  //MAIN LOOP
  Vertex_iterator vitb = P.vertices_begin(), vite = P.vertices_end();
  for (; vitb != vite; vitb++) {
    //initialize
    Vertex* v = &(*vitb);
    unsigned int nbp = 0, //current nb of collected points
      ith = 0;	//i-th ring index
    current_ring.clear();
    next_ring.clear();
    p_current_ring = &current_ring;
    p_next_ring = &next_ring;
    gathered.clear();
    in_points.clear();  
    My_Monge_rep monge_rep;
    My_Monge_info monge_info;
      
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
	nbp += T_PolyhedralSurf_rings<PolyhedralSurf>::
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

    //skip if the nb of gathered points is to small 
    if ( gathered.size() < min_nb_points ) continue;

    //store the gathered points
    std::vector<Vertex*>::iterator itb = gathered.begin(),
      ite = gathered.end();
    CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
 
    // run the main fct : perform the fitting
    My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				   d_fitting, d_monge, 
				   monge_rep, monge_info);
 
    //switch min-max ppal curv/dir wrt the mesh orientation
    DVector normal_mesh = P.computeFacetsAverageUnitNormal(v);
    DVector normal_monge = monge_rep.n();
    if ( normal_mesh*normal_monge < 0 )
      {
	monge_rep.n() = -monge_rep.n();
	DVector temp = monge_rep.d1();
	monge_rep.d1() = monge_rep.d2();
	monge_rep.d2() = temp;
	DFT temp_curv = monge_rep.coefficients()[0];
	monge_rep.coefficients()[0] = -monge_rep.coefficients()[1];
	monge_rep.coefficients()[1] = -temp_curv;
      }
      
    //OpenGL output
    //scaling for ppal dir,
    // may be optimized with a global mean edges length computed only once
    // on all edges of P
    DFT scale_ppal_dir = P.compute_mean_edges_length_around_vertex(v)/2;
    (*out_4ogl) << v->point()  << " "
		<< monge_rep.origin_pt()  << " "
		<< monge_rep.d1() * scale_ppal_dir << " "
		<< monge_rep.d2() * scale_ppal_dir << " "
		<< monge_rep.coefficients()[0] << " "
		<< monge_rep.coefficients()[1] << " "
		<< std::endl;

    //verbose txt output 
    if (verbose){     
      (*out_verbose) << "--- vertex " <<  ++nb_vertices_considered 
		     <<	" : " << v->point() << std::endl
		     << "number of points used : " << in_points.size() << std::endl
		     << "number of rings used : " << ith << std::endl
		     << "origin : " << monge_rep.origin_pt() << std::endl
		     << "d1 : " << monge_rep.d1() << std::endl 
		     << "d2 : " << monge_rep.d2() << std::endl
		     << "n : " << monge_rep.n() << std::endl
		     << "cond_nb : " << monge_info.cond_nb() << std::endl 
		     << "pca_eigen_vals " << monge_info.pca_eigen_vals()[0] 
		     << " " << monge_info.pca_eigen_vals()[2] 
		     << " " << monge_info.pca_eigen_vals()[3]  << std::endl 
		     << "pca_eigen_vecs : " << std::endl 
		     << monge_info.pca_eigen_vecs()[0] << std::endl 
		     << monge_info.pca_eigen_vecs()[1] << std::endl 
		     << monge_info.pca_eigen_vecs()[2] << std::endl;
      if ( d_monge >= 2) 
	(*out_verbose) << "k1 : " << monge_rep.coefficients()[0] << std::endl 
		       << "k2 : " << monge_rep.coefficients()[1] << std::endl
		       << std::endl ;
    } //END IF 
  } //END FOR LOOP
  //  std::cout << "ok 1 " << std::endl;
  //cleanup filenames
  delete res4openGL_fname; 
  out_4ogl->close(); 
  delete out_4ogl;
  if(verbose) {
    delete verbose_fname;
    out_verbose->close(); 
    delete out_verbose;
  }
  //  std::cout << "ok 2 " << std::endl ;
  return 1;
}
 
