/*
 *  mesher_tester.cpp
 *  Mesh_3_applications
 *
 *  Created by Jane Tournois on 02/12/09.
 *  Copyright 2009 INRIA. All rights reserved.
 */

//#include <debug.h>
#include <CGAL/AABB_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
// Implicit domain
#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>

#include <CGAL/make_mesh_3.h>
#include "../examples/Mesh_3/implicit_functions.h"
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Mesh_3/Mesh_global_optimizer.h>

/* INPUTS */
// Polyhedral domain
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Implicit domain is above
// 3D Image
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

/* OUTPUT */
#include <CGAL/IO/File_medit.h>

/* tools */
#include <sstream>
#include <stdlib.h>
#include <algorithm>

/* OPTIONS and PARAMETERS */
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace CGAL::parameters; //to avoid verbose function and named parameters call

/* FILE SYSTEM */
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem; 

/* DOMAIN */
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
// Polyhedral domain
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits; // exact constructions here
typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> Polyhedral_domain; 
// Implicit domain
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> I_Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<I_Function> I_Function_wrapper;
typedef I_Function_wrapper::Function_vector I_Function_vector;
typedef CGAL::Labeled_mesh_domain_3<I_Function_wrapper, K> Implicit_domain;
// 3D Image
typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_image_mesh_domain_3<Image,K> Image_domain;

/* TRIANGULATION */
//Polyhedral domain
typedef CGAL::Mesh_triangulation_3<Polyhedral_domain>::type Tr_polyhedron;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_polyhedron> C3t3_polyhedron;
//Implicit domain
typedef CGAL::Mesh_triangulation_3<Implicit_domain>::type Tr_implicit;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_implicit> C3t3_implicit;
//3D Image
typedef CGAL::Mesh_triangulation_3<Image_domain>::type Tr_image;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_image> C3t3_image;

/* MESHING CRITERIA */
//Polyhedral domain
typedef CGAL::Mesh_criteria_3<Tr_polyhedron>       Mesh_criteria_polyhedron;
//typedef Mesh_criteria_polyhedron::Facet_criteria   Facet_criteria_polyhedron;
//typedef Mesh_criteria_polyhedron::Cell_criteria    Cell_criteria_polyhedron;
//Implicit domain
typedef CGAL::Mesh_criteria_3<Tr_implicit>         Mesh_criteria_implicit;
//typedef Mesh_criteria_implicit::Facet_criteria     Facet_criteria_implicit;
//typedef Mesh_criteria_implicit::Cell_criteria      Cell_criteria_implicit;
//3D Image
typedef CGAL::Mesh_criteria_3<Tr_image>            Mesh_criteria_image;
//typedef Mesh_criteria_image::Facet_criteria        Facet_criteria_image;
//typedef Mesh_criteria_image::Cell_criteria         Cell_criteria_image;


template <typename T>
T set_arg(const std::string& param_name,
		  const std::string& param_string,
		  const po::variables_map& vm)
{
	if(vm.count(param_name))
	{
		T param_value = vm[param_name].as<T>();
		//std::cout << param_string << ": " << param_value << "\n";
		return param_value;
	}
	else
	{
		//std::cout << param_string << " ignored.\n";
		return T();
	}
}

void set_implicit_function(I_Function_vector& v,
						   I_Function& f,
						   const std::string& function_name,
						   const po::variables_map& vm)
{
	if(vm.count(function_name))
	{
		v.push_back(f);
		std::cout << function_name << " ";
	}
}

std::vector<std::string> split_line(const std::string& str)
{
	std::vector<std::string> args;
	
	std::string::size_type lastPos = str.find_first_not_of(" ", 0);
	std::string::size_type pos = str.find_first_of(" ", lastPos);
	while(pos != std::string::npos || lastPos != std::string::npos)
	{
		args.push_back(str.substr(lastPos, pos-lastPos));
		lastPos = str.find_first_not_of(" ", pos);
		pos = str.find_first_of(" ", lastPos);
	}
	return args;
}

template<typename C3t3>
void save_histogram(std::string& filename, 
					const C3t3& c3t3)
{
	std::vector<int> histo(181,0);

	for (typename C3t3::Cell_iterator cit = c3t3.cells_begin() ;
		 cit != c3t3.cells_end() ;
		 ++cit)
	{
		if( !c3t3.is_in_complex(cit))
			continue;
		
		typedef typename K::Point_3 Point_3;
		const Point_3& p0 = cit->vertex(0)->point();
		const Point_3& p1 = cit->vertex(1)->point();
		const Point_3& p2 = cit->vertex(2)->point();
		const Point_3& p3 = cit->vertex(3)->point();
		
		double a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0,p1,p2,p3)));
		histo[std::floor(a)] += 1;
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p2, p1, p3)));
		histo[std::floor(a)] += 1;
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p3, p1, p2)));
		histo[std::floor(a)] += 1;
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p2, p0, p3)));
		histo[std::floor(a)] += 1;
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p3, p0, p2)));
		histo[std::floor(a)] += 1;
		a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p2, p3, p0, p1)));
		histo[std::floor(a)] += 1;
	}
	std::ofstream file(filename.c_str());
	std::copy(histo.begin(), histo.end(), std::ostream_iterator<int>(file, "\n"));	
}


template<typename MeshingCriteria>
MeshingCriteria get_parameters(const std::string& param_line,
							   po::variables_map& vm)
{
	po::options_description mesh("Mesh generation parameters");
	mesh.add_options()
  ("mesh", "Generate mesh")
	("facet_angle", po::value<double>(), "Set facet angle bound")
	("facet_size", po::value<double>(), "Set facet size bound")
	("facet_error", po::value<double>(), "Set approximation error bound")
	("tet_shape", po::value<double>(), "Set tet radius-edge bound")
	("tet_size", po::value<double>(), "Set tet size bound");
	
	po::options_description implicit_functions("Implicit functions");
	implicit_functions.add_options()
	("torus", "Mesh torus function")
	("sphere", "Mesh sphere function")
	("chair", "Mesh chair function")
	("tanglecube", "Mesh tanglecube function")
	("cube", "Mesh cube function")
	("ellipsoid", "Mesh ellipsoid function")
	("heart", "Mesh heart function")
	("octic", "Mesh octic function");
	
	po::options_description optim("Optimization parameters");
	optim.add_options()
	("exude", "Exude mesh after refinement")
	("perturb", po::value<double>(), "Perturb mesh after refinement (sliver removal)")
	("lloyd", po::value<int>(), "Lloyd smoothing after refinement. arg is max nb iterations")
	("odt", po::value<int>(), "ODT smoothing after refinement. arg is max nb iterations");
	
	po::options_description additional_options("Options");
	additional_options.add_options()
	("off_vertices", po::value<int>(), "Use polyhedron vertices as initialization step")
	("no_label_rebind", "Don't rebind cell labels in medit output")
	("show_patches", "Show surface patches in medit output");
	
	po::options_description cmdline_options("Usage", 1);
	cmdline_options.add(mesh).add(implicit_functions).add(optim).add(additional_options);
	
	std::vector<std::string> args = split_line(param_line);
	po::store(po::command_line_parser(args).options(cmdline_options).run(), vm);
	po::notify(vm);

	double facet_angle = set_arg<double>("facet_angle", "Facet angle", vm);
	double facet_size  = set_arg<double>("facet_size", "Facet size", vm);
	double facet_error = set_arg<double>("facet_error", "Facet approximation error", vm);
	double tet_shape = set_arg<double>("tet_shape","Tet shape (radius-edge)", vm);
	double tet_size = set_arg<double>("tet_size","Tet size", vm);

	typename MeshingCriteria::Facet_criteria fc(facet_angle, facet_size, facet_error);
	typename MeshingCriteria::Cell_criteria cc(tet_shape, tet_size);
	return MeshingCriteria(fc, cc);	
}

void get_implicit_function(const po::variables_map& vm,
						   I_Function_vector& fv)
{
	// Define functions
	I_Function f1( &torus_function);
	I_Function f2( &sphere_function<3>);
	I_Function f3( &chair_function);
	I_Function f4( &tanglecube_function);
	I_Function f5( &cube_function);
	I_Function f6( &ellipsoid_function);
	I_Function f7( &heart_function);
	I_Function f8( &octic_function);
	
	std::cout << "Function(s): ";
	set_implicit_function(fv, f1, "torus", vm);
	set_implicit_function(fv, f2, "sphere", vm);
	set_implicit_function(fv, f3, "chair", vm);
	set_implicit_function(fv, f4, "tanglecube", vm);
	set_implicit_function(fv, f5, "cube", vm);
	set_implicit_function(fv, f6, "ellipsoid", vm);
	set_implicit_function(fv, f7, "heart", vm);
	set_implicit_function(fv, f8, "octic", vm);
	std::cout << "\n\n";
	
	if(fv.empty())
		std::cout << "Warning: No implicit function set.\n";	
}


template <class Domain>
class Domain_builder
{
public:
  /*void build(const std::string& str);
  Domain domain();*/
};

template <>
class Domain_builder<Polyhedral_domain>
{
  typedef Polyhedral_domain Domain;
public:
  Domain_builder(const std::string& str)
  : domain_(NULL)
  {
    std::ifstream input(str.c_str());
    input >> polyhedron_;
    domain_ = new Domain(polyhedron_);
  }
  
  ~Domain_builder() { delete domain_; }
  
  Domain& domain() { return *domain_; }

private:
  Domain* domain_;
  Polyhedron polyhedron_;
};

template <>
class Domain_builder<Image_domain>
{
  typedef Image_domain Domain;
public:
  Domain_builder(const std::string& str)
  : domain_(NULL)
  {
    image_.read(str.c_str());
    domain_ = new Domain(image_, 1e-6);
  }
  
  ~Domain_builder() { delete domain_; }
  
  Domain& domain() { return *domain_; }
  
private:
  Domain* domain_;
  Image image_;
};

//template <>
//class Domain_builder<Implicit_domain>
//{
//  typedef Implicit_domain Domain;
//public:
//  Domain_builder()
//  :	f1_( &torus_function);
//	, f2_( &sphere_function<3>);
//	, f3_( &chair_function);
//	, f4_( &tanglecube_function);
//	, f5_( &cube_function);
//  , f6_( &ellipsoid_function);
//	, f7_( &heart_function);
//	, f8_( &octic_function);
//  
//  void build(const std::string& str)
//  {
//    Image image;
//    image.read(str.c_str());
//    delete domain_;
//    domain_ = new Domain(image, 1e-6);
//  }
//  
//  Domain& domain() { return *domain_; }
//  
//private:
//  Domain* domain_;
//  
//  // Define functions
//	I_Function f1_;
//	I_Function f2_;
//	I_Function f3_;
//	I_Function f4_;
//	I_Function f5_;
//	I_Function f6_;
//	I_Function f7_;
//	I_Function f8_;
//};



template <class C3T3, class Mesh_criteria, class Domain>
void mesh(const std::string& data, const po::variables_map& vm)
{
  if(!fs::is_directory(data))
    std::cout << "!! Problem while reading " << data << "\n";
	
  std::string output_dir;
  if(vm.count("outdir"))
    output_dir = vm["outdir"].as<std::string>();
  else output_dir = data + "output";
  if(!fs::is_directory(output_dir) && !fs::create_directory(output_dir))
    std::cout << "!! Problem while creating " << output_dir << "\n";
  
  fs::path path(data);
  for(fs::directory_iterator it(path); it != fs::directory_iterator(); ++it)
  {
    if(fs::is_directory(*it)
       || (fs::extension(*it) != ".off" && (fs::extension(*it) != ".inr" && ( fs::extension(fs::basename(*it)) != ".inr" || fs::extension(*it) != ".gz"))) )
      continue;
    
    std::string line_param;
    std::string filename(fs::basename(*it));
    if ( fs::extension(*it) == ".gz" )
      filename = fs::basename(fs::basename(*it));
    std::string filename_param(data + filename + ".txt");
    
    std::ifstream file_param(filename_param.data()); //parameters
    if(!file_param) 
    {
      std::cout << "Could not read parameters in : '" << filename_param << "'. Next file." << std::endl;
      continue;
    }
    unsigned int i = 1;
    
    // Timer
    CGAL::Timer timer;
    timer.start();
    
    // we keep c3t3 between lines
    C3T3 c3t3_save;
    
    //Load the domain
    std::cout << "****** [" << filename << "] Create domain...";
    std::flush(std::cout);
    Domain_builder<Domain> domain_builder(it->path().string());
    std::cout <<"done (" << timer.time() << "s) ******\n\n";
    
    while(std::getline(file_param,line_param))
    {
      std::cout << "*** Meshing " << filename << "[" << i << "] with : " << line_param << std::endl;

      po::variables_map vm_p;
      Mesh_criteria mcp = get_parameters<Mesh_criteria>(line_param, vm_p);
      
      //Mesh generation (reload domain and rebuild c3t3)
      if ( vm_p.count("mesh") )
      {
        timer.stop();
        timer.reset();
        timer.start();
        std::cout << "  Generate mesh...";
        std::flush(std::cout);
        c3t3_save = CGAL::make_mesh_3<C3T3>(domain_builder.domain(), mcp, no_exude(), no_perturb());
        std::cout << "done (" << timer.time() << "s - "
                  << c3t3_save.triangulation().number_of_vertices() << " vertices, "
                  << c3t3_save.number_of_cells() << " cells)\n";
      }
      
      C3T3 c3t3 = c3t3_save;
      
      //Optimization
      timer.stop();
      timer.reset();
      timer.start();
      if(vm_p.count("lloyd"))
      {
        std::cout << "  Lloyd optimization...";
        std::flush(std::cout);
        CGAL::lloyd_optimize_mesh_3(c3t3, domain_builder.domain(), max_iteration_number=vm_p["lloyd"].as<int>());
        std::cout << "done (" << timer.time() << "s)\n";
      }
      timer.stop();
      timer.reset();
      timer.start();
      if(vm_p.count("odt"))
      {
        std::cout << "  Odt optimization...";
        std::flush(std::cout);
        CGAL::odt_optimize_mesh_3(c3t3, domain_builder.domain(), max_iteration_number=vm_p["odt"].as<int>());
        std::cout << "done (" << timer.time() << "s)\n";
      }	
      timer.stop();
      timer.reset();
      timer.start();
      if(vm_p.count("perturb"))
      {
        std::cout << "  Perturbation...";
        std::flush(std::cout);
        CGAL::perturb_mesh_3(c3t3, domain_builder.domain(), time_limit=vm_p["perturb"].as<double>());
        std::cout << "done (" << timer.time() << "s)\n";
      }
      timer.stop();
      timer.reset();
      timer.start();
      if(vm_p.count("exude"))
      {
        std::cout << "  Exudation...";
        std::flush(std::cout);
        CGAL::exude_mesh_3(c3t3);
        std::cout << "done (" << timer.time() << "s)\n";
      }
      timer.stop();
      timer.reset();
      timer.start();
      
      //save mesh
      std::cout << "  Save mesh...";
      std::stringstream ssout;
      ssout << i;				
      std::string output_filename = output_dir +"/" + filename + "-out-" + ssout.str().c_str() + ".mesh";
      std::ofstream medit_file(output_filename.c_str());
      c3t3.output_to_medit(medit_file, !vm_p.count("no_label_rebind"), vm_p.count("show_patches"));
      
      //save histogram
      std::cout << "done. \n  Save histogram...";
      std::string histo_filename = output_dir +"/" + filename + "-histo-" + ssout.str().c_str() + ".txt";
      save_histogram<C3T3>(histo_filename, c3t3);
      i++;
      std::cout << "done.\n\n\n";
    }
  }
}


int main(int argc, char** argv)
{
	// options
	po::options_description generic("Options");
	generic.add_options()("help", "Produce help message")
	("polyhedron", "Test the polyhedral domain mesher")
	("implicit",   "Test the implicit domain mesher")
	("image",      "Test the 3D image domain mesher")
	("outdir", po::value<std::string>(), "Output directory. arg is location");
	
	po::options_description cmdline_options("Usage", 1);
	cmdline_options.add(generic);//.add(others);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
	po::notify(vm);
	
	if(vm.count("help") || argc < 2)
	{
		std::cout << cmdline_options << std::endl;
		std::cout << "* Polyhedra: .off files should be in applications/data/Polyhedra\n";
		std::cout << "* Images:    .inr.gz files should be in applications/data/3D_images\n";
		std::cout << "  Note: for each file toto.domain, add a toto.txt file with meshing parameters\n\n";
		return 1;
	}
	
	//what are we testing
	bool mesh_polyhedra = vm.count("polyhedron");
	bool mesh_implicit = vm.count("implicit");
	bool mesh_images = vm.count("image");
	
	// iterate on data files
	std::string data_dir("data"); // or a user defined path...
	fs::path data_path(data_dir); 
	
	if(mesh_polyhedra)
	{
		std::string data_poly = data_dir + "/Polyhedra/";
    mesh<C3t3_polyhedron,Mesh_criteria_polyhedron,Polyhedral_domain>(data_poly,vm); 
	}

	if(mesh_images)
	{
		std::string data_img = data_dir + "/3D_images/";
    mesh<C3t3_image,Mesh_criteria_image,Image_domain>(data_img,vm); 
  }
	
	if(mesh_implicit)
	{
		std::string data_imp = data_dir + "/Implicit/";
		//mesh<C3t3_implicit,Mesh_criteria_implicit>(data_imp,vm); 
	}
	return 0;
}
