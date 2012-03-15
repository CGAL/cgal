#include "mesher_tester.h"
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Image_3.h>


/* DOMAIN */
typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_image_mesh_domain_3<Image,K> Image_domain;

/* TRIANGULATION */
typedef CGAL::Mesh_triangulation_3<Image_domain>::type Tr_image;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_image> C3t3_image;

/* MESHING CRITERIA */
typedef CGAL::Mesh_criteria_3<Tr_image>            Mesh_criteria_image;

/** Builds the image domain from the image name */
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
  
  std::vector<Tr_image::Point>::iterator points_begin()
  { CGAL_error_msg("No input point in 3D image (wrong --off_vertices option?)"); } 
  
  std::vector<Tr_image::Point>::iterator points_end()
  { CGAL_error_msg("No input point in 3D image (wrong --off_vertices option?)"); } 

  C3t3_image::Index points_index() { return C3t3_image::Index(); }

private:
  Domain* domain_;
  Image image_;
};



int main(int argc, char** argv)
{
  int nb_threads;
  std::string outdir;
  
	// options
	po::options_description generic("Options");
  generic.add_options()
  ("help", "Produce help message")
  ("threads", po::value<int>(&nb_threads)->default_value(2), "Run <arg> threads")
	("outdir", po::value<std::string>(&outdir)->default_value("tester_output"), "Output directory. <arg> is location");
	
	po::options_description cmdline_options("Usage", 1);
	cmdline_options.add(generic);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
	po::notify(vm);
	
	if(vm.count("help"))
	{
		std::cout << cmdline_options << std::endl;
		std::cout << "* Images:    .inr.gz files should be in data/3D_images\n";
		std::cout << "  Note: for each file toto.domain, add a toto.txt file with meshing parameters\n\n";
		return 1;
	}
	
	// iterate on data files
	std::string data_img("data/3D_images/"); // or a user defined path...
  mesh<C3t3_image,Mesh_criteria_image,Image_domain>(data_img,outdir,nb_threads); 

	return 0;
}
