#include <mesher_tester.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

/* DOMAIN */
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits; // exact constructions here
typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> Polyhedral_domain; 

/* TRIANGULATION */
typedef CGAL::Mesh_triangulation_3<Polyhedral_domain>::type Tr_polyhedron;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_polyhedron> C3t3_polyhedron;

/* MESHING CRITERIA */
typedef CGAL::Mesh_criteria_3<Tr_polyhedron>       Mesh_criteria_polyhedron;


/** Builds the polyhedron domain from the off file name */
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
		std::cout << "* Polyhedron: .off files should be in data/Polyhedra\n";
		std::cout << "  Note: for each file toto.off, add a toto.txt file with meshing parameters\n\n";
		return 1;
	}
	
	// iterate on data files
	std::string data_img("data/Polyhedra/"); // or a user defined path...
  mesh<C3t3_polyhedron,Mesh_criteria_polyhedron,Polyhedral_domain>(data_img,outdir,nb_threads); 
  
	return 0;
}
