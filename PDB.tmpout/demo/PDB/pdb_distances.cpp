#include <CGAL/PDB/PDB.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <CGAL/PDB/distance.h>

#ifdef PDB_USE_MAGICK
#include <Magick++.h>
#endif

#include <boost/program_options.hpp>


int main(int argc, char *argv[]){
 std::vector<std::string> inputs;
  bool print_help=false;
  bool crms=false, drms=false;
  bool warn=false;
  bool all_atoms=false;
  std::string image_name;
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("crms,c", boost::program_options::bool_switch(&crms),
       "Output the cRMS between the two pdbs (after alignment).")
      ("drms,d", boost::program_options::bool_switch(&drms),
       "Output the dRMS between the two pdbs.")
      ("all-atoms,a", boost::program_options::bool_switch(&all_atoms),
       "Output the distances between all atoms, not just the C_alphas.")
      ("verbose,v", boost::program_options::bool_switch(&warn),
       "Warn about errors parsing pdb file.")
      ("image-file,i", boost::program_options::value<std::string>(&image_name),
       "Output the max distance difference between pairwise distances.")
      ("help", boost::program_options::bool_switch(&print_help), "Produce help message");

    po.add_options()("input-pdbs",
		     boost::program_options::value< std::vector<std::string> >(&inputs)->composing(),
		     "The input files. Names with a % in their name are assumed to have printf style converters for ints and will be expanded to the first n integers that correspond to names of actual files.");
    ao.add(o).add(po);

    boost::program_options::positional_options_description p;
    p.add("input-pdbs", -1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);


    if (inputs.empty() || print_help) {
      std::cout << "This program computes the distances between a collection of pdb files.\n";
      std::cout << "usage: " << argv[0] << " file1.pdb file%03d.pdb ...\n\n";
      std::cout << o << "\n";
      return EXIT_SUCCESS;
    }
  }


  std::vector<std::string> names;
  std::vector<CGAL_PDB_NS::Protein> pdbs;
  for (unsigned int i=0; i < inputs.size(); ++i){
    if (inputs[i].find('%') == std::string::npos){
      std::ifstream in(inputs[i].c_str());
      if (!in) {
	std::cerr << "Error opening file " << inputs[i] << std::endl;
      } else {
	pdbs.push_back(CGAL_PDB_NS::Protein(in));
	names.push_back(inputs[i]);
      }
    } else {
      for (unsigned int j=0; ; ++j){
	char buf[5000];
	sprintf(buf, inputs[i].c_str(), j);
	std::ifstream in(buf);
	if (!in) {
	  break;
	} else {
	  pdbs.push_back(CGAL_PDB_NS::Protein(in));
	  names.push_back(buf);
	}
      }
    }
  }

  CGAL_TNT_NS::Array2D<double> dists(pdbs.size(), pdbs.size(), 0.0);

  double max=0;
  for (unsigned int i=0; i< pdbs.size(); ++i){
    for (unsigned int j=0; j<i; ++j){
      double d;
      if (drms) {
	if (all_atoms) {
	  d= CGAL_PDB_NS::dRMS(pdbs[i], pdbs[j]);
	} else {
	  d= CGAL_PDB_NS::ca_dRMS(pdbs[i], pdbs[j]);
	}
      } else {
	if (all_atoms) {
	  d= CGAL_PDB_NS::cRMS(pdbs[i], pdbs[j]);
	}  else {
	  d= CGAL_PDB_NS::ca_cRMS(pdbs[i], pdbs[j]);
	}
      }
      dists[i][j]=dists[j][i]=d;
      max= (std::max)(d, max);
    }
  }


  std::cout << "The maximum error is " << max << std::endl;

#ifdef PDB_USE_MAGICK
  if (!image_name.empty()){
    Magick::Geometry geom(dists.dim1(),dists.dim2());

    Magick::Image im(geom, "red");
    for (int j=0; j< dists.dim1(); ++j){
      for (int k=0; k< dists.dim2(); ++k){
	double v= (std::min) (1.0, dists[j][k]/6.0);
	im.pixelColor(j, k, Magick::ColorGray(1-v));

      }
    }

    im.write(image_name.c_str());
  }

#endif

  return EXIT_SUCCESS;
}
