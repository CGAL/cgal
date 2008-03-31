#include <CGAL/ImageIO.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>\

#define SHOW(attribut) "\n  "#attribut": " << image->attribut
#define SHOWENUM(enumitem) #enumitem"=" << enumitem

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  if(argc == 2)
  {
    const std::string input_filename = argv[1];
    std::stringstream sub_filename;
    sub_filename << input_filename.substr(0, input_filename.size()-4) << ".txt";

    std::ifstream input_aux_file(sub_filename.str().c_str());
    int xdim, ydim, zdim;
    input_aux_file >> xdim >> ydim >> zdim;
    if(!input_aux_file) 
    {
      cout << "\nERROR: cannot read the auxiliary file '" << sub_filename.str() << "'!\n";
      return 1;
    }

    _image* image = ::_readImage_raw(argv[1],
                                     xdim, ydim, zdim);
    if(image)
    {
      input_aux_file >> image->vx >> image->vy >> image->vz;
      if(!input_aux_file)
        image->vx = image->vy = image->vz = 1.;

      std::stringstream output_filename;
      output_filename << input_filename.substr(0, input_filename.size()-4)
#ifdef CGAL_USE_ZLIB
                      << ".inr.gz";
#else
                      << ".inr";
#endif
      return ::_writeImage(image, output_filename.str().c_str());
    }
  }
  else
    std::cerr << "Usage:\n"
              << "  " << argv[0] << " filename.raw\n";
}
