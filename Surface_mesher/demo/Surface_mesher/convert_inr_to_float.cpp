#include "../../../CGALimageIO/src/CGALimageIO/imageio/ImageIO.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{
  if(argc > 2)
  {
    std::cerr << "Openning \"" << argv[1] << "\"... ";
    _image* image = ::_readImage(argv[1]);
    std::cerr << "done.\n";
    if(image)
    {
      ::convertImageTypeToFloat(image);
      std::cerr << "Writing \"" << argv[2] << "\"... ";
      int result = ::_writeImage(image, argv[2]);
      std::cerr << "done. Result=" << result << ".\n";
      ::_freeImage(image);
      if( result != ImageIO_NO_ERROR )
	std::cerr << "Error while writing the image!\n";
      return result;
    }
    else
      std::cerr << "Cannot read the image!\n";
    return EXIT_FAILURE;
  }
  else
    std::cerr << "Usage:\n  " << argv[0] << " <input filename> <output filename>\n";
  return EXIT_FAILURE;
}
