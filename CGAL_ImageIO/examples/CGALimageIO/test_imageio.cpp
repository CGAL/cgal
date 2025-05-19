#include <CGAL/ImageIO.h>
#include <iostream>

#define SHOW(attribut) "\n  "#attribut": " << image->attribut
#define SHOWENUM(enumitem) #enumitem"=" << enumitem

int main(int argc, char** argv)
{
  if(argc > 1)
  {
    _image* image = ::_readImage(argv[1]);
    if(image)
    {
      std::cerr
        << "Image infos:"
        << "\ndimensions"
        << SHOW(xdim)
        << SHOW(ydim)
        << SHOW(zdim)
        << SHOW(vdim)
        << "\nvoxel size"
        << SHOW(vx)
        << SHOW(vy)
        << SHOW(vz)
        <<"\nimage offset"
        << SHOW(tx)
        << SHOW(ty)
        << SHOW(tz)
        <<"\nrotation vector"
        << SHOW(rx)
        << SHOW(ry)
        << SHOW(rz)
        <<"\nimage center"
        << SHOW(cx)
        << SHOW(cy)
        << SHOW(cz)
        << "\nword size (in bytes)"
        << SHOW(wdim)
        << "\nimage format"
        << "\n  " << image->imageFormat->realName
        << " (extension list: " << image->imageFormat->fileExtension << ")"
        << "\nvectors interlaced or not ("
        << SHOWENUM(VM_INTERLACED) << ", "
        << SHOWENUM(VM_NON_INTERLACED) << ", "
        << SHOWENUM(VM_SCALAR) << ")"
        << SHOW(vectMode)
        << "\nword kind ("
        << SHOWENUM(WK_FIXED) << ", "
        << SHOWENUM(WK_FLOAT) << ", "
        << SHOWENUM(WK_UNKNOWN)
        << ")"
        << SHOW(wordKind)
        << "\nword sign ("
        << SHOWENUM(SGN_SIGNED) << ", "
        << SHOWENUM(SGN_UNSIGNED) << ", "
        << SHOWENUM(SGN_UNKNOWN)
        << ")"
        << SHOW(sign)
        << "\n";

      ::_freeImage(image);
    }
    else
      std::cerr << "\"" << argv[1] << "\" is not a supported file.\n";
  }
  else
    ::printSupportedFileFormat();
}
