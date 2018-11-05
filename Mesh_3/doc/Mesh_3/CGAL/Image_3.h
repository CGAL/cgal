namespace CGAL {

/*!
\ingroup PkgMesh3Domains

The class `Image_3` is a C++ wrapper around the <a
href="http://inrimage.gforge.inria.fr/">InrImage library</a>. It holds a
shared pointer to a 3D image buffer.

*/
class Image_3 {
public:

  /// The default-constructor. The object is invalid until a call to `read()`.
  Image_3();

  /// Open an 3D image file.
  ///
  /// Returns `true` if the file was sucessfully loaded.
  bool read(const char* file);

};
} /* end namespace CGAL */
