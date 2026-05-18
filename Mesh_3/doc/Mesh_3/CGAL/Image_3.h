namespace CGAL {

/*!
\ingroup PkgMesh3Domains

The class `Image_3` holds a shared pointer to a 3D image buffer.
The implementation is based on a fork of the <a
href="https://gitlab.lip6.fr/inrimage/inrimage">InrImage library</a>.

*/
class Image_3 {
public:

  /// The default-constructor. The object is invalid until a call to `read()`.
  Image_3();

  /// Open a 3D image file.
  ///
  /// Returns `true` if the file was successfully loaded.
  bool read(const char* file);

};
} /* end namespace CGAL */
