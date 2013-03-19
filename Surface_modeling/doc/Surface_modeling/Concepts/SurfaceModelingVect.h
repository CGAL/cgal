 /// \ingroup PkgSurfaceModeling
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for vector to be used in CGAL::Deform_mesh::rotate
 ///
 /// \code
 /// // a simple model to Vect
 /// class VectModel
 /// {
 /// public:
 ///   VectModel(double x, double y, double z) : x(x), y(y), z(z) { }
 ///   double operator[](int i) { return coors[i]; }
 ///   union
 ///   {
 ///     struct{ double coors[3]; };
 ///     struct{ double x; double y; double z; };      
 ///   }; 
 /// };
 /// \endcode
class SurfaceModelingVect
{
public:
  /// Constructor accepting three double
  SurfaceModelingVect(double x, double y, double z);
  /// return value at position i
  double operator[](int i);
};


