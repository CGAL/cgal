/*!
  \ingroup PkgJet_fitting_3Concepts
  \cgalconcept

  The concept `SvdTraits` describes the set of requirements to be 
  fulfilled by any class used to instantiate the third template 
  parameter of the class 
  `CGAL::Monge_via_jet_fitting<DataKernel,LocalKernel,SvdTraits>`. 

  It describes the linear algebra types and algorithms needed by the 
  class `CGAL::Monge_via_jet_fitting`. 

  ### Requirements ###

  The scalar type, `SvdTraits::FT`, must be the same as that of 
  the `LocalKernel` concept : `LocalKernel::FT`. 

  \hasModel `CGAL::Eigen_svd`

  \sa `LocalKernel` 

*/
class SvdTraits {
public:

  /// \name Types 
  /// @{

  /*! 
    The scalar type. 
  */ 
  typedef Hidden_type FT; 

  /*! 
    The vector type. 
  */ 
  struct Vector {
    /*! 
      initialize all the elements of the vector to zero. 
    */ 
    Vector(size_t n); 
    /*! 

     */ 
    size_t size(); 

    /*!
      return the \f$ i^{th}\f$ entry, \f$ i\f$ from \f$ 0\f$ to \f$ size()-1\f$. 
    */ 
    FT operator()(size_t i); 

    /*! 
      set the \f$ i^{th}\f$ entry to \f$ value\f$. 
    */ 
    void set(size_t i, const FT value); 

    /*! 
      return the vector as an array. 
    */ 
    FT* vector(); 
  };


  /*! 
    The matrix type. 
  */ 
  struct Matrix {
    /*! 
      initialize all the entries of the matrix to zero. 
    */ 
    Matrix(size_t n1, size_t n2); 

    /*! 

     */ 
    size_t number_of_rows(); 

    /*! 

     */ 
    size_t number_of_columns(); 

    /*! 
      return the entry at row \f$ i\f$ and column \f$ j\f$, \f$ i\f$ from \f$ 0\f$ to `number_of_rows - 1`, 
      \f$ j\f$ from \f$ 0\f$ to `number_of_columns - 1`. 
    */ 
    FT operator()(size_t i, size_t j); 

    /*! 
      set the entry at row \f$ i\f$ and column \f$ j\f$ to \f$ value\f$. 
    */ 
    void set(size_t i, size_t j, const FT value); 

    
  };
  
  /// @} 

  /// \name Operations 
  /// The concept `SvdTraits` has a linear solver using a
  /// singular value decomposition algorithm.
  /// @{

  /*! 
    Solves the system \f$ MX=B\f$ (in the least square sense if \f$ M\f$ is not 
    square) using a singular value decomposition and returns the condition 
    number of \f$ M\f$. The solution is stored in \f$ B\f$. 
  */ 
  FT solve(Matrix& M, Vector& B); 

  /// @}

}; /* end SvdTraits */

