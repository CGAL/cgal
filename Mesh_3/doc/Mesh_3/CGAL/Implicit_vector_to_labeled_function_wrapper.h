namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

\deprecated
Use `Implicit_multi_domain_to_labeled_function_wrapper` instead.

The class `Implicit_vector_to_labeled_function_wrapper` wraps a set of implicit function [f1,f2,...] to one function F which
takes its values into N.

Let p be a point, f1 and f2 two functions:
F(p) = 0b000000(f2(p)<0)(f1(p)<0)
(Functions must take a point as input parameter and return an arithmetic value.)

It can handle at most 8 functions.

This wrapper class can be passed to `Labeled_mesh_domain_3` as first template parameter.

\tparam Function is the type of the input function.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\sa `Implicit_mesh_domain_3`.

*/
template<class Function, class BGT>
class Implicit_vector_to_labeled_function_wrapper
{
public:
  /// \name Types
  /// @{
  //!
  typedef std::vector<Function*>   Function_vector;
  //!
  typedef typename BGT::Point_3     Point_3;
  /// @}

  /// \name Creation
  /// @{
  /*!
   * \brief Construction from a vector of implicit functions.
   * \param implicit_functions the vector of implicit functions.
   */
  Implicit_vector_to_labeled_function_wrapper(const std::vector<Function*>& implicit_functions);

/// @}

}; /* end Labeled_mesh_domain_3 */
} /* end namespace CGAL */
