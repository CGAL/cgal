namespace Kinetic {
/*!
  \ingroup PkgKdsFrameworkConcepts
  \cgalConcept

  The concept `Kinetic::FunctionKernel` encapsulates all the methods for representing 
  and handing functions. The set is kept deliberately small to easy use 
  of new `Kinetic::FunctionKernel`s, but together these operations are sufficient to 
  allow the correct processing of events, handling of degeneracies, 
  usage of static data structures, run-time error checking as well as 
  run-time verification of the correctness of kinetic data structures. 
  The computation of a polynomial with the variable negated is used for 
  reversing time in kinetic data structures and can be omitted if that 
  capability is not needed. 

  \cgalHasModel `POLYNOMIAL::Kernel<RootStack>`
  \cgalHasModel `POLYNOMIAL::Filtered_kernel<RootStack>`

  \sa `Kinetic::RootEnumerator`

  \cgalHeading{Example}

  We provide several models of the concept, which are not documented 
  separately. The models of `Kinetic::SimulationTraits` all choose 
  appropriate models. However, if 
  more control is desired, we here provide examples of how to create the 
  various supported `Kinetic::FunctionKernel`. 

  A Sturm sequence based kernel which supports exact comparisons of roots of polynomials (certificate failure times): 

  \code{.cpp} 

  typedef CGAL::POLYNOMIAL::Polynomial<CGAL::Gmpq> Function; 
  typedef CGAL::POLYNOMIAL::Sturm_root_stack_traits<Function> Root_stack_traits; 
  typedef CGAL::POLYNOMIAL::Sturm_root_stack<Root_stack_traits> Root_stack; 
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel; 

  \endcode 

  A wrapper for `CORE::Expr` which implements the necessary 
  operations: 

  \code{.cpp} 

  typedef CGAL::POLYNOMIAL::CORE_kernel Function_kernel; 

  \endcode 

  A function kernel which computes approximations to the roots of the polynomials: 

  \code{.cpp} 

  typedef CGAL::POLYNOMIAL::Polynomial<double> Function; 
  typedef CGAL::POLYNOMIAL::Root_stack_default_traits<Function> Root_stack_traits; 
  typedef CGAL::POLYNOMIAL::Numeric_root_stack<Root_stack_traits> Root_stack; 
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel; 

  \endcode 

  When using the function kernel in kinetic data structures, especially 
  one that is in exact, it is useful to wrap the root stack. The wrapper 
  checks the sign of the certificate function being solved and uses that 
  to handle degeneracies. This is done by, for the inexact solvers 

  \code{.cpp} 

  typedef Kinetic::Derivitive_filter_function_kernel<Function_kernel> KDS_function_kernel; 

  \endcode 

  and for exact solvers 

  \code{.cpp} 

  typedef Kinetic::Handle_degeneracy_function_kernel<Function_kernel> KDS_function_kernel; 

  \endcode 

  For exact computations, the primary representation for roots is the 
  now standard choice of a polynomial with an associated isolating 
  interval (and interval containing exactly one distinct root of a 
  polynomial) along with whether the root has odd or even multiplicity 
  and, if needed, the Sturm sequence of the polynomial. Two intervals 
  can be compared by first seeing if the isolating intervals are 
  disjoint. If they are, then we know the ordering of the respective 
  roots. If not we can subdivide each of the intervals (using the 
  endpoints of the other interval) and repeat. In order to avoid 
  subdividing endlessly when comparing equal roots, once we subdivide a 
  constant number of times, we use the Sturm sequence of \f$ p\f$ and \f$ p'q\f$ 
  (where \f$ p\f$ and \f$ q\f$ are the two polynomials and \f$ p'\f$ is the derivative 
  of \f$ p\f$) to evaluate the sign of the second at the root of the first 
  one directly (note that this Sturm sequence is applied to a common 
  isolating interval of the roots of interest of both polynomials). 

*/

class FunctionKernel {
public:

  /// \name Types 
  /// @{

  /*!
    \ingroup PkgKdsFrameworkOtherConcepts
    \cgalConcept

    The concept `Function` represents a function. 

    \sa `FunctionKernel`
    \sa `FunctionKernel::ConstructFunction` 

    \cgalHeading{Example}

    Several ways to create functions: 

    Using `Kinetic::ConstructFunction`: 

    \code{.cpp} 

    Traits::Function_kernel::Construct_function cf= traits.function_kernel_object().construct_function_object(); 
    Traits::Kinetic_kernel::Motion_function x= cf(0.0,1.0,2.0); 
    Traits::Kinetic_kernel::Motion_function y= cf(0.0,1.0,2.0); 
    Traits::Kinetic_kernel::Point_2 pt(x,y); 

    \endcode 

    Using the constructor: 

    \code{.cpp} 

    double coefs[]={1.0, 2.0, 3.0}; 
    Traits::Kinetic_kernel::Motion_function z(coefs, coefs+3); 

    \endcode 

    Using ring operations: 

    \code{.cpp} 

    Traits::Kinetic_kernel::Motion_function z= x*z+y; 

    \endcode 

  */
  class Function {
  public:

    /// \name Types 
    /// @{

    /*!
      The number type used in describing the function. 
    */ 
    typedef unspecified_type NT; 

    /*!
      Construct a constant function from a number. 
    */ 
    Function(NT); 

    /// @} 

    /// \name Operations 
    /// @{

    /*!
      Evaluate the function at an `NT`. 
    */ 
    NT operator()(NT); 

    /// @}

  }; /* end Function */

  /*!
    The basic representational number type. 
  */ 
  typedef unspecified_type NT; 

  /*!
    A type representing the roots of a `Function`. 
  */ 
  typedef unspecified_type Root; 

  /*!
    A model of `RootStack`. These objects can be created by calling the `root_stack_object` method with a `Function` and two (optional) `Root` objects. The enumerator then enumerates all roots of the function in the open inverval defined by the two root arguments. They optional arguments default to positive and negative infinity. 
  */ 
  typedef unspecified_type Root_stack; 

  /*!
    The traits for the `Root_enumerator` class. 
  */ 
  typedef unspecified_type Root_enumerator_traits; 

  /// @}

  /// \name 
  /// Each of the following types has a corresponding `type_object`
  /// method (not explicitly documented) which takes a `Function` as an
  /// argument.
  /// @{

  /*!
    A functor which returns the sign of a 
    `Function` at a `NT` or `Root`. 
  */ 
  typedef unspecified_type Sign_at; 

  /*!
    A functor which returns sign of a function 
    immediately after a root. 
  */ 
  typedef unspecified_type Sign_after; 

  /// @}
  
  /// \name
  /// The following type is used to construct functions from a list of
  /// coefficients. To get an instance use the
  /// `construct_function_object()` method.
  /// @{  

  /*!
    \ingroup PkgKdsFrameworkOtherConcepts
    \cgalConcept

    The concept `ConstructFunction` is used to construct functions. 

    \sa `FunctionKernel` 

    \cgalHeading{Example}

    \code{.cpp} 

    Function_kernel fk; 
    Function_kernel::Construct_function cf= fk.construct_function_object(); 
    Function_kernel::Function f= cf(0,1,2,3,4,5); 

    \endcode 

  */

  class ConstructFunction {
  public:

    /// \name Operations 
    /// @{

    /*!
      This family of methods 
      takes a list of coefficients and returns a function. There can be 
      any number of coeffients passed as arguments (up to about 25 in the 
      current implementations). 
    */ 
    Function operator()(NT a, ...); 

    /// @}

  }; /* end ConstructFunction */


  /// @}
  
  /// \name
  /// The following functor likewise have a `type_object` method, but
  /// these take arguments other than a `Function`. The arguments are
  /// given below.
  /// @{  

  /*!
    This functor, creation of which 
    requires two `Root`s, returns the sign of a passed function 
    between the pair of roots. 
  */ 
  typedef unspecified_type Sign_between_roots; 

  /*!
    This functor computes the derivitive of a `Function`. Construction takes no arguments. 
  */ 
  typedef unspecified_type Differentiate; 

  /// @}
  
  /// \name
  /// The following methods do not require any arguments to get the
  /// functor and take one `Function` as a functor argument.
  /// @{  

  /*!
    Map \f$ f(x)\f$ to \f$ f(-x)\f$. 
  */ 
  typedef unspecified_type Negate_variable; 

  /// @}

}; /* end Kinetic::FunctionKernel */

} /* end namespace KineticConcepts */
