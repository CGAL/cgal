namespace CGAL {
namespace cpp11 {

  /*!
  \ingroup PkgStlExtension
  Alias to the tr1 implementation from boost of the `result_of` mechanism.
  When all compilers supported by %CGAL will have a Standard compliant implemention of the
  the \cpp11 `decltype` feature, it will become an alias to <code>std::result_of</code>.
  
  \sa <a href=http://www.boost.org/libs/utility/utility.htm#result_of><code>boost::result_of</code></a>
  */
  template <typename F>
  struct result_of{
    /// starting from boost version 1.44, it is `boost::tr1_result_of<F>::%type`, and 
    /// `boost::result_of<F>::%type` otherwise.
    typedef unspecified_type type;
  };
  
} } //end of namespace CGAL::cpp11
