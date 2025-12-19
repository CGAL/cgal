namespace CGAL {
namespace cpp11 {

  /*!
  \ingroup PkgSTLExtensionRef
  Alias to the implementation of the `std::result_of` mechanism. When all compilers
  supported by \cgal have a Standard compliant implementation of the `std::invoke_result`
  mechanism, it will become an alias to the <code>std::invoke_result</code>.

  \sa <a href=https://en.cppreference.com/w/cpp/types/result_of><code>std::result_of</code></a>
  */
  template <typename F>
  struct result_of {

    /*!
      It is a type `std::result_of<F>::%type`.
    */
    typedef unspecified_type type;
  };

} } //end of namespace CGAL::cpp11
