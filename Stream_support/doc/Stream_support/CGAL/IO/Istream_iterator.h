
namespace CGAL {

/*!
\ingroup PkgIOstreams

The class `Istream_iterator` is an input iterator adaptor for the 
input stream class `Stream` and value type `T`. It is particularly 
useful for classes that are similar but not compatible to `std::istream`. 

\cgalModels `InputIterator`


*/
  template< typename T, typename Stream >
class Istream_iterator {
public:

/// \name Creation 
/// @{

/*!
creates an end-of-stream iterator. 
This is a past-the-end iterator, and it is useful 
when constructing a range. 
*/ 
Istream_iterator(); 

/*!
creates an input iterator reading from `s`. 
When `s` reaches end of stream, 
this iterator will compare equal to an end-of-stream iterator 
created using the default constructor. 
*/ 
Istream_iterator( Stream& s); 

/// @}

}; /* end Istream_iterator */
} /* end namespace CGAL */
