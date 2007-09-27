// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-sb.mpg.de 
//
// ============================================================================

#ifndef CGAL_EXPONENT_VECTOR_H
#define CGAL_EXPONENT_VECTOR_H

#include <deque>
#include <iterator>
#include <algorithm>
#include <vector>
#include <CGAL/assertions.h>



CGAL_BEGIN_NAMESPACE

class Exponent_vector : public std::vector<int> {
    typedef std::vector<int> Base;
public:
    Exponent_vector(): Base(){};
    Exponent_vector(Base::size_type i): Base(i){};
    Exponent_vector(Base::size_type i, Base::value_type x): Base(i,x){};
    Exponent_vector(const Exponent_vector& x): Base ( x ){};

    template <class InputIterator>
    Exponent_vector(InputIterator, InputIterator){
        typedef typename InputIterator::value_type value_type;
        BOOST_STATIC_ASSERT(( ::boost::is_same<value_type, int>::value));
    }
    
    void push_front( int exponent ) {
        this->insert(this->begin(), exponent );
    }
    
    bool operator<( const Exponent_vector& ev ) const {
        CGAL_precondition(this->size()==ev.size());
        for (unsigned int i = this->size()-1; i >= 0; i--){
            if((*this)[i] < ev[i]) return true;
            if((*this)[i] > ev[i]) return false;
        }
        return false; 
    }
    
    void output_benchmark( std::ostream& os ) const {
        os << "( ";
        for( unsigned i = 0; i < size(); ++i ) {
            if( i != 0 )
                os << ", ";
            os << at(i); 
        }
        os << " )";
    }
};


/*
class Exponent_vector {
  public:  
    int size() { return exponents.size(); }
    
    void push_front( int exponent ) {
      exponents.push_front( exponent );
    }
    
    int& operator[]( int i ) {
      CGAL_assertion( i >= 0 && i < (int)exponents.size() );
      return exponents[i];
    }
    
    void swap( int i, int j ) {
      CGAL_assertion( i >= 0 && i < (int)exponents.size() );
      CGAL_assertion( j >= 0 && j < (int)exponents.size() );
      int tmp = exponents[i];
      exponents[i] = exponents[j];
      exponents[j] = tmp;
    }
    
    void output_exponents( std::ostream& out ) const {
      std::copy( exponents.begin(), exponents.end(), 
                 std::ostream_iterator< int > ( out, " " ) );
    }
    
    bool operator<( const Exponent_vector& ev ) const {
      return exponents < ev.exponents;
    }
     
  private:
    std::deque<int> exponents;
};

// TODO: Debug-Output
std::ostream& operator<<(std::ostream& out, const Exponent_vector& ev) {
  ev.output_exponents( out );    
  return out;
}
*/




CGAL_END_NAMESPACE

#endif // CGAL_EXPONENT_VECTOR_H
