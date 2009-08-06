// Boost. Iterable Range Library (rangelib)
//
// Copyright 2003-2004 John Torjo (john@torjo.com) and Matthew Wilson (matthew@synesis.com.au)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef CGAL_PDB_BOOST_RTL_ALGO_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_ALGO_HPP_INCLUDED

#include <algorithm>
#include <iterator>
#include <iosfwd>
#include <boost/config.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>

/*
    Notes:
    1. For the purpose of algorithms, EVERY container IS a range.
    
    2. At this time (06th Nov 2004), all STL algorithms are implemented, except for:
    - partial_sort_copy


*/


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib { 


template< class value_type, class char_type, class traits_type>
irange< ::std::istream_iterator< value_type, char_type, traits_type> > istream_range( ::std::basic_istream<char_type,traits_type> & s) {
    typedef ::std::istream_iterator< value_type, char_type, traits_type> it_type;
    typedef irange< it_type > range_type;
    return range_type( it_type(s), it_type() );
}









// helper - for passing an iterator instead of a range as a second parameter,
//          to an algorithm that needs two ranges (example: mismatch)
//
// You can have 
// - rng::mismatch( rng1, rng2); // passing two ranges
// - rng::mismatch( rng, i_(iterator) ); // passing the second param as an iterator
// - rng::mismatch( rng, i_param(iterator) ); // same as above, but more descriptive
template< class it>
struct i_param_t {
    typedef it iterator;
    i_param_t( iterator first) : m_first(first) {}
    iterator begin() const { return m_first;}

    // provided for copy_backward, which if a range is the second param, will copy starting from its end
    iterator end() const { return m_first;}
private:
    iterator m_first;
};

// shorter form: 'i_'
template<class it> inline i_param_t<it> i_( it val) { return i_param_t<it>(val); }
// longer and more descriptable/readable form 'i_param'
template<class it> inline i_param_t<it> i_param( it val) { return i_param_t<it>(val); }










// for algorithms returning a range
enum range_return_type {
    iter,
    from_beg,
    to_end
};



namespace rng {







// easy access to .find() from associative containers
//
// Example:
// if ( rng::coll_find(m_cached_days, current_day) ) {
// }
//
// equivalent to:
// if ( m_cached_days.find(current_day) != m_cached_days.end() ) {
// }

template< class container, class value_type> inline crange<const container> 
coll_find( const container & c, const value_type & val) {
    return crange<const container>( c.find(val), c.end() );
}

template< class container, class value_type> inline crange<container>
coll_find( container & c, const value_type & val) {
    //typedef typename container::iterator iterator;
    //typedef irange<iterator> range_type;
    return crange<container>( c.find(val), c.end() );
}




template<class container, class r> inline void
erase( container & c, const crange<r> & rng) {
    c.erase( rng.begin(), rng.end() );
}

namespace detail {
    template<int i> struct int_keeper {};

    template<class container, class r> inline crange<r>
    erase_current_impl(container & c, const crange<r> & rng, int_keeper<0>) {
        // non-associative container
        return crange<r>( c.erase( rng.begin() ), rng.end() );
    }

    template<class container, class r> inline crange<r>
    erase_current_impl(container & c, const crange<r> & rng, int_keeper<1>) {
        // associative container
        typename container::iterator next = rng.begin();
        ++next;
        c.erase( rng.begin() );
        return crange<r>( next, rng.end() );
    }
}

template<class container, class r> inline crange<r>
erase_current(container & c, const crange<r> & rng) {
    return detail::erase_current_impl(c, rng, 
        detail::int_keeper< ::CGAL::PDB::internal::rangelib::detail::has_key_type<container>::value>() );
}








//template<class in_it> typename iterator_traits<in_it>::difference_type
//distance(in_it first, in_it last);
template<class r> inline typename range_finder<const r>::difference_type distance(const r& val1) { 
    return ::std::distance( val1.begin(), val1.end() ); 
} 

template<class r> inline typename range_finder<r>::difference_type distance(r& val1) { 
    return ::std::distance( val1.begin(), val1.end() ); 
}




// advance an iterator by a certain distance
template< class iterator, class Distance> inline iterator advanced( iterator i, Distance n) {
    ::std::advance(i, n);
    return i;
}




//template<class in_it, class Function>
//Function for_each(in_it first, in_it last, Function f);
template<class r, class function> inline function for_each(const r& val1, function val2) { 
    return ::std::for_each( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class function> inline function for_each(r& val1, function val2) { 
    return ::std::for_each( val1.begin(), val1.end(), val2 ); 
}



// template<class in_it, class T> typename iterator_traits<in_it>::difference_type
// count(in_it first, in_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::difference_type count(const r& val1, const T& val2) { 
    return ::std::count( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class T> inline typename range_finder<r>::difference_type count(r& val1, const T& val2) { 
    return ::std::count( val1.begin(), val1.end(), val2 ); 
}




// template<class in_it, class Predicate> typename iterator_traits<in_it>::difference_type
// count_if(in_it first, in_it last, Predicate pred);
template<class r, class T> inline typename range_finder<const r>::difference_type count_if(const r& val1, T val2) { 
    return ::std::count_if( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class T> inline typename range_finder<r>::difference_type count_if(r& val1, T val2) { 
    return ::std::count_if( val1.begin(), val1.end(), val2 ); 
}





// template<class in_it1, class in_it2> bool 
// equal(in_it1 first1, in_it1 last1, in_it2 first2);
template<class r1, class r2> inline bool equal(r1& val1, const r2& val2) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline bool equal(r1& val1, r2& val2) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline bool equal(const r1& val1, const r2& val2) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline bool equal(const r1& val1, r2& val2) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin() ); 
}




// template<class in_it1, class in_it2, class BinaryPredicate> bool 
// equal(in_it1 first1, in_it1 last1, in_it2 first2, BinaryPredicate pred);
template<class r1, class r2, class pred> inline bool equal(r1& val1, const r2& val2, pred val3) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool equal(r1& val1, r2& val2, pred val3) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool equal(const r1& val1, const r2& val2, pred val3) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool equal(const r1& val1, r2& val2, pred val3) { 
    return ::std::equal( val1.begin(), val1.end(), val2.begin(), val3 ); 
}




// template<class in_it, class out_it> out_it 
// copy(in_it first, in_it last, out_it result);
template<class r, class out_it> inline out_it copy(const r& val1, out_it val2) { 
    return ::std::copy( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class out_it> inline out_it copy(r& val1, out_it val2) { 
    return ::std::copy( val1.begin(), val1.end(), val2 ); 
}








// template<class bid_it1, class bid_it2> bid_it2 
// copy_backward(bid_it1 first, bid_it1 last, bid_it2 result);
//     FIXME - here, we could actually return a range
template<class r1, class r2> inline typename range_finder<const r2>::iterator_type copy_backward(r1& val1, const r2& val2) { 
    return ::std::copy_backward( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<r2>::iterator_type copy_backward(r1& val1, r2& val2) { 
    return ::std::copy_backward( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<const r2>::iterator_type copy_backward(const r1& val1, const r2& val2) { 
    return ::std::copy_backward( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<r2>::iterator_type copy_backward(const r1& val1, r2& val2) { 
    return ::std::copy_backward( val1.begin(), val1.end(), val2.begin() ); 
}





// template<class fwd_it1, class fwd_it2> fwd_it2 
// swap_ranges(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2);
//     FIXME - here, we could actually return a range
template<class r1, class r2> inline typename range_finder<const r2>::iterator_type swap_ranges(r1& val1, const r2& val2) { 
    return ::std::swap_ranges( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<r2>::iterator_type swap_ranges(r1& val1, r2& val2) { 
    return ::std::swap_ranges( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<const r2>::iterator_type swap_ranges(const r1& val1, const r2& val2) { 
    return ::std::swap_ranges( val1.begin(), val1.end(), val2.begin() ); 
} 
template<class r1, class r2> inline typename range_finder<r2>::iterator_type swap_ranges(const r1& val1, r2& val2) { 
    return ::std::swap_ranges( val1.begin(), val1.end(), val2.begin() ); 
}






// template<class fwd_it1, class fwd_it2>
// void iter_swap(fwd_it1 a, fwd_it2 b);
template<class r1, class r2> inline void iter_swap(r1& val1, const r2& val2) { 
    return ::std::iter_swap( val1.begin(), val2.begin() ); 
} 
template<class r1, class r2> inline void iter_swap(r1& val1, r2& val2) { 
    return ::std::iter_swap( val1.begin(), val2.begin() ); 
} 
template<class r1, class r2> inline void iter_swap(const r1& val1, const r2& val2) { 
    return ::std::iter_swap( val1.begin(), val2.begin() ); 
} 
template<class r1, class r2> inline void iter_swap(const r1& val1, r2& val2) { 
    return ::std::iter_swap( val1.begin(), val2.begin() ); 
}







// template<class in_it, class out_it, class UnaryOperation> out_it 
// transform(in_it first, in_it last, out_it result, UnaryOperation op);
template<class r, class out_it, class op> inline out_it transform(const r& val1, out_it val2, op val3) { 
    return ::std::transform( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class op> inline out_it transform(r& val1, out_it val2, op val3) { 
    return ::std::transform( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class BinaryOperation> out_it 
// transform(in_it1 first1, in_it1 last1, in_it2 first2, out_it result, BinaryOperation binary_op);
template<class r1, class r2, class out_it, class op> inline out_it transform(r1& val1, const r2& val2, out_it val3, op val4) { 
    return ::std::transform( val1.begin(), val1.end(), val2.begin(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class op> inline out_it transform(r1& val1, r2& val2, out_it val3, op val4) { 
    return ::std::transform( val1.begin(), val1.end(), val2.begin(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class op> inline out_it transform(const r1& val1, const r2& val2, out_it val3, op val4) { 
    return ::std::transform( val1.begin(), val1.end(), val2.begin(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class op> inline out_it transform(const r1& val1, r2& val2, out_it val3, op val4) { 
    return ::std::transform( val1.begin(), val1.end(), val2.begin(), val3, val4 ); 
}




// template<class fwd_it, class T> void 
// replace(fwd_it first, fwd_it last, const T& old_value, const T& new_value);
template<class r, class T, class T_same> inline void replace(const r& val1, const T& val2, const T_same& val3) { 
    ::std::replace( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class T, class T_same> inline void replace(r& val1, const T& val2, const T_same& val3) { 
    ::std::replace( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class fwd_it, class Predicate, class T> void 
// replace_if(fwd_it first, fwd_it last, Predicate pred, const T& new_value);
template<class r, class pred, class T> inline void replace_if(const r& val1, pred val2, const T& val3) { 
    ::std::replace_if( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class pred, class T> inline void replace_if(r& val1, pred val2, const T& val3) { 
    ::std::replace_if( val1.begin(), val1.end(), val2, val3 ); 
}





// template<class in_it, class out_it, class T> out_it 
// replace_copy(in_it first, in_it last, out_it result, const T& old_value, const T& new_value);
template<class r, class out_it, class T, class T_same> inline out_it replace_copy(const r& val1, out_it val2, const T& val3, const T_same& val4) { 
    return ::std::replace_copy( val1.begin(), val1.end(), val2, val3, val4 ); 
} 
template<class r, class out_it, class T, class T_same> inline out_it replace_copy(r& val1, out_it val2, const T& val3, const T_same& val4) { 
    return ::std::replace_copy( val1.begin(), val1.end(), val2, val3, val4 ); 
}




// template<class Iterator, class out_it, class Predicate, class T> out_it 
// replace_copy_if(Iterator first, Iterator last, out_it result, Predicate pred, const T& new_value);
template<class r, class out_it, class pred, class T> inline out_it replace_copy_if(const r& val1, out_it val2, pred val3, const T& val4) { 
    return ::std::replace_copy_if( val1.begin(), val1.end(), val2, val3, val4 ); 
} 
template<class r, class out_it, class pred, class T> inline out_it replace_copy_if(r& val1, out_it val2, pred val3, const T& val4) { 
    return ::std::replace_copy_if( val1.begin(), val1.end(), val2, val3, val4 ); 
}




// template<class fwd_it, class T> void 
// fill(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline void fill(const r& val1, const T& val2) { 
    (void) ::std::fill( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class T> inline void fill(r& val1, const T& val2) { 
    (void) ::std::fill( val1.begin(), val1.end(), val2 ); 
}





// template<class out_it, class Size, class T> void 
// fill_n(out_it first, Size n, const T& value);
//useless - 





// template<class fwd_it, class Generator> void 
// generate(fwd_it first, fwd_it last, Generator gen);
template<class r, class generator> inline void generate(const r& val1, generator val2) { 
    (void) ::std::generate( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class generator> inline void generate(r& val1, generator val2) { 
    (void) ::std::generate( val1.begin(), val1.end(), val2 ); 
}





// 
// template<class out_it, class Size, class Generator> void 
// generate_n(out_it first, Size n, Generator gen);
// useless - 





// template<class in_it, class out_it> out_it 
// unique_copy(in_it first, in_it last, out_it result);
template<class r, class out_it> inline out_it unique_copy(const r& val1, out_it val2) { 
    return ::std::unique_copy( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class out_it> inline out_it unique_copy(r& val1, out_it val2) { 
    return ::std::unique_copy( val1.begin(), val1.end(), val2 ); 
}




// template<class in_it, class out_it, class BinaryPredicate> out_it 
// unique_copy(in_it first, in_it last, out_it result, BinaryPredicate pred);
template<class r, class out_it, class pred> inline out_it unique_copy(const r& val1, out_it val2, pred val3) { 
    return ::std::unique_copy( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class pred> inline out_it unique_copy(r& val1, out_it val2, pred val3) { 
    return ::std::unique_copy( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class bid_it> void 
// reverse(bid_it first, bid_it last);
template<class r> inline void reverse(const r& val1) { 
    (void) ::std::reverse( val1.begin(), val1.end() ); 
} 
template<class r> inline void reverse(r& val1) { 
    (void) ::std::reverse( val1.begin(), val1.end() ); 
}




// template<class bid_it, class out_it> out_it 
// reverse_copy(bid_it first, bid_it last, out_it result);
template<class r, class out_it> inline out_it reverse_copy(const r& val1, out_it val2) { 
    return ::std::reverse_copy( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class out_it> inline out_it reverse_copy(r& val1, out_it val2) { 
    return ::std::reverse_copy( val1.begin(), val1.end(), val2 ); 
}




// template<class fwd_it> void 
// rotate(fwd_it first, fwd_it middle, fwd_it last);
template<class r, class fwd_it> inline void rotate(const r& val1, fwd_it val2) { 
    (void) ::std::rotate( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class fwd_it> inline void rotate(r& val1, fwd_it val2) { 
    (void) ::std::rotate( val1.begin(), val1.end(), val2 ); 
}





// template<class rand_it> void 
// random_shuffle(rand_it first, rand_it last);
template<class r> inline void random_shuffle(const r& val1) { 
    (void) ::std::random_shuffle( val1.begin(), val1.end() ); 
} 
template<class r> inline void random_shuffle(r& val1) { 
    (void) ::std::random_shuffle( val1.begin(), val1.end() ); 
}




// template<class rand_it, class RandomNumberGenerator> void 
// random_shuffle(rand_it first, rand_it last, RandomNumberGenerator& rand);
template<class r, class generator> inline void random_shuffle(const r& val1, generator& val2) { 
    (void) ::std::random_shuffle( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class generator> inline void random_shuffle(r& val1, generator& val2) { 
    (void) ::std::random_shuffle( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// sort(rand_it first, rand_it last);
template<class r> inline void sort(const r& val1) { 
    (void) ::std::sort( val1.begin(), val1.end() ); 
} 
template<class r> inline void sort(r& val1) { 
    (void) ::std::sort( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare>
// void sort(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void sort(const r& val1, pred val2) { 
    (void) ::std::sort( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void sort(r& val1, pred val2) { 
    (void) ::std::sort( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// stable_sort(rand_it first, rand_it last);
template<class r> inline void stable_sort(const r& val1) { 
    (void) ::std::stable_sort( val1.begin(), val1.end() ); 
} 
template<class r> inline void stable_sort(r& val1) { 
    (void) ::std::stable_sort( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare> void 
// stable_sort(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void stable_sort(const r& val1, pred val2) { 
    (void) ::std::stable_sort( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void stable_sort(r& val1, pred val2) { 
    (void) ::std::stable_sort( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// partial_sort(rand_it first, rand_it middle, rand_it last);
template<class r, class rand_it> inline void partial_sort(const r& val1, rand_it val2) { 
    (void) ::std::partial_sort( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class rand_it> inline void partial_sort(r& val1, rand_it val2) { 
    (void) ::std::partial_sort( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it, class Compare> void 
// partial_sort(rand_it first, rand_it middle, rand_it last, Compare comp);
template<class r, class rand_it, class pred> inline void partial_sort(const r& val1, rand_it val2, pred val3) { 
    (void) ::std::partial_sort( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class rand_it, class pred> inline void partial_sort(r& val1, rand_it val2, pred val3) { 
    (void) ::std::partial_sort( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class rand_it> void 
// nth_element(rand_it first, rand_it nth, rand_it last);
template<class r, class rand_it> inline void nth_element(const r& val1, rand_it val2) { 
    (void) ::std::nth_element( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class rand_it> inline void nth_element(r& val1, rand_it val2) { 
    (void) ::std::nth_element( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it, class Compare> void 
// nth_element(rand_it first, rand_it nth, rand_it last, Compare comp);
template<class r, class rand_it, class pred> inline void nth_element(const r& val1, rand_it val2, pred val3) { 
    (void) ::std::nth_element( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class rand_it, class pred> inline void nth_element(r& val1, rand_it val2, pred val3) { 
    (void) ::std::nth_element( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class fwd_it, class T> bool 
// binary_search(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline bool binary_search(const r& val1, const T& val2) { 
    return ::std::binary_search( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class T> inline bool binary_search(r& val1, const T& val2) { 
    return ::std::binary_search( val1.begin(), val1.end(), val2 ); 
}




// template<class fwd_it, class T, class Compare> bool 
// binary_search(fwd_it first, fwd_it last, const T& value, Compare comp);
template<class r, class T, class pred> inline bool binary_search(const r& val1, const T& val2, pred val3) { 
    return ::std::binary_search( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class T, class pred> inline bool binary_search(r& val1, const T& val2, pred val3) { 
    return ::std::binary_search( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class in_it1, class in_it2, class out_it> out_it 
// merge(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
template<class r1, class r2, class out_it> inline out_it merge(r1& val1, const r2& val2, out_it val3) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it merge(r1& val1, r2& val2, out_it val3) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it merge(const r1& val1, const r2& val2, out_it val3) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it merge(const r1& val1, r2& val2, out_it val3) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// merge(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
template<class r1, class r2, class out_it, class pred> inline out_it merge(r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it merge(r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it merge(const r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it merge(const r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::merge( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
}




// template<class bid_it> void 
// inplace_merge(bid_it first, bid_it middle, bid_it last);
template<class r, class bid_it> inline void inplace_merge(const r& val1, bid_it val2) { 
    (void) ::std::inplace_merge( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class bid_it> inline void inplace_merge(r& val1, bid_it val2) { 
    (void) ::std::inplace_merge( val1.begin(), val1.end(), val2 ); 
}




// template<class bid_it, class Compare> void 
// inplace_merge(bid_it first, bid_it middle, bid_it last, Compare comp);
template<class r, class bid_it, class pred> inline void inplace_merge(const r& val1, bid_it val2, pred val3) { 
    (void) ::std::inplace_merge( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class bid_it, class pred> inline void inplace_merge(r& val1, bid_it val2, pred val3) { 
    (void) ::std::inplace_merge( val1.begin(), val1.end(), val2, val3 ); 
}





// template<class in_it1, class in_it2> bool 
// includes(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2);
template<class r1, class r2> inline bool includes(r1& val1, const r2& val2) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool includes(r1& val1, r2& val2) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool includes(const r1& val1, const r2& val2) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool includes(const r1& val1, r2& val2) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
}




// template<class in_it1, class in_it2, class Compare> bool 
// includes(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, Compare comp);
template<class r1, class r2, class comp> inline bool includes(r1& val1, const r2& val2, comp val3) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class comp> inline bool includes(r1& val1, r2& val2, comp val3) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class comp> inline bool includes(const r1& val1, const r2& val2, comp val3) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class comp> inline bool includes(const r1& val1, r2& val2, comp val3) { 
    return ::std::includes( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it> out_it 
// set_union(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
template<class r1, class r2, class out_it> inline out_it set_union(r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_union(r1& val1, r2& val2, out_it val3) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_union(const r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_union(const r1& val1, r2& val2, out_it val3) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class Compare>
// out_it set_union(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
template<class r1, class r2, class out_it, class pred> inline out_it set_union(r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_union(r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_union(const r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_union(const r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_union( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
}




// template<class in_it1, class in_it2, class out_it> out_it 
// set_intersection(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
template<class r1, class r2, class out_it> inline out_it set_intersection(r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_intersection(r1& val1, r2& val2, out_it val3) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_intersection(const r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_intersection(const r1& val1, r2& val2, out_it val3) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// set_intersection(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
template<class r1, class r2, class out_it, class pred> inline out_it set_intersection(r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_intersection(r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_intersection(const r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_intersection(const r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_intersection( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
}




// template<class in_it1, class in_it2, class out_it> out_it 
// set_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
template<class r1, class r2, class out_it> inline out_it set_difference(r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_difference(r1& val1, r2& val2, out_it val3) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_difference(const r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_difference(const r1& val1, r2& val2, out_it val3) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// set_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
template<class r1, class r2, class out_it, class pred> inline out_it set_difference(r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_difference(r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_difference(const r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_difference(const r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
}




// template<class in_it1, class in_it2, class out_it> out_it
// set_symmetric_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
template<class r1, class r2, class out_it> inline out_it set_symmetric_difference(r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_symmetric_difference(r1& val1, r2& val2, out_it val3) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_symmetric_difference(const r1& val1, const r2& val2, out_it val3) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class out_it> inline out_it set_symmetric_difference(const r1& val1, r2& val2, out_it val3) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class in_it1, class in_it2, class out_it, class Compare> out_it
// set_symmetric_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
template<class r1, class r2, class out_it, class pred> inline out_it set_symmetric_difference(r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_symmetric_difference(r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_symmetric_difference(const r1& val1, const r2& val2, out_it val3, pred val4) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
} 
template<class r1, class r2, class out_it, class pred> inline out_it set_symmetric_difference(const r1& val1, r2& val2, out_it val3, pred val4) { 
    return ::std::set_symmetric_difference( val1.begin(), val1.end(), val2.begin(), val2.end(), val3, val4 ); 
}





// template<class rand_it> void 
// push_heap(rand_it first, rand_it last);
template<class r> inline void push_heap(const r& val1) { 
    (void) ::std::push_heap( val1.begin(), val1.end() ); 
} 
template<class r> inline void push_heap(r& val1) { 
    (void) ::std::push_heap( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare> void 
// push_heap(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void push_heap(const r& val1, pred val2) { 
    (void) ::std::push_heap( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void push_heap(r& val1, pred val2) { 
    (void) ::std::push_heap( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// pop_heap(rand_it first, rand_it last);
template<class r> inline void pop_heap(const r& val1) { 
    (void) ::std::pop_heap( val1.begin(), val1.end() ); 
} 
template<class r> inline void pop_heap(r& val1) { 
    (void) ::std::pop_heap( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare> void 
// pop_heap(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void pop_heap(const r& val1, pred val2) { 
    (void) ::std::pop_heap( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void pop_heap(r& val1, pred val2) { 
    (void) ::std::pop_heap( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// make_heap(rand_it first, rand_it last);
template<class r> inline void make_heap(const r& val1) { 
    (void) ::std::make_heap( val1.begin(), val1.end() ); 
} 
template<class r> inline void make_heap(r& val1) { 
    (void) ::std::make_heap( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare> void 
// make_heap(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void make_heap(const r& val1, pred val2) { 
    (void) ::std::make_heap( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void make_heap(r& val1, pred val2) { 
    (void) ::std::make_heap( val1.begin(), val1.end(), val2 ); 
}




// template<class rand_it> void 
// sort_heap(rand_it first, rand_it last);
template<class r> inline void sort_heap(const r& val1) { 
    (void) ::std::sort_heap( val1.begin(), val1.end() ); 
} 
template<class r> inline void sort_heap(r& val1) { 
    (void) ::std::sort_heap( val1.begin(), val1.end() ); 
}




// template<class rand_it, class Compare> void 
// sort_heap(rand_it first, rand_it last, Compare comp);
template<class r, class pred> inline void sort_heap(const r& val1, pred val2) { 
    (void) ::std::sort_heap( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline void sort_heap(r& val1, pred val2) { 
    (void) ::std::sort_heap( val1.begin(), val1.end(), val2 ); 
}




// template<class in_it1, class in_it2> bool 
// lexicographical_compare(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2);
template<class r1, class r2> inline bool lexicographical_compare(r1& val1, const r2& val2) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool lexicographical_compare(r1& val1, r2& val2) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool lexicographical_compare(const r1& val1, const r2& val2) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
} 
template<class r1, class r2> inline bool lexicographical_compare(const r1& val1, r2& val2) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end() ); 
}




// template<class in_it1, class in_it2, class Compare> bool 
// lexicographical_compare(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, Compare comp);
template<class r1, class r2, class pred> inline bool lexicographical_compare(r1& val1, const r2& val2, pred val3) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool lexicographical_compare(r1& val1, r2& val2, pred val3) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool lexicographical_compare(const r1& val1, const r2& val2, pred val3) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
} 
template<class r1, class r2, class pred> inline bool lexicographical_compare(const r1& val1, r2& val2, pred val3) { 
    return ::std::lexicographical_compare( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ); 
}




// template<class bid_it> bool 
// next_permutation(bid_it first, bid_it last);
template<class r> inline bool next_permutation(const r& val1) { 
    return ::std::next_permutation( val1.begin(), val1.end() ); 
} 
template<class r> inline bool next_permutation(r& val1) { 
    return ::std::next_permutation( val1.begin(), val1.end() ); 
}




// template<class bid_it> bool 
// next_permutation(bid_it first, bid_it last);
template<class r, class pred> inline bool next_permutation(const r& val1, pred val2) { 
    return ::std::next_permutation( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline bool next_permutation(r& val1, pred val2) { 
    return ::std::next_permutation( val1.begin(), val1.end(), val2 ); 
}




// template<class bid_it> bool 
// prev_permutation(bid_it first, bid_it last);
template<class r> inline bool prev_permutation(const r& val1) { 
    return ::std::prev_permutation( val1.begin(), val1.end() ); 
} 
template<class r> inline bool prev_permutation(r& val1) { 
    return ::std::prev_permutation( val1.begin(), val1.end() ); 
}




// template<class bid_it, class Compare> bool 
// prev_permutation(bid_it first, bid_it last, Compare comp);
template<class r, class pred> inline bool prev_permutation(const r& val1, pred val2) { 
    return ::std::prev_permutation( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class pred> inline bool prev_permutation(r& val1, pred val2) { 
    return ::std::prev_permutation( val1.begin(), val1.end(), val2 ); 
}



// template<class in_it, class out_it, class T> out_it 
// remove_copy(in_it first, in_it last, out_it result, const T& value);
template<class r, class out_it, class T> inline out_it remove_copy(const r& val1, out_it val2, const T& val3) { 
    return ::std::remove_copy( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class T> inline out_it remove_copy(r& val1, out_it val2, const T& val3) { 
    return ::std::remove_copy( val1.begin(), val1.end(), val2, val3 ); 
}




// template<class in_it, class out_it, class Predicate> out_it 
// remove_copy_if(in_it first, in_it last, out_it result, Predicate pred);
template<class r, class out_it, class pred> inline out_it remove_copy_if(const r& val1, out_it val2, pred val3) { 
    return ::std::remove_copy_if( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class pred> inline out_it remove_copy_if(r& val1, out_it val2, pred val3) { 
    return ::std::remove_copy_if( val1.begin(), val1.end(), val2, val3 ); 
}










////////////////////////////////////////////////////////////////////////////
// from <numerics>



// template <class in_it, class T> T 
// accumulate(in_it first, in_it last, T init);
template<class r, class T> inline T accumulate(const r& val1, T val2) { 
    return ::std::accumulate( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class T> inline T accumulate(r& val1, T val2) { 
    return ::std::accumulate( val1.begin(), val1.end(), val2 ); 
}


// template <class in_it, class T, class Func> T 
// accumulate(in_it first, in_it last, T init, Func binary_op);
template<class r, class T, class Func> inline T accumulate(const r& val1, T val2, Func val3) { 
    return ::std::accumulate( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class T, class Func> inline T accumulate(r& val1, T val2, Func val3) { 
    return ::std::accumulate( val1.begin(), val1.end(), val2, val3 ); 
}

// template <class in_it1, class in_it2, class T> T 
// inner_product(in_it1 first1, in_it1 last1, in_it2 first2, T init);
template<class r1, class r2, class T> inline T inner_product(r1& val1, const r2& val2, T init) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init ); 
} 
template<class r1, class r2, class T> inline T inner_product(r1& val1, r2& val2, T init) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init ); 
} 
template<class r1, class r2, class T> inline T inner_product(const r1& val1, const r2& val2, T init) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init ); 
} 
template<class r1, class r2, class T> inline T inner_product(const r1& val1, r2& val2, T init) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init ); 
}

// template <class in_it1, class in_it2, class T, class Func1, class Func2> T 
// inner_product(in_it1 first1, in_it1 last1, in_it2 first2, T init, Func1 binary_op1, Func2 binary_op2);
template<class r1, class r2, class T, class Func1, class Func2> inline T inner_product(r1& val1, const r2& val2, T init, Func1 f1, Func2 f2) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init, f1, f2 ); 
} 
template<class r1, class r2, class T, class Func1, class Func2> inline T inner_product(r1& val1, r2& val2, T init, Func1 f1, Func2 f2) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init, f1, f2 ); 
} 
template<class r1, class r2, class T, class Func1, class Func2> inline T inner_product(const r1& val1, const r2& val2, T init, Func1 f1, Func2 f2) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init, f1, f2 ); 
} 
template<class r1, class r2, class T, class Func1, class Func2> inline T inner_product(const r1& val1, r2& val2, T init, Func1 f1, Func2 f2) { 
    return ::std::inner_product( val1.begin(), val1.end(), val2.begin(), init, f1, f2 ); 
}

// template <class in_it, class out_it> out_it 
// partial_sum(in_it first, in_it last, out_it result);
template<class r, class out_it> inline out_it partial_sum(const r& val1, out_it val2) { 
    return ::std::partial_sum( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class out_it> inline out_it partial_sum(r& val1, out_it val2) { 
    return ::std::partial_sum( val1.begin(), val1.end(), val2 ); 
}

// template <class in_it, class out_it, class Func> out_it 
// partial_sum(in_it first, in_it last, out_it result, Func binary_op);
template<class r, class out_it, class Func> inline out_it partial_sum(const r& val1, out_it val2, Func val3) { 
    return ::std::partial_sum( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class Func> inline out_it partial_sum(r& val1, out_it val2, Func val3) { 
    return ::std::partial_sum( val1.begin(), val1.end(), val2, val3 ); 
}

// template <class in_it, class out_it> out_it 
// adjacent_difference(in_it first, in_it last, out_it result);
template<class r, class out_it> inline out_it adjacent_difference(const r& val1, out_it val2) { 
    return ::std::adjacent_difference( val1.begin(), val1.end(), val2 ); 
} 
template<class r, class out_it> inline out_it adjacent_difference(r& val1, out_it val2) { 
    return ::std::adjacent_difference( val1.begin(), val1.end(), val2 ); 
}

// template <class in_it, class out_it, class Func> out_it 
// adjacent_difference(in_it first, in_it last, out_it result, Func binary_op);
template<class r, class out_it, class Func> inline out_it adjacent_difference(const r& val1, out_it val2, Func val3) { 
    return ::std::adjacent_difference( val1.begin(), val1.end(), val2, val3 ); 
} 
template<class r, class out_it, class Func> inline out_it adjacent_difference(r& val1, out_it val2, Func val3) { 
    return ::std::adjacent_difference( val1.begin(), val1.end(), val2, val3 ); 
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////


// 
// 
// 
// 
// * These functions return a range

/*
    Right now, we always return a range, from the element "found",
    up to the end of the original range.

    This might not necessary be what you want.

    Therefore, you can prefix an extra template parameter, which can be any of:
    - iter      - return only the found iterator, just like the std:: algorithm
    - from_beg  - returns a range: [begin,found)
    - to_end    - (default) returns a range: [found, end)

    Example:

    std::vector<int> v;
    find(v, 1); // default
    find<iter>(v, 1); // return the iterator at which we've found 1
    find<from_beg>(v, 1); // return range( v.begin(), rng::find(v,1));
    find<to_end>(v, 1); // default


*/



// template<class in_it, class T> in_it 
// find(in_it first, in_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::range_type find(const r& val1, const T& val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::find( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class T> inline typename range_finder<r>::range_type find(r& val1, const T& val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::find( val1.begin(), val1.end(), val2 ), val1.end() ); 
}



// template<class in_it, class Predicate> in_it 
// find_if(in_it first, in_it last, Predicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type find_if(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::find_if( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type find_if(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::find_if( val1.begin(), val1.end(), val2 ), val1.end() ); 
}




// template<class fwd_it1, class fwd_it2> fwd_it1 
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<class r1, class r2> inline typename range_finder<const r1>::range_type find_end(const r1& val1, r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type find_end(r1& val1, r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<const r1>::range_type find_end(const r1& val1, const r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type find_end(r1& val1, const r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
}




// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type find_end(const r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type find_end(r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type find_end(const r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type find_end(r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
}




// template<class fwd_it1, class fwd_it2> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<class r1, class r2> inline typename range_finder<const r1>::range_type find_first_of(const r1& val1, r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type find_first_of(r1& val1, r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<const r1>::range_type find_first_of(const r1& val1, const r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type find_first_of(r1& val1, const r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
}



// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type find_first_of(const r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type find_first_of(r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type find_first_of(const r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type find_first_of(r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
}




// template<class fwd_it> fwd_it 
// adjacent_find(fwd_it first, fwd_it last);
template<class r> inline typename range_finder<const r>::range_type adjacent_find(const r& val1) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::adjacent_find( val1.begin(), val1.end() ), val1.end() ); 
} 
template<class r> inline typename range_finder<r>::range_type adjacent_find(r& val1) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::adjacent_find( val1.begin(), val1.end() ), val1.end() ); 
}



// template<class fwd_it, class BinaryPredicate> fwd_it 
// adjacent_find(fwd_it first, fwd_it last, BinaryPredicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type adjacent_find(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::adjacent_find( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type adjacent_find(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::adjacent_find( val1.begin(), val1.end(), val2 ), val1.end() ); 
}





// template<class fwd_it1, class fwd_it2> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<class r1, class r2> inline typename range_finder<const r1>::range_type search(const r1& val1, r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type search(r1& val1, r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<const r1>::range_type search(const r1& val1, const r2& val2) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
} 
template<class r1, class r2> inline typename range_finder<r1>::range_type search(r1& val1, const r2& val2) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end() ); 
}



// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type search(const r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type search(r1& val1, r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<const r1>::range_type search(const r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<const r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
} 
template<class r1, class r2, class pred> inline typename range_finder<r1>::range_type search(r1& val1, const r2& val2, pred val3) { 
    typedef typename range_finder<r1>::range_type range_type; 
    return range_type( ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end() ); 
}




// template<class fwd_it, class Size, class T> fwd_it 
// search_n(fwd_it first, fwd_it last, Size count, const T& value);
template<class r, class Size, class T> inline typename range_finder<const r>::range_type search_n(const r& val1, Size val2, const T& val3) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::search_n( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
} 
template<class r, class Size, class T> inline typename range_finder<r>::range_type search_n(r& val1, Size val2, const T& val3) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::search_n( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
}




// template<class fwd_it, class Size, class T, class BinaryPredicate> fwd_it1 
// search_n(fwd_it first, fwd_it last, Size count, const T& value, BinaryPredicate pred);
template<class r, class Size, class T, class pred> inline typename range_finder<const r>::range_type search_n(const r& val1, Size val2, const T& val3, pred val4) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::search_n( val1.begin(), val1.end(), val2, val3, val4 ), val1.end() ); 
} 
template<class r, class Size, class T, class pred> inline typename range_finder<r>::range_type search_n(r& val1, Size val2, const T& val3, pred val4) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::search_n( val1.begin(), val1.end(), val2, val3, val4 ), val1.end() ); 
}


// 
// ***** return the range of elements to be removed;
//        use them in correlation with rng::erase.
//






// template<class fwd_it, class T> fwd_it 
// remove(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::range_type remove(const r& val1, const T& val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::remove( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class T> inline typename range_finder<r>::range_type remove(r& val1, const T& val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::remove( val1.begin(), val1.end(), val2 ), val1.end() ); 
}



// template<class fwd_it, class Predicate> fwd_it 
// remove_if(fwd_it first, fwd_it last, Predicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type remove_if(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::remove_if( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type remove_if(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::remove_if( val1.begin(), val1.end(), val2 ), val1.end() ); 
}



// template<class fwd_it> fwd_it 
// unique(fwd_it first, fwd_it last);
template<class r> inline typename range_finder<const r>::range_type unique(const r& val1) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::unique( val1.begin(), val1.end() ), val1.end() ); 
} 
template<class r> inline typename range_finder<r>::range_type unique(r& val1) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::unique( val1.begin(), val1.end() ), val1.end() ); 
}




// template<class fwd_it, class BinaryPredicate> fwd_it 
// unique(fwd_it first, fwd_it last, BinaryPredicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type unique(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::unique( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type unique(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::unique( val1.begin(), val1.end(), val2 ), val1.end() ); 
}






// // return the range that satisfies the predicate
// template<class bid_it, class Predicate> bid_it 
// partition(bid_it first, bid_it last, Predicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type partition(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::partition( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type partition(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::partition( val1.begin(), val1.end(), val2 ), val1.end() ); 
}




// template<class bid_it, class Predicate> bid_it 
// stable_partition(bid_it first, bid_it last, Predicate pred);
template<class r, class pred> inline typename range_finder<const r>::range_type stable_partition(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::stable_partition( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type stable_partition(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::stable_partition( val1.begin(), val1.end(), val2 ), val1.end() ); 
}


// 
// template<class in_it, class rand_it> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last);
// FIXME - watch out what we return!!!

// 
// template<class in_it, class rand_it, class Compare> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last, Compare comp);
// FIXME - watch out what we return!!!





// template<class fwd_it, class T> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::range_type lower_bound(const r& val1, const T& val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::lower_bound( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class T> inline typename range_finder<r>::range_type lower_bound(r& val1, const T& val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::lower_bound( val1.begin(), val1.end(), val2 ), val1.end() ); 
}




// template<class fwd_it, class T, class Compare> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
template<class r, class T, class pred> inline typename range_finder<const r>::range_type lower_bound(const r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::lower_bound( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
} 
template<class r, class T, class pred> inline typename range_finder<r>::range_type lower_bound(r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::lower_bound( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
}




// template<class fwd_it, class T> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::range_type upper_bound(const r& val1, const T& val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::upper_bound( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class T> inline typename range_finder<r>::range_type upper_bound(r& val1, const T& val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::upper_bound( val1.begin(), val1.end(), val2 ), val1.end() ); 
}




// template<class fwd_it, class T, class Compare> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
template<class r, class T, class pred> inline typename range_finder<const r>::range_type upper_bound(const r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::upper_bound( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
} 
template<class r, class T, class pred> inline typename range_finder<r>::range_type upper_bound(r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::upper_bound( val1.begin(), val1.end(), val2, val3 ), val1.end() ); 
}




// template<class fwd_it, class T> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value);
template<class r, class T> inline typename range_finder<const r>::range_type equal_range(const r& val1, const T& val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2 ) ); 
} 
template<class r, class T> inline typename range_finder<r>::range_type equal_range(r& val1, const T& val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2 ) ); 
}





// template<class fwd_it, class T, class Compare> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value, Compare comp);
template<class r, class T, class pred> inline typename range_finder<const r>::range_type equal_range(const r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<const r>::range_type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2, val3 ) ); 
} 
template<class r, class T, class pred> inline typename range_finder<r>::range_type equal_range(r& val1, const T& val2, pred val3) { 
    typedef typename range_finder<r>::range_type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2, val3 ) ); 
}





// template<class fwd_it> fwd_it 
// min_element(fwd_it first, fwd_it last);
template<class r> inline typename range_finder<const r>::range_type min_element(const r& val1) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::min_element( val1.begin(), val1.end() ), val1.end() ); 
} 
template<class r> inline typename range_finder<r>::range_type min_element(r& val1) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::min_element( val1.begin(), val1.end() ), val1.end() ); 
}




// template<class fwd_it, class Compare> fwd_it 
// min_element(fwd_it first, fwd_it last, Compare comp);
template<class r, class pred> inline typename range_finder<const r>::range_type min_element(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::min_element( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type min_element(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::min_element( val1.begin(), val1.end(), val2 ), val1.end() ); 
}




// template<class fwd_it> fwd_it 
// max_element(fwd_it first, fwd_it last);
template<class r> inline typename range_finder<const r>::range_type max_element(const r& val1) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::max_element( val1.begin(), val1.end() ), val1.end() ); 
} 
template<class r> inline typename range_finder<r>::range_type max_element(r& val1) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::max_element( val1.begin(), val1.end() ), val1.end() ); 
}




// template<class fwd_it, class Compare> fwd_it 
// max_element(fwd_it first, fwd_it last, Compare comp);
template<class r, class pred> inline typename range_finder<const r>::range_type max_element(const r& val1, pred val2) { 
    typedef typename range_finder<const r>::range_type range_type; 
    return range_type( ::std::max_element( val1.begin(), val1.end(), val2 ), val1.end() ); 
} 
template<class r, class pred> inline typename range_finder<r>::range_type max_element(r& val1, pred val2) { 
    typedef typename range_finder<r>::range_type range_type; 
    return range_type( ::std::max_element( val1.begin(), val1.end(), val2 ), val1.end() ); 
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
    IMPLEMENTATION of

    prefix an extra template parameter, which can be any of:
    - iter      - return only the found iterator, just like the std:: algorithm
    - from_beg  - returns a range: [begin,found)
    - to_end    - (default) returns a range: [found, end)

    Example:

    std::vector<int> v;
    find(v, 1); // default
    find<iter>(v, 1); // return the iterator at which we've found 1
    find<from_beg>(v, 1); // return range( v.begin(), rng::find(v,1));
    find<to_end>(v, 1); // default


*/



namespace detail {
    // find the result of the range
    template<class r, range_return_type ret> struct range_result {
        typedef typename range_finder<r>::range_type range_type;
        typedef range_type type;
    };
    template<class r> struct range_result<r, iter> {
        typedef typename range_finder<r>::range_type range_type;
        typedef typename range_type::iterator type;
    };

    template<range_return_type ret> struct range_result_tag {};

    template<class iterator> inline iterator return_range(iterator begin, iterator found, iterator end, range_result_tag<iter>) {
        return found;
    }

    template<class iterator> inline ::std::pair<iterator,iterator> return_range(iterator begin, iterator found, iterator end, range_result_tag<from_beg>) {
        return std::make_pair(begin, found);
    }

    template<class iterator> inline ::std::pair<iterator,iterator> return_range(iterator begin, iterator found, iterator end, range_result_tag<to_end>) {
        return ::std::make_pair(found, end);
    }

}




/*
// template<class in_it, class T> in_it 
// find(in_it first, in_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type find(const r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type find(r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}
*/


// template<class in_it, class T> in_it 
// find(in_it first, in_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type find(const r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type find(r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class in_it, class Predicate> in_it 
// find_if(in_it first, in_it last, Predicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type find_if(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::find_if( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type find_if(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::find_if( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it1, class fwd_it2> fwd_it1 
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type find_end(const r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  find_end(r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type find_end(const r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  find_end(r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type find_end(const r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  find_end(r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type find_end(const r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  find_end(r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_end( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it1, class fwd_it2> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type find_first_of(const r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  find_first_of(r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type find_first_of(const r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  find_first_of(r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type find_first_of(const r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  find_first_of(r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type find_first_of(const r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  find_first_of(r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::find_first_of( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it> fwd_it 
// adjacent_find(fwd_it first, fwd_it last);
template<range_return_type result_tag, class r> inline typename detail::range_result<const r,result_tag>::type adjacent_find(const r& val1) { 
    return detail::return_range( val1.begin(), ::std::adjacent_find( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r> inline typename detail::range_result<r,result_tag>::type adjacent_find(r& val1) { 
    return detail::return_range( val1.begin(), ::std::adjacent_find( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class fwd_it, class BinaryPredicate> fwd_it 
// adjacent_find(fwd_it first, fwd_it last, BinaryPredicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type adjacent_find(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::adjacent_find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type adjacent_find(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::adjacent_find( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}





// template<class fwd_it1, class fwd_it2> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type search(const r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  search(r1& val1, r2& val2) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<const r1,result_tag>::type search(const r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2> inline typename detail::range_result<r1,result_tag>::type  search(r1& val1, const r2& val2) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type search(const r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  search(r1& val1, r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<const r1,result_tag>::type search(const r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r1, class r2, class pred> inline typename detail::range_result<r1,result_tag>::type  search(r1& val1, const r2& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::search( val1.begin(), val1.end(), val2.begin(), val2.end(), val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class Size, class T> fwd_it 
// search_n(fwd_it first, fwd_it last, Size count, const T& value);
template<range_return_type result_tag, class r, class Size, class T> inline typename detail::range_result<const r,result_tag>::type search_n(const r& val1, Size val2, const T& val3) { 
    return detail::return_range( val1.begin(), ::std::search_n( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class Size, class T> inline typename detail::range_result<r,result_tag>::type search_n(r& val1, Size val2, const T& val3) { 
    return detail::return_range( val1.begin(), ::std::search_n( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class Size, class T, class BinaryPredicate> fwd_it1 
// search_n(fwd_it first, fwd_it last, Size count, const T& value, BinaryPredicate pred);
template<range_return_type result_tag, class r, class Size, class T, class pred> inline typename detail::range_result<const r,result_tag>::type search_n(const r& val1, Size val2, const T& val3, pred val4) { 
    return detail::return_range( val1.begin(), ::std::search_n( val1.begin(), val1.end(), val2, val3, val4 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class Size, class T, class pred> inline typename detail::range_result<r,result_tag>::type search_n(r& val1, Size val2, const T& val3, pred val4) { 
    return detail::return_range( val1.begin(), ::std::search_n( val1.begin(), val1.end(), val2, val3, val4 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}


// 
// ***** return the range of elements to be removed;
//        use them in correlation with rng::erase.
//






// template<class fwd_it, class T> fwd_it 
// remove(fwd_it first, fwd_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type remove(const r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::remove( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type remove(r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::remove( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class fwd_it, class Predicate> fwd_it 
// remove_if(fwd_it first, fwd_it last, Predicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type remove_if(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::remove_if( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type remove_if(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::remove_if( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



// template<class fwd_it> fwd_it 
// unique(fwd_it first, fwd_it last);
template<range_return_type result_tag, class r> inline typename detail::range_result<const r,result_tag>::type unique(const r& val1) { 
    return detail::return_range( val1.begin(), ::std::unique( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r> inline typename detail::range_result<r,result_tag>::type unique(r& val1) { 
    return detail::return_range( val1.begin(), ::std::unique( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class BinaryPredicate> fwd_it 
// unique(fwd_it first, fwd_it last, BinaryPredicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type unique(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::unique( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type unique(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::unique( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}






// // return the range that satisfies the predicate
// template<class bid_it, class Predicate> bid_it 
// partition(bid_it first, bid_it last, Predicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type partition(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::partition( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type partition(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::partition( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class bid_it, class Predicate> bid_it 
// stable_partition(bid_it first, bid_it last, Predicate pred);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type stable_partition(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::stable_partition( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type stable_partition(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::stable_partition( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}


// 
// template<class in_it, class rand_it> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last);
// FIXME - watch out what we return!!!

// 
// template<class in_it, class rand_it, class Compare> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last, Compare comp);
// FIXME - watch out what we return!!!





// template<class fwd_it, class T> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type lower_bound(const r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::lower_bound( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type lower_bound(r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::lower_bound( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class T, class Compare> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<const r,result_tag>::type lower_bound(const r& val1, const T& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::lower_bound( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<r,result_tag>::type lower_bound(r& val1, const T& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::lower_bound( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class T> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type upper_bound(const r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::upper_bound( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type upper_bound(r& val1, const T& val2) { 
    return detail::return_range( val1.begin(), ::std::upper_bound( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class T, class Compare> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<const r,result_tag>::type upper_bound(const r& val1, const T& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::upper_bound( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<r,result_tag>::type upper_bound(r& val1, const T& val2, pred val3) { 
    return detail::return_range( val1.begin(), ::std::upper_bound( val1.begin(), val1.end(), val2, val3 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}



/*
nonsense

// template<class fwd_it, class T> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value);
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<const r,result_tag>::type equal_range(const r& val1, const T& val2) { 
    typedef typename detail::range_result<const r,result_tag>::type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2 ) ); 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T> inline typename detail::range_result<r,result_tag>::type equal_range(r& val1, const T& val2) { 
    typedef typename detail::range_result<r,result_tag>::type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2 ) ); 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class T, class Compare> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value, Compare comp);
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<const r,result_tag>::type equal_range(const r& val1, const T& val2, pred val3) { 
    typedef typename detail::range_result<const r,result_tag>::type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2, val3 ) ); 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class T, class pred> inline typename detail::range_result<r,result_tag>::type equal_range(r& val1, const T& val2, pred val3) { 
    typedef typename detail::range_result<r,result_tag>::type range_type; 
    typedef typename range_type::iterator iterator; 
    return ::CGAL::PDB::internal::rangelib::irange<iterator>( ::std::equal_range( val1.begin(), val1.end(), val2, val3 ) ); 
        detail::range_result_tag<result_tag>() ); 
}
*/




// template<class fwd_it> fwd_it 
// min_element(fwd_it first, fwd_it last);
template<range_return_type result_tag, class r> inline typename detail::range_result<const r,result_tag>::type min_element(const r& val1) { 
    return detail::return_range( val1.begin(), ::std::min_element( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r> inline typename detail::range_result<r,result_tag>::type min_element(r& val1) { 
    return detail::return_range( val1.begin(), ::std::min_element( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class Compare> fwd_it 
// min_element(fwd_it first, fwd_it last, Compare comp);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type min_element(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::min_element( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type min_element(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::min_element( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it> fwd_it 
// max_element(fwd_it first, fwd_it last);
template<range_return_type result_tag, class r> inline typename detail::range_result<const r,result_tag>::type max_element(const r& val1) { 
    return detail::return_range( val1.begin(), ::std::max_element( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r> inline typename detail::range_result<r,result_tag>::type max_element(r& val1) { 
    return detail::return_range( val1.begin(), ::std::max_element( val1.begin(), val1.end() ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}




// template<class fwd_it, class Compare> fwd_it 
// max_element(fwd_it first, fwd_it last, Compare comp);
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<const r,result_tag>::type max_element(const r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::max_element( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
} 
template<range_return_type result_tag, class r, class pred> inline typename detail::range_result<r,result_tag>::type max_element(r& val1, pred val2) { 
    return detail::return_range( val1.begin(), ::std::max_element( val1.begin(), val1.end(), val2 ), val1.end(), 
        detail::range_result_tag<result_tag>() ); 
}






} // namespace rng
}}}}

#endif
/*
Compiles under VC6 !!!

#include <iostream>
typedef enum range_return_type {
    iter,
    from_beg,
    to_end
};


template< class T>
void test( const T& val){
    std::cout << "general" << std::endl;
}

template< range_return_type r, class T>
void test( const T& val, range_return_type ignore_vc_bug = r){
    std::cout << "specific" << std::endl;
}

struct try_it {};
int main(int argc, char* argv[])
{
    test( 1);
    test( try_it() );
    test<it>( 1 );
    test<from_beg>( 2 );
    test<it>( try_it() );
    test<to_end>( 2 );
}
*/
