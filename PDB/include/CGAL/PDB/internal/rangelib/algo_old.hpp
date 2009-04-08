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

    2. All predicates are passed by value, since VC6 does not hadle 'const Pred&' well,
    when they are functions.
    
    3. At this time (14th Nov 2003), all STL algorithms are implemented, except for:
    - partial_sort_copy
    - mismatch


// FIXME - implement <numeric> algorithms as well
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




namespace rng {




template< class r1, class r2>
struct mismatch_return_finder {
    typedef typename range_finder<r1>::range_type r1_range_type;
    typedef typename r1_range_type::iterator it1;
    
    typedef typename r2::iterator it2;
    typedef ::std::pair<it1, it2> return_type;
};



/* 
    Algorithms that DO NOT return a range.
*/

//template<class in_it> typename iterator_traits<in_it>::difference_type
//distance(in_it first, in_it last);
// FIXME - return difference_type
CGAL_PDB_BOOST_RTL_RNG_ALGO10( distance, int)

//template <class in_it, class Distance>
//void advance(in_it& i, Distance n);



template< class iterator, class Distance> inline
iterator advanced( iterator i, Distance n) {
    ::std::advance(i, n);
    return i;
}



//template<class in_it, class Function>
//Function for_each(in_it first, in_it last, Function f);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( for_each, function, function, function)

// template<class in_it, class T> typename iterator_traits<in_it>::difference_type
// count(in_it first, in_it last, const T& value);
//     FIXME - in the future, use difference_type
CGAL_PDB_BOOST_RTL_RNG_ALGO11( count, int, T, const T&)

// template<class in_it, class Predicate> typename iterator_traits<in_it>::difference_type
// count_if(in_it first, in_it last, Predicate pred);
//     FIXME - in the future, use difference_type
CGAL_PDB_BOOST_RTL_RNG_ALGO11( count_if, int, T, T)

// 
// template<class in_it1, class in_it2> pair<in_it1, in_it2> 
// mismatch(in_it1 first1, in_it1 last1, in_it2 first2);
template< class r1, class r2>
struct mismatch_return_finder {
    typedef typename range_finder<r1>::range_type r1_range_type;
    typedef typename r1_range_type::iterator it1;
    // FIXME - if r2 is a const container, use const_iterator!!!
    typedef typename r2::iterator it2;
    typedef ::std::pair<it1, it2> return_type;
};

//FIXME - fails under VC6
//CGAL_PDB_BOOST_RTL_RNG_ALGO20( mismatch, (BOOST_DEDUCED_TYPENAME mismatch_return_finder<r1,r2>::return_type))
mismatch
// 
// template<class in_it1, class in_it2, class BinaryPredicate> pair<in_it1, in_it2> 
// mismatch(in_it1 first1, in_it1 last1, in_it2 first2, BinaryPredicate pred);
//FIXME - fails under VC6
//CGAL_PDB_BOOST_RTL_RNG_ALGO21( mismatch, (BOOST_DEDUCED_TYPENAME mismatch_return_finder<r1,r2>::return_type), pred, pred)



// 
// template<class in_it1, class in_it2> bool 
// equal(in_it1 first1, in_it1 last1, in_it2 first2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20( equal, bool)

// 
// template<class in_it1, class in_it2, class BinaryPredicate> bool 
// equal(in_it1 first1, in_it1 last1, in_it2 first2, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO21( equal, bool, pred, pred)

// 
// template<class in_it, class out_it> out_it 
// copy(in_it first, in_it last, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( copy, out_it, out_it, out_it)




// 
// template<class bid_it1, class bid_it2> bid_it2 
// copy_backward(bid_it1 first, bid_it1 last, bid_it2 result);
//     FIXME - here, we could actually return a range
CGAL_PDB_BOOST_RTL_RNG_ALGO20( copy_backward, BOOST_DEDUCED_TYPENAME r2::iterator)

// 
// template<class fwd_it1, class fwd_it2> fwd_it2 
// swap_ranges(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2);
//     FIXME - here, we could actually return a range
CGAL_PDB_BOOST_RTL_RNG_ALGO20( swap_ranges, BOOST_DEDUCED_TYPENAME r2::iterator)

// 
// template<class fwd_it1, class fwd_it2>
// void iter_swap(fwd_it1 a, fwd_it2 b);
//CGAL_PDB_BOOST_RTL_RNG_ALGO20_void( iter_swap)

// 
// template<class in_it, class out_it, class UnaryOperation> out_it 
// transform(in_it first, in_it last, out_it result, UnaryOperation op);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( transform, out_it, out_it, out_it, op, op)

// 
// template<class in_it1, class in_it2, class out_it, class BinaryOperation> out_it 
// transform(in_it1 first1, in_it1 last1, in_it2 first2, out_it result, BinaryOperation binary_op);
CGAL_PDB_BOOST_RTL_RNG_ALGO22( transform, out_it, out_it, out_it, op, op)

// 
// template<class fwd_it, class T> void 
// replace(fwd_it first, fwd_it last, const T& old_value, const T& new_value);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( replace, T, const T&, T_same, const T_same&)

// 
// template<class fwd_it, class Predicate, class T> void 
// replace_if(fwd_it first, fwd_it last, Predicate pred, const T& new_value);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( replace_if, pred, pred, T, const T&)


// 
// template<class in_it, class out_it, class T> out_it 
// replace_copy(in_it first, in_it last, out_it result, const T& old_value, const T& new_value);
CGAL_PDB_BOOST_RTL_RNG_ALGO13( replace_copy, out_it, out_it, out_it, T, const T&, T_same, const T_same&)

// 
// template<class Iterator, class out_it, class Predicate, class T> out_it 
// replace_copy_if(Iterator first, Iterator last, out_it result, Predicate pred, const T& new_value);
CGAL_PDB_BOOST_RTL_RNG_ALGO13( replace_copy_if, out_it, out_it, out_it, pred, pred, T, const T&)

// 
// template<class fwd_it, class T> void 
// fill(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( fill, T, const T&)

// 
// template<class out_it, class Size, class T> void 
// fill_n(out_it first, Size n, const T& value);
//useless - CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( fill_n, Size, Size, T, const T&)

// 
// template<class fwd_it, class Generator> void 
// generate(fwd_it first, fwd_it last, Generator gen);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( generate, generator, generator)


// 
// template<class out_it, class Size, class Generator> void 
// generate_n(out_it first, Size n, Generator gen);
// useless - CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( generate_n, Size, Size, generator, generator)

// 
// template<class in_it, class out_it> out_it 
// unique_copy(in_it first, in_it last, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( unique_copy, out_it, out_it, out_it)

// 
// template<class in_it, class out_it, class BinaryPredicate> out_it 
// unique_copy(in_it first, in_it last, out_it result, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( unique_copy, out_it, out_it, out_it, pred, pred)

// 
// template<class bid_it> void 
// reverse(bid_it first, bid_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( reverse)

// 
// template<class bid_it, class out_it> out_it 
// reverse_copy(bid_it first, bid_it last, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( reverse_copy, out_it, out_it, out_it)

// 
// template<class fwd_it> void 
// rotate(fwd_it first, fwd_it middle, fwd_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( rotate, fwd_it, fwd_it)

// 
// template<class fwd_it, class out_it> out_it 
// rotate_copy(fwd_it first, fwd_it middle, fwd_it last, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( rotate_coy, out_it, fwd_it, fwd_it, out_it, out_it)

// 
// template<class rand_it> void 
// random_shuffle(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( random_shuffle)

// 
// template<class rand_it, class RandomNumberGenerator> void 
// random_shuffle(rand_it first, rand_it last, RandomNumberGenerator& rand);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( random_shuffle, generator, generator&)

// 
// template<class rand_it> void 
// sort(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( sort)

// 
// template<class rand_it, class Compare>
// void sort(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( sort, pred, pred)

// 
// template<class rand_it> void 
// stable_sort(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( stable_sort)

// 
// template<class rand_it, class Compare> void 
// stable_sort(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( stable_sort, pred, pred)

// 
// template<class rand_it> void 
// partial_sort(rand_it first, rand_it middle, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( partial_sort, rand_it, rand_it)

// 
// template<class rand_it, class Compare> void 
// partial_sort(rand_it first, rand_it middle, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( partial_sort, rand_it, rand_it, pred, pred)

// 
// template<class rand_it> void 
// nth_element(rand_it first, rand_it nth, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( nth_element, rand_it, rand_it)

// 
// template<class rand_it, class Compare> void 
// nth_element(rand_it first, rand_it nth, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( nth_element, rand_it, rand_it, pred, pred)

// 
// template<class fwd_it, class T> bool 
// binary_search(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( binary_search, bool, T, const T&)

// 
// template<class fwd_it, class T, class Compare> bool 
// binary_search(fwd_it first, fwd_it last, const T& value, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( binary_search, bool, T, const T&, pred, pred)

// 
// template<class in_it1, class in_it2, class out_it> out_it 
// merge(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( merge, out_it, out_it, out_it)

// 
// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// merge(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO22_2param( merge, out_it, out_it, out_it, pred, pred)

// 
// template<class bid_it> void 
// inplace_merge(bid_it first, bid_it middle, bid_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( inplace_merge, bid_it, bid_it)

// 
// template<class bid_it, class Compare> void 
// inplace_merge(bid_it first, bid_it middle, bid_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_void( inplace_merge, bid_it, bid_it, pred, pred)


// 
// template<class in_it1, class in_it2> bool 
// includes(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20_2param( includes, bool)

// 
// template<class in_it1, class in_it2, class Compare> bool 
// includes(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( includes, bool, comp, comp)

// 
// template<class in_it1, class in_it2, class out_it> out_it 
// set_union(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( set_union, out_it, out_it, out_it)

// 
// template<class in_it1, class in_it2, class out_it, class Compare>
// out_it set_union(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO22_2param( set_union, out_it, out_it, out_it, pred, pred)

// 
// template<class in_it1, class in_it2, class out_it> out_it 
// set_intersection(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( set_intersection, out_it, out_it, out_it)

// 
// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// set_intersection(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO22_2param( set_intersection, out_it, out_it, out_it, pred, pred)

// 
// template<class in_it1, class in_it2, class out_it> out_it 
// set_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( set_difference, out_it, out_it, out_it)

// 
// template<class in_it1, class in_it2, class out_it, class Compare> out_it 
// set_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO22_2param( set_difference, out_it, out_it, out_it, pred, pred)

// 
// template<class in_it1, class in_it2, class out_it> out_it
// set_symmetric_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( set_symmetric_difference, out_it, out_it, out_it)

// 
// template<class in_it1, class in_it2, class out_it, class Compare> out_it
// set_symmetric_difference(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, out_it result, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO22_2param( set_symmetric_difference, out_it, out_it, out_it, pred, pred)

// 
// 
// template<class rand_it> void 
// push_heap(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( push_heap)

// 
// template<class rand_it, class Compare> void 
// push_heap(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( push_heap, pred, pred)

// 
// template<class rand_it> void 
// pop_heap(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( pop_heap)

// 
// template<class rand_it, class Compare> void 
// pop_heap(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( pop_heap, pred, pred)

// 
// template<class rand_it> void 
// make_heap(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( make_heap)

// 
// template<class rand_it, class Compare> void 
// make_heap(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( make_heap, pred, pred)

// 
// template<class rand_it> void 
// sort_heap(rand_it first, rand_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_void( sort_heap)

// 
// template<class rand_it, class Compare> void 
// sort_heap(rand_it first, rand_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_void( sort_heap, pred, pred)

// 
// template<class in_it1, class in_it2> bool 
// lexicographical_compare(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20_2param( lexicographical_compare, bool)

// 
// template<class in_it1, class in_it2, class Compare> bool 
// lexicographical_compare(in_it1 first1, in_it1 last1, in_it2 first2, in_it2 last2, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_2param( lexicographical_compare, bool, pred, pred)

// 
// template<class bid_it> bool 
// next_permutation(bid_it first, bid_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10( next_permutation, bool)

// 
// template<class bid_it, class Compare> bool 
// next_permutation(bid_it first, bid_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( next_permutation, bool, pred, pred)

// 
// template<class bid_it> bool 
// prev_permutation(bid_it first, bid_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10( prev_permutation, bool)

// 
// template<class bid_it, class Compare> bool 
// prev_permutation(bid_it first, bid_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( prev_permutation, bool, pred, pred)

// template<class in_it, class out_it, class T> out_it 
// remove_copy(in_it first, in_it last, out_it result, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( remove_copy, out_it, out_it, out_it, T, const T&)

// 
// template<class in_it, class out_it, class Predicate> out_it 
// remove_copy_if(in_it first, in_it last, out_it result, Predicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( remove_copy_if, out_it, out_it, out_it, pred, pred)




////////////////////////////////////////////////////////////////////////////
// from <numerics>


// template <class in_it, class T> T 
// accumulate(in_it first, in_it last, T init);
CGAL_PDB_BOOST_RTL_RNG_ALGO11( accumulate, T, T, T)

// template <class in_it, class T, class Func> T 
// accumulate(in_it first, in_it last, T init, Func binary_op);
CGAL_PDB_BOOST_RTL_RNG_ALGO12( accumulate, T, T, T, Func, Func)


// template <class in_it1, class in_it2, class T> T 
// inner_product(in_it1 first1, in_it1 last1, in_it2 first2, T init);

// template <class in_it1, class in_it2, class T, class Func1, class Func2> T 
// inner_product(in_it1 first1, in_it1 last1, in_it2 first2, T init, Func1 binary_op1, Func2 binary_op2);

// template <class in_it, class out_it> out_it 
// partial_sum(in_it first, in_it last, out_it result);

// template <class in_it, class out_it, class Func> out_it 
// partial_sum(in_it first, in_it last, out_it result, Func binary_op);

// template <class in_it, class out_it> out_it 
// adjacent_difference(in_it first, in_it last, out_it result);

// template <class in_it, class out_it, class Func> out_it 
// adjacent_difference(in_it first, in_it last, out_it result, Func binary_op);









// 
// 
// 
// 
// * These functions returns a range

/*
    Right now, we always return a range, from the element "found",
    up to the end of the original range.

    This might not necessary be what you want.

    Therefore, FIXME (to do).

    Allow for algorithm of this type:
    find_if<what_to_return>(range);

    What to return is any of:
    1. return the iterator itself ("r_it")
    2. (default - as is now) return the range from the 
       found iterator up to the end ("r_to_end")
    3. return the range starting from begin, up to the found iterator ("r_from_beg")

    Note: see end of file.
*/
// template<class in_it, class T> in_it 
// find(in_it first, in_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( find, T, const T&)

// 
// template<class in_it, class Predicate> in_it 
// find_if(in_it first, in_it last, Predicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( find_if, pred, pred)

// 
// template<class fwd_it1, class fwd_it2> fwd_it1 
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20_EXT( find_end)

// 
// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_end(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_EXT( find_end, pred, pred)

// 
// template<class fwd_it1, class fwd_it2> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20_EXT( find_first_of)

// 
// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1
// find_first_of(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_EXT( find_first_of, pred, pred)

// 
// template<class fwd_it> fwd_it 
// adjacent_find(fwd_it first, fwd_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_EXT( adjacent_find)

// template<class fwd_it, class BinaryPredicate> fwd_it 
// adjacent_find(fwd_it first, fwd_it last, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( adjacent_find, pred, pred)

// 
// 
// template<class fwd_it1, class fwd_it2> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2);
CGAL_PDB_BOOST_RTL_RNG_ALGO20_EXT( search)

// 
// template<class fwd_it1, class fwd_it2, class BinaryPredicate> fwd_it1 
// search(fwd_it1 first1, fwd_it1 last1, fwd_it2 first2, fwd_it2 last2, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO21_EXT( search, pred, pred)

// 
// template<class fwd_it, class Size, class T> fwd_it 
// search_n(fwd_it first, fwd_it last, Size count, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_EXT( search_n, Size, Size, T, const T&)

// 
// template<class fwd_it, class Size, class T, class BinaryPredicate> fwd_it1 
// search_n(fwd_it first, fwd_it last, Size count, const T& value, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO13_EXT( search_n, Size, Size, T, const T&, pred, pred)

// 
// ***** return the range of elements to be removed;
//        use them in correlation with r_erase.
//

// template<class fwd_it, class T> fwd_it 
// remove(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( remove, T, const T&)

// 
// template<class fwd_it, class Predicate> fwd_it 
// remove_if(fwd_it first, fwd_it last, Predicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( remove_if, pred, pred)

// 
// template<class fwd_it> fwd_it 
// unique(fwd_it first, fwd_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_EXT( unique)

// 
// template<class fwd_it, class BinaryPredicate> fwd_it 
// unique(fwd_it first, fwd_it last, BinaryPredicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( unique, pred, pred)

// 
// 
// // return the range that satisfies the predicate
// template<class bid_it, class Predicate> bid_it 
// partition(bid_it first, bid_it last, Predicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( partition, pred, pred)

// 
// template<class bid_it, class Predicate> bid_it 
// stable_partition(bid_it first, bid_it last, Predicate pred);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( stable_partition, pred, pred)

// 
// template<class in_it, class rand_it> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last);
// FIXME - watch out what we return!!!

// 
// template<class in_it, class rand_it, class Compare> rand_it
// partial_sort_copy(in_it first, in_it last, rand_it result_first, rand_it result_last, Compare comp);
// FIXME - watch out what we return!!!


// 
// template<class fwd_it, class T> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( lower_bound, T, const T&)

// 
// template<class fwd_it, class T, class Compare> fwd_it 
// lower_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_EXT( lower_bound, T, const T&, pred, pred)

// 
// template<class fwd_it, class T> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( upper_bound, T, const T&)

// 
// template<class fwd_it, class T, class Compare> fwd_it 
// upper_bound(fwd_it first, fwd_it last, const T& value, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_EXT( upper_bound, T, const T&, pred, pred)

// 
// template<class fwd_it, class T> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT2( equal_range, T, const T&)


// 
// template<class fwd_it, class T, class Compare> pair<fwd_it, fwd_it>
// equal_range(fwd_it first, fwd_it last, const T& value, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO12_EXT2( equal_range, T, const T&, pred, pred)


// 
// template<class fwd_it> fwd_it 
// min_element(fwd_it first, fwd_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_EXT( min_element)

// 
// template<class fwd_it, class Compare> fwd_it 
// min_element(fwd_it first, fwd_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( min_element, pred, pred)

// 
// template<class fwd_it> fwd_it 
// max_element(fwd_it first, fwd_it last);
CGAL_PDB_BOOST_RTL_RNG_ALGO10_EXT( max_element)

// 
// template<class fwd_it, class Compare> fwd_it 
// max_element(fwd_it first, fwd_it last, Compare comp);
CGAL_PDB_BOOST_RTL_RNG_ALGO11_EXT( max_element, pred, pred)






#if 0


I just thought about it, and one could actually say:
r.begin( v.erase(r.begin()) );

So I could add a function called errase_current which would implement this idiom. Does that make sense?
(and it could work for collections as well)


must be done - already written in article listing


 FIXME
Hi John and Mathew,

Amid my exam confusion, I came to think on this example.

In normal iterator mode we have

typedef vector<int> vec_t;
vec_t v;
for( vec_t::iterator i = v.begin(); i != v.end(); ++i )
   if( something( *i ) )
     i = v.erase( i );

notice that we need to recompute the end iterator since erase can invalidate it.

How can we do that with ranges? I guess this would be wrong:

crange<vec_t> r( v );
for( ; r ; ++r ) // not good
  if( something(*r) )
    r = v.erase( r.begin() );

Do your free-stading erase work:

crange<vec_t> r( v );
for( ; r ; ++r ) /
  if( something(*r) )
    r = erase( v, r );

? ie, is erase defined as

template< class C >
crange<C> erase( C& c, const crange<C>& r )
{
   returrn crange<C>( c.erase( r.begin() ), c.end() );
}

?

btw, I checked the problem with retur values of map/set erase() being void. It turns out that we decided it was a defect 
in the standard that the return value is void; hence standard libraries are requiered to change the implementation ASAP, ie, before
C++Ox. I think some libraries already has a non-void return.


br

Thorsten




































As was mentioned in the first article, algorithms that return an iterator are indeed a different breed from the rest. You use such an algorithm to:
	1. do something to the iterator you get in return, or
	2. walk from that iterator to the end or
	3. walk from the beginning up to the iterator 

In the first two cases, you'll typically be happy with the default range algorithm behavior - you do need to test if the iterator is valid, and if so, use it (or walk to the end). In case you want the algorithm to:
	- return only the iterator, and behave like the STL algorithm: use rng::algorithm<iter>(…) instead of rng::algorithm(…)
 	return the [begin, found) range (case 3.): use rng::algorithm<from_beg>(…) instead of rng::algorithm(…)
[for searching and find, there might be a couple of other policies: maybe it could also be possible to say
vector<int> v;
find<or_throw>( v, 3 ); // throws an exception on failure
int i = find<or_return_default>( v, is_even() ); // finds first even number or 0 otherwise
]






[Another issue might be to have different overloads for eg. Copy: copy( range, iterator ); copy( range, range);] 
This is already solved, just mention it in the docs.

#endif

















} // namespace rng
}}}}

#endif
/*
Compiles under VC6 !!!

#include <iostream>
typedef enum range_return_type {
    it,
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
