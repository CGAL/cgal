//  Boost string_algo library iterator_range.hpp header file  ---------------------------//

//  Copyright Pavol Droba 2002-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_STRING_ITERATOR_RANGE_HPP
#define BOOST_STRING_ITERATOR_RANGE_HPP

#include <boost/algorithm/string/config.hpp>
#include <utility>
#include <iterator>
#include <algorithm>
#include <ostream>
#include <boost/detail/iterator.hpp>

/*! \file
    Defines the \c iterator_class and related functions. 
    \c iterator_range is a simple wrapper of the iterator pair idiom. It provides
    a rich subset of the Container interface.
*/

namespace boost {
    namespace algorithm {

//  iterator range template class -----------------------------------------//

        //! iterator_range class
        /*!
            An \c iterator_range delimits a range in a sequence by beginning and ending iterators. 
            An iterator_range can be passed to an algorithm which requires a sequence as an input. 
            For example, the \c toupper() function may most frequently be used on strings, 
            but can also be used on iterator_ranges: 
            
            \code
                boost::tolower( find( s, "UPPERCASE STRING" ) );
            \endcode

            Many algorithms working with sequences take a pair of iterators, 
            delimiting a working range, as arguments. The \c iterator_range class is an 
            encapsulation of a range identified by a pair of iterators. 
            It provides a collection interface, 
            so it is possible to pass an instance to an algorithm requiring a collection as an input. 
        */
        template<typename IteratorT> 
        class iterator_range
        {
        public:
            //! this type
            typedef iterator_range<IteratorT> type;
            //! Encapsulated value type
            typedef BOOST_STRING_TYPENAME boost::detail::
                iterator_traits<IteratorT>::value_type value_type;
            //! Reference type
            typedef BOOST_STRING_TYPENAME boost::detail::
                iterator_traits<IteratorT>::reference reference;
            //! Difference type
            typedef BOOST_STRING_TYPENAME boost::detail::
                iterator_traits<IteratorT>::difference_type difference_type;
            //! Size type
            typedef BOOST_STRING_TYPENAME boost::detail::
                iterator_traits<IteratorT>::difference_type size_type;
            
            //! const_iterator type
            /*! 
                There is no distinction between const_iterator and iterator.
                These typedefs are provides to fulfill container interface
            */ 
            typedef IteratorT const_iterator;
            //! iterator type
            typedef IteratorT iterator;

            //! Empty constructor
            iterator_range() : 
                m_Begin(), m_End() {}

            //! Constructor from a pair of iterators
            iterator_range( iterator Begin, iterator End ) : 
                m_Begin(Begin), m_End(End) {}

            //! Constructor from a std::pair
            iterator_range( const std::pair<IteratorT,IteratorT>& Range ) : 
                m_Begin(Range.first), m_End(Range.second) {}

            //! Copy constructor
            iterator_range( const iterator_range& Other ) :
                m_Begin(Other.begin()), m_End(Other.end()) {}

            //! Templated copy constructor
            /*!
                This constructor is provided to allow conversion between
                const and mutable iterator instances of this class template
            */
            template< typename OtherItT >
            iterator_range( const iterator_range<OtherItT>& Other ) :
                m_Begin(Other.begin()), m_End(Other.end()) {}

            //! Assignment operator
            iterator_range& operator=( const iterator_range& Other )
            {
                m_Begin=Other.begin(); m_End=Other.end();
                return *this;
            }

            //! Assignment operator ( templated version )
            template< typename OtherItT >
            iterator_range& operator=( const iterator_range<OtherItT>& Other )
            {
                m_Begin=Other.begin(); m_End=Other.end();
                return *this;
            }
            
            //! Comparison operator ( equal )
            /*! 
                Compare operands for equality
            */
            template< typename OtherItT > 
            bool operator==( const iterator_range<OtherItT>& Other ) const
            {
                return m_Begin==Other.begin() && m_End==Other.end();
            }

            //! Comparison operator ( not-equal )
            /*! 
                Compare operands for non-equality
            */
            template< typename OtherItT > 
            bool operator!=( const iterator_range<OtherItT>& Other ) const
            {
                return m_Begin!=Other.begin() || m_End!=Other.end();
            }

            //! begin access
            /*! 
                Retrieve the begin iterator
            */
            IteratorT begin() const 
            { 
                return m_Begin; 
            }

            //! end access
            /*!
                Retrieve the end iterator
            */
            IteratorT end() const 
            { 
                return m_End; 
            } 

            //! Empty container test
            /*!
                Test whether the range is empty
            */
            bool empty() const 
            { 
                return m_Begin==m_End; 
            }
            
            //! Size of the range
            /*!
                Retrieve the size of the range
            */
            difference_type size() const
            { 
                return std::distance( m_Begin, m_End ); 
            }

            //! Swap
            /*!
                Swap two ranges
            */
            void swap( iterator_range& Other )
            {
                std::swap( m_Begin, Other.begin() );
                std::swap( m_End, Other.end() );
            }
            
            //! Safe bool conversion
            /*!
                Check whether the range is empty.
                Allows to use construction like this:
                \code
                    iterator_range r;
                    if (!r)
                    {
                    ...
                    }
                \endcode
            */
            typedef iterator (iterator_range::*unspecified_bool_type) () const;
            operator unspecified_bool_type() const
            {
                return empty()? 0: &iterator_range::end;
            }

        private:
            // begin and end iterators
            IteratorT m_Begin;
            IteratorT m_End;
        };

//  iterator range free-standing operators ---------------------------//

        //! iterator_range output operator
        /*!
            Output the range to an ostream. Elements are outputed
            in a sequence without separators.
        */
        template< typename IteratorT, typename Elem, typename Traits >
        std::basic_ostream<Elem,Traits>& operator<<( 
            std::basic_ostream<Elem, Traits>& Os,
            const iterator_range<IteratorT>& Range )
        {
            std::copy(Range.begin(), Range.end(), std::ostream_iterator<Elem>(Os));

            return Os;
        }


//  iterator range utilities -----------------------------------------//

        //! iterator_range construct helper 
        /*!
            Construct an \c iterator_range from a pair of iterators

            \param Begin A begin iterator
            \param End An end iterator
            \return iterator_range object
        */
        template< typename IteratorT >
        inline iterator_range< IteratorT > make_iterator_range( IteratorT Begin, IteratorT End ) 
        {   
            return iterator_range<IteratorT>( Begin, End );
        }
      
        //! iterator_range construct helper
        /*!
            Construct an \c iterator_range from a \c std::pair<> containing the begin
            and end iterators.

            \param Pair A \c std::pair<> with begin and end iterators
            \return \c iterator_range object
        */
        template< typename IteratorT >
        inline iterator_range< IteratorT > make_iterator_range( const std::pair<IteratorT,IteratorT>& Pair ) 
        {   
            return iterator_range<IteratorT>( Pair.first, Pair.second );
        }

        //! copy a range into a sequence
        /*!
            Construct a new sequence of the specified type from the elements
            in the given range

            \param Range An input range
            \return New sequence
        */
        template< typename SeqT, typename IteratorT >
        inline SeqT copy_iterator_range( const iterator_range<IteratorT>& Range )
        {
            return SeqT( Range.begin(), Range.end() );
        }

        //! transform a range into a sequence
        /*!
            Create a new sequence from the elements in the range, transformed
            by a function

            \param Range An input range
            \param Func Transformation function
            \return New sequence
        */
        template< typename SeqT, typename IteratorT, typename FuncT >
        inline SeqT transform_iterator_range( const iterator_range<IteratorT>& Range, FuncT Func )
        {
            SeqT Seq;
            std::transform( Range.begin(), Range.end(), std::back_inserter(Seq), Func );
            return Seq;
        }

    } // namespace algorithm
    
    // pull names to the namespace boost
    using algorithm::iterator_range;
    using algorithm::make_iterator_range;
    using algorithm::copy_iterator_range;
    using algorithm::transform_iterator_range;

} // namespace boost


#endif  // BOOST_STRING_ITERATOR_RANGE_HPP
