//  Boost string_algo library erase.hpp header file  ---------------------------//

//  Copyright Pavol Droba 2002-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_STRING_ERASE_HPP
#define BOOST_STRING_ERASE_HPP

#include <boost/algorithm/string/config.hpp>
#include <boost/algorithm/string/collection_traits.hpp>
#include <boost/algorithm/string/iterator_range.hpp>
#include <boost/algorithm/string/find_format.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/formatter.hpp>

/*! \file
    Defines various erase algorithms. Each algorithm removes
    part(s) of the input according to a searching criteria.
*/

namespace boost {
    namespace algorithm {

//  erase_range -------------------------------------------------------//

        //! Erase range algorithm
        /*!
            Remove the given range from the input. The result is a modified copy of 
            the input. It is returned as a sequence or copied to the output iterator.
    
            \param Output An output iterator to which the result will be copied
            \param Input An input sequence
            \param SearchRange A range in the input to be removed
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<typename OutputIteratorT, typename CollectionT>
        inline OutputIteratorT erase_range_copy(
            OutputIteratorT Output,
            const CollectionT& Input,
            const iterator_range<
                BOOST_STRING_TYPENAME 
                    const_iterator_of<CollectionT>::type>& SearchRange )
        {
            return find_format_copy(
                Output,
                Input,
                range_finder(SearchRange),
                empty_formatter(Input) );
        }

        //! Erase range algorithm
        /*!
            \overload
        */
        template<typename SequenceT>
        inline SequenceT erase_range_copy( 
            const SequenceT& Input,
            const iterator_range<
                BOOST_STRING_TYPENAME 
                    const_iterator_of<SequenceT>::type>& SearchRange )
        {
            return find_format_copy( 
                Input,
                range_finder(SearchRange),
                empty_formatter(Input) );
        }

        //! Erase range algorithm
        /*!
            Remove the given range from the input.
            The input sequence is modified in-place.

            \param Input An input sequence
            \param SearchRange A range in the input to be removed
        */
        template<typename SequenceT>
        inline void erase_range( 
            SequenceT& Input,
            const iterator_range<
                BOOST_STRING_TYPENAME 
                    iterator_of<SequenceT>::type>& SearchRange )
        {
            find_format( 
                Input, 
                range_finder(SearchRange),
                empty_formatter(Input) );
        }

//  erase_first  --------------------------------------------------------//

        //! Erase first algorithm
        /*!
            Remove the first occurence of the substring from the input.
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for 
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input
            
            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT erase_first_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search )
        {
            return find_format_copy(
                Output,
                Input,
                first_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase first algorithm
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT erase_first_copy( 
            const SequenceT& Input,
            const CollectionT& Search )
        {
            return find_format_copy( 
                Input, 
                first_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase first algorithm
        /*!
            Remove the first occurence of the substring from the input. 
            The input sequence is modified in-place.

            \param Input An input string
            \param Search A substring to be searched for. 
        */
        template<typename SequenceT, typename CollectionT>
        inline void erase_first( 
            SequenceT& Input,
            const CollectionT& Search )
        {
            find_format( 
                Input, 
                first_finder(Search),
                empty_formatter(Input) );
        }

//  erase_first ( case insensitive ) ------------------------------------//

        //! Erase first algorithm ( case insensitive )
        /*!
            Remove the first occurence of the substring from the input. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.
            Searching is case insensitive.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for 
            \param Loc A locale used for case insensitive comparison
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT ierase_first_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy(
                Output,
                Input,
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase first algorithm ( case insensitive )
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT ierase_first_copy( 
            const SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy( 
                Input, 
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase first algorithm ( case insensitive )
        /*!
            Remove the first occurence of the substring from the input. 
            The input sequence is modified in-place. Searching is case insensitive.

            \param Input An input string
            \param Search A substring to be searched for
            \param Loc A locale used for case insensitive comparison
        */
        template<typename SequenceT, typename CollectionT>
        inline void ierase_first( 
            SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            find_format( 
                Input, 
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

//  erase_last  --------------------------------------------------------//

        //! Erase last algorithm
        /*!
            Remove the last occurence of the substring from the input. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for.
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

             \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT erase_last_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search )
        {
            return find_format_copy(
                Output,
                Input,
                last_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase last algorithm
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT erase_last_copy( 
            const SequenceT& Input,
            const CollectionT& Search )
        {
            return find_format_copy( 
                Input, 
                last_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase last algorithm
        /*!
            Remove the last occurence of the substring from the input. 
            The input sequence is modified in-place.

            \param Input An input string
            \param Search A substring to be searched for 
        */
        template<typename SequenceT, typename CollectionT>
        inline void erase_last( 
            SequenceT& Input,
            const CollectionT& Search )
        {
            find_format( 
                Input, 
                last_finder(Search),
                empty_formatter(Input) );
        }

//  erase_last ( case insensitive ) ------------------------------------//

        //! Erase last algorithm ( case insensitive )
        /*!
            Remove the last occurence of the substring from the input. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.
            Searching is case insensitive.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for
            \param Loc A locale used for case insensitive comparison
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

             \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT ierase_last_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy(
                Output,
                Input,
                last_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase last algorithm ( case insensitive )
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT ierase_last_copy( 
            const SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy( 
                Input, 
                last_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase last algorithm ( case insensitive )
        /*!
            Remove the last occurence of the substring from the input. 
            The input sequence is modified in-place. Searching is case insensitive.

            \param Input An input string
            \param Search A substring to be searched for
            \param Loc A locale used for case insensitive comparison
        */
        template<typename SequenceT, typename CollectionT>
        inline void ierase_last( 
            SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            find_format( 
                Input, 
                last_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

//  erase_nth --------------------------------------------------------------------//

        //! Erase nth algorithm
        /*!
            Remove the Nth occurence of the substring in the input.
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.
            

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for
            \param Nth An index of the match to be replaced. The index is 0-based.
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

             \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT erase_nth_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search,
            unsigned int Nth )
        {
            return find_format_copy(
                Output,
                Input,
                nth_finder(Search, Nth),
                empty_formatter(Input) );
        }

        //! Erase nth algorithm
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT erase_nth_copy( 
            const SequenceT& Input,
            const CollectionT& Search,
            unsigned int Nth )
        {
            return find_format_copy( 
                Input, 
                nth_finder(Search, Nth),
                empty_formatter(Input) );
        }

        //! Erase nth algorithm
        /*!
            Remove the Nth occurence of the substring in the input.
            The input sequence is modified in-place.

            \param Input An input string
            \param Search A substring to be searched for. 
            \param Nth An index of the match to be replaced. The index is 0-based.
        */
        template<typename SequenceT, typename CollectionT>
        inline void erase_nth( 
            SequenceT& Input,
            const CollectionT& Search,
            unsigned int Nth )
        {
            find_format( 
                Input, 
                nth_finder(Search, Nth),
                empty_formatter(Input) );
        }

//  erase_nth ( case insensitive ) ---------------------------------------------//

        //! Erase nth algorithm ( case insensitive )
        /*!
            Remove the Nth occurence of the substring in the input.
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator. 
            Searching is case insensitive.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for.
            \param Nth An index of the match to be replaced. The index is 0-based.
            \param Loc A locale used for case insensitive comparison
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT ierase_nth_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search,
            unsigned int Nth,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy(
                Output,
                Input,
                nth_finder(Search, Nth, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase nth algorithm
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT ierase_nth_copy( 
            const SequenceT& Input,
            const CollectionT& Search,
            unsigned int Nth,
            const std::locale& Loc=std::locale() )
        {
            return find_format_copy( 
                Input, 
                nth_finder(Search, Nth, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase nth algorithm
        /*!
            Remove the Nth occurence of the substring in the input.
            The input sequence is modified in-place. Searching is case insensitive.

            \param Input An input string
            \param Search A substring to be searched for. 
            \param Nth An index of the match to be replaced. The index is 0-based.
            \param Loc A locale used for case insensitive comparison
        */
        template<typename SequenceT, typename CollectionT>
        inline void ierase_nth( 
            SequenceT& Input,
            const CollectionT& Search,
            unsigned int Nth,
            const std::locale& Loc=std::locale() )
        {
            find_format( 
                Input, 
                nth_finder(Search, Nth, is_iequal(Loc)),
                empty_formatter(Input) );
        }


//  erase_all  --------------------------------------------------------//

        //! Erase all algorithm
        /*!
            Remove all the occurrences of the string from the input. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.
                        

            \param Output An output iterator to which the result will be copied
            \param Input An input sequence
            \param Search A substring to be searched for. 
            \return An output iterator pointing just after the last inserted character or
                    a modified copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT erase_all_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search )
        {
            return find_format_all_copy(
                Output,
                Input,
                first_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase all algorithm
        /*!
            \overload
        */  
        template<typename SequenceT, typename CollectionT>
        inline SequenceT erase_all_copy( 
            const SequenceT& Input,
            const CollectionT& Search )
        {
            return find_format_all_copy( 
                Input, 
                first_finder(Search),
                empty_formatter(Input) );
        }

        //! Erase all algorithm
        /*!
            Remove all the occurrences of the string from the input. 
            The input sequence is modified in-place.

            \param Input An input string
            \param Search A substring to be searched for. 
        */
        template<typename SequenceT, typename CollectionT>
        inline void erase_all( 
            SequenceT& Input,
            const CollectionT& Search )
        {
            find_format_all( 
                Input, 
                first_finder(Search),
                empty_formatter(Input) );
        }

//  erase_all ( case insensitive ) ------------------------------------//

        //! Erase all algorithm ( case insensitive )
        /*!
            Remove all the occurrences of the string from the input. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator. 
            Searching is case insensitive.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param Search A substring to be searched for
            \param Loc A locale used for case insensitive comparison
            \return An output iterator pointing just after the last inserted character or
                    a modified copy of the input

              \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename Collection1T, 
            typename Collection2T>
        inline OutputIteratorT ierase_all_copy(
            OutputIteratorT Output,
            const Collection1T& Input,
            const Collection2T& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_all_copy(
                Output,
                Input,
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase all algorithm ( case insensitive )
        /*!
            \overload
        */
        template<typename SequenceT, typename CollectionT>
        inline SequenceT ierase_all_copy( 
            const SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            return find_format_all_copy( 
                Input, 
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

        //! Erase all algorithm ( case insensitive )
        /*!
            Remove all the occurrences of the string from the input. 
            The input sequence is modified in-place. Searching is case insensitive.

            \param Input An input string
            \param Search A substring to be searched for. 
            \param Loc A locale used for case insensitive comparison
        */
        template<typename SequenceT, typename CollectionT>
        inline void ierase_all( 
            SequenceT& Input,
            const CollectionT& Search,
            const std::locale& Loc=std::locale() )
        {
            find_format_all( 
                Input, 
                first_finder(Search, is_iequal(Loc)),
                empty_formatter(Input) );
        }

//  erase_head --------------------------------------------------------------------//

        //! Erase head algorithm
        /*!
            Remove the head from the input. The head is a prefix of a seqence of given size. 
            If the sequence is shorter then required, the whole string is 
            considered to be the head. The result is a modified copy of the input. 
            It is returned as a sequence or copied to the output iterator.
            

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param N Length of the head
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input

             \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename CollectionT>
        inline OutputIteratorT erase_head_copy(
            OutputIteratorT Output,
            const CollectionT& Input,
            unsigned int N )
        {
            return find_format_copy(
                Output,
                Input,
                head_finder(N),
                empty_formatter( Input ) );
        }

        //! Erase head algorithm
        /*!
            \overload
        */
        template<typename SequenceT>
        inline SequenceT erase_head_copy( 
            const SequenceT& Input,
            unsigned int N )
        {
            return find_format_copy( 
                Input,
                head_finder(N),
                empty_formatter( Input ) );
        }

        //! Erase head algorithm
        /*!
            Remove the head from the input. The head is a prefix of a seqence of given size. 
            If the sequence is shorter then required, the whole string is 
            considered to be the head. The input sequence is modified in-place.

            \param Input An input string
            \param N Length of the head
        */
        template<typename SequenceT>
        inline void erase_head( 
            SequenceT& Input,
            unsigned int N )
        {
            find_format( 
                Input, 
                head_finder(N),
                empty_formatter( Input ) );
        }

//  erase_tail --------------------------------------------------------------------//

        //! Erase tail algorithm
        /*!
            Remove the tail from the input. The tail is a suffix of a seqence of given size. 
            If the sequence is shorter then required, the whole string is 
            considered to be the tail. 
            The result is a modified copy of the input. It is returned as a sequence 
            or copied to the output iterator.

            \param Output An output iterator to which the result will be copied
            \param Input An input string
            \param N Length of the head
            \return An output iterator pointing just after the last inserted character or
                a modified copy of the input
            
             \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<
            typename OutputIteratorT,
            typename CollectionT>
        inline OutputIteratorT erase_tail_copy(
            OutputIteratorT Output,
            const CollectionT& Input,
            unsigned int N )
        {
            return find_format_copy(
                Output,
                Input,
                tail_finder(N),
                empty_formatter( Input ) );
        }

        //! Erase tail algorithm
        /*!
            \overload
        */
        template<typename SequenceT>
        inline SequenceT erase_tail_copy( 
            const SequenceT& Input,
            unsigned int N )
        {
            return find_format_copy( 
                Input,
                tail_finder(N),
                empty_formatter( Input ) );
        }

        //! Erase tail algorithm
        /*!
            Remove the tail from the input. The tail is a suffix of a seqence of given size. 
            If the sequence is shorter then required, the whole string is
            considered to be the tail. The input sequence is modified in-place.

            \param Input An input string
            \param N Length of the head
        */
        template<typename SequenceT>
        inline void erase_tail( 
            SequenceT& Input,
            unsigned int N )
        {
            find_format( 
                Input, 
                tail_finder(N),
                empty_formatter( Input ) );
        }

    } // namespace algorithm

    // pull names into the boost namespace
    using algorithm::erase_range_copy;
    using algorithm::erase_range;
    using algorithm::erase_first_copy;
    using algorithm::erase_first;
    using algorithm::ierase_first_copy;
    using algorithm::ierase_first;
    using algorithm::erase_last_copy;
    using algorithm::erase_last;
    using algorithm::ierase_last_copy;
    using algorithm::ierase_last;
    using algorithm::erase_nth_copy;
    using algorithm::erase_nth;
    using algorithm::ierase_nth_copy;
    using algorithm::ierase_nth;
    using algorithm::erase_all_copy;
    using algorithm::erase_all;
    using algorithm::ierase_all_copy;
    using algorithm::ierase_all;
    using algorithm::erase_head_copy;
    using algorithm::erase_head;
    using algorithm::erase_tail_copy;
    using algorithm::erase_tail;

} // namespace boost


#endif  // BOOST_ERASE_HPP
