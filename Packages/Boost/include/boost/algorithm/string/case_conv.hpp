//  Boost string_algo library case_conv.hpp header file  ---------------------------//

//  Copyright Pavol Droba 2002-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_STRING_CASE_CONV_HPP
#define BOOST_STRING_CASE_CONV_HPP

#include <boost/algorithm/string/config.hpp>
#include <algorithm>
#include <locale>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/algorithm/string/collection_traits.hpp>
#include <boost/algorithm/string/detail/case_conv.hpp>

/*! \file
    Defines sequence case-conversion algorithms.
    Algorithms convert each element in the input sequence to the
    desired case using provided locales.
*/

namespace boost {
    namespace algorithm {

//  to_lower  -----------------------------------------------//

        //! Convert to lower case
        /*!
            Each element of the input sequence is converted to lower
            case. The result is a copy of the input converted to lower case.
            It is returned as a sequence or copied to the output iterator.

            \param Output An output iterator to which the result will be copied
            \param Input An input collection
            \param Loc A locale used for conversion
            \return 
                An output iterator pointing just after the last inserted character or
                a copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
                
        */
        template<typename OutputIteratorT, typename CollectionT>
        inline OutputIteratorT 
        to_lower_copy(
            OutputIteratorT Output,
            const CollectionT& Input,
            const std::locale& Loc=std::locale())
        {
            return std::transform( 
                begin(Input), 
                end(Input), 
                Output,
                detail::to_lowerF<
                    typename value_type_of<CollectionT>::type >(Loc));
        }

        //! Convert to lower case
        /*!
            \overload
        */
        template<typename SequenceT>
        inline SequenceT to_lower_copy( 
            const SequenceT& Input, 
            const std::locale& Loc=std::locale())
        {
            return SequenceT(
                make_transform_iterator(
                    begin(Input),
                    detail::to_lowerF<
                        typename value_type_of<SequenceT>::type >(Loc)),
                make_transform_iterator(
                    end(Input), 
                    detail::to_lowerF<
                        typename value_type_of<SequenceT>::type >(Loc)));
        }

        //! Convert to lower case
        /*!
            Each element of the input sequence is converted to lower
            case. The input sequence is modified in-place.

            \param Input A collection
            \param Loc a locale used for conversion
        */
        template<typename MutableCollectionT>
        inline void to_lower( 
            MutableCollectionT& Input, 
            const std::locale& Loc=std::locale())
        {
            std::transform( 
                begin(Input), 
                end(Input), 
                begin(Input), 
                detail::to_lowerF<
                    typename value_type_of<MutableCollectionT>::type >(Loc));
        }
        
//  to_upper  -----------------------------------------------//

        //! Convert to upper case
        /*!
            Each element of the input sequence is converted to upper
            case. The result is a copy of the input converted to upper case.
            It is returned as a sequence or copied to the output iterator

            \param Output An output iterator to which the result will be copied
            \param Input An input collection
            \param Loc A locale used for conversion
            \return 
                An output iterator pointing just after the last inserted character or
                a copy of the input

            \note The second variant of this function provides the strong exception-safety guarantee
        */
        template<typename OutputIteratorT, typename CollectionT>
        inline OutputIteratorT 
        to_upper_copy(
            OutputIteratorT Output,
            const CollectionT& Input,
            const std::locale& Loc=std::locale())
        {
            return std::transform( 
                begin(Input), 
                end(Input), 
                Output,
                detail::to_upperF<
                    typename value_type_of<CollectionT>::type >(Loc));
        }

        //! Convert to upper case
        /*!
            \overload
        */
        template<typename SequenceT>
        inline SequenceT to_upper_copy( 
            const SequenceT& Input, 
            const std::locale& Loc=std::locale())
        {
            return SequenceT(
                make_transform_iterator(
                    begin(Input),
                    detail::to_upperF<
                        typename value_type_of<SequenceT>::type >(Loc)),
                make_transform_iterator(
                    end(Input), 
                    detail::to_upperF<
                        typename value_type_of<SequenceT>::type >(Loc)));

        }

        //! Convert to upper case
        /*!
            Each element of the input sequence is converted to upper
            case. The input sequence is modified in-place.

            \param Input An input collection
            \param Loc a locale used for conversion
        */
        template<typename MutableCollectionT>
        inline void to_upper( 
            MutableCollectionT& Input, 
            const std::locale& Loc=std::locale())
        {
            std::transform( 
                begin(Input), 
                end(Input), 
                begin(Input), 
                detail::to_upperF<
                    typename value_type_of<MutableCollectionT>::type >(Loc));
        }

    } // namespace algorithm

    // pull names to the boost namespace
    using algorithm::to_lower;
    using algorithm::to_lower_copy;
    using algorithm::to_upper;
    using algorithm::to_upper_copy;

} // namespace boost

#endif  // BOOST_STRING_CASE_CONV_HPP
