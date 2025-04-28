/**
 * \file Bitboard.hpp
 * \brief Arbitrary long bitboard objects.
 * \author Fedyna K.
 * \version 0.2.0
 * \date 23/05/2024
 * 
 * Define everything for the Bitboard class.
 * Bitboards are used with SparseMatrix to quickly find and set non empty columns.
 */

#ifndef __OSM_BITBOARD_HPP__
#define __OSM_BITBOARD_HPP__


#include "__base.hpp"
#include <iterator>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <ostream>
#include <string>
#include <limits>
const size_t size_t_maxvalue = std::numeric_limits<size_t>::max() ;


#define BITBOARD_INT_SIZE 64
#define VECTOR_SIZE(size) size / BITBOARD_INT_SIZE + (size % BITBOARD_INT_SIZE != 0)




/** \brief Table used for fast 64bit log2. */
const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5
};

/** \brief Table used for De Bruijn multiplication. */
const int index64[64] = {
    0, 47,  1, 56, 48, 27,  2, 60,
   57, 49, 41, 37, 28, 16,  3, 61,
   54, 58, 35, 52, 50, 42, 21, 44,
   38, 32, 29, 23, 17, 11,  4, 62,
   46, 55, 26, 59, 40, 36, 15, 53,
   34, 51, 20, 43, 31, 22, 10, 45,
   25, 39, 14, 33, 19, 30,  9, 24,
   13, 18,  8, 12,  7,  6,  5, 63
};


namespace OSM {

/**
 * \class Bitboard
 * \brief Bitboards are long bitset with specific operations.
 * 
 * Bitboards are used with SparseMatrix to quickly find and set non empty columns.
 * 
 * \author Fedyna K.
 * \version 0.2.0
 * \date 23/05/2024
 */
class Bitboard {
private:
    /** \brief The inner data storage (a contiguous array of 64 unsigned bit integers). */
    std::vector<std::uint64_t> boards;
    
    /** \brief The number of bits we want to store. */
    std::size_t _size;
    
public:
    /**
     * \brief Compute fast 64bit int log2.
     * \pre The value must have a population count of 1.
     *
     * \param[in] _value The value we want to get the bit index (i.e. log2).
     * \returns The bit position.
     */
    static std::size_t bitIndex(std::uint64_t _value) {
        if (_value == 1) return 0;
        return tab64[((std::uint64_t)((_value - (_value >> 1)) * 0x07EDD5E59A4E28C2)) >> 58] + 1;
    }
    
    /**
     * \struct Bitboard::iterator
     * \brief The iterator over bitboards.
     *
     * Return all indexes with a bit set to 1.
     *
     * \note The iterator is constant.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    struct iterator {
        /** C++ mandatory type definitions. */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = std::uint64_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        
        /**
         * \brief Iterator constructor.
         * \param[in] _index The initial index.
         * \param[in] _size The bitboard size.
         * \param[in] _boards The reference to the bitboards.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        iterator(const std::size_t _index, const std::size_t _size, const std::vector<std::uint64_t> &_boards):
        index(_index),
        size(_size),
        iteratedBoards(_boards),
        arrayIndex(0)
        {
            if (index == 0) ++(*this);
        }
        
        /**
         * \brief Iterator dereference.
         * \returns A no-null index on the billboard or past-the-end index.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        const std::size_t operator*() const { return index; }
        
        /**
         * \brief Prefix incrementation. Finds the next not-null index.
         * \returns The reference to the current iterator.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        iterator& operator++() {
            // The iterator rely on bit manipulation that are quite fast on the CPU.
            //
            // The main ideas come from the chess programming community that deals with 8x8 sparse matrices.
            // The LSB1 (Less Significant Bit On) is given by X & -X.
            // This same bit can be "reset" using X & (X - 1), which means that all ones remain except the LSB1.
            //
            // We can use this to iterate through all ones in order.
            // The following code scales up this principle for multiple contiguous 64bit bitboards.
            std::size_t innerIndex ; //= bitIndex(iteratedBoards[arrayIndex] & -iteratedBoards[arrayIndex]);
            if (arrayIndex < iteratedBoards.size())
            {
                do {
                    innerIndex = bitIndex(iteratedBoards[arrayIndex] & -iteratedBoards[arrayIndex]);
                    index = arrayIndex * BITBOARD_INT_SIZE + innerIndex;
                    iteratedBoards[arrayIndex] &= iteratedBoards[arrayIndex] - 1;
                    
                    arrayIndex += innerIndex / 64;
                } while (innerIndex == 64 && arrayIndex < iteratedBoards.size());
            }
            
            index = index > size ? size : index;
            
            return *this;
        }
        
        /**
         * \brief Postfix incrementation. Finds the next not-null index.
         * \returns The pre-incremented iterator.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
        
        /**
         * \brief Equality check.
         * \returns True if the indices are equal.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        friend bool operator==(const iterator &a, const iterator &b) { return a.index == b.index; }
        
        /**
         * \brief Inequality check.
         * \returns True if the indices are different.
         *
         * \author Fedyna K.
         * \version 0.2.0
         * \date 23/05/2024
         */
        friend bool operator!=(const iterator &a, const iterator &b) { return a.index != b.index; }
        
    private:
        std::size_t index;
        std::size_t arrayIndex;
        std::size_t size;
        std::vector<std::uint64_t> iteratedBoards;
    };
    
    /**
     * \struct Bitboard::reverse_iterator
     * \brief The reverse iterator over bitboards.
     *
     * Return all indexes with a bit set to 1 in reverse order.
     *
     * \note The reverse iterator is constant.
     *
     * \author Bac A.
     * \version 0.2.0
     * \date 19/09/2024
     */
    struct reverse_iterator {
        /** C++ mandatory type definitions. */
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = std::uint64_t;
        using pointer           = value_type*;
        using reference         = value_type&;
        
        /**
         * \brief Reverse iterator constructor.
         * \param[in] _index The initial index.
         * \param[in] _size The bitboard size.
         * \param[in] _boards The reference to the bitboards.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        reverse_iterator(const std::size_t _index, const std::size_t _size, const std::vector<std::uint64_t> &_boards):
        index(_index),
        size(_size),
        iteratedBoards(_boards),
        arrayIndex(index/64)
        {
            // set 0 for all bits larger than index
            if (index < size_t_maxvalue)
            {
                do {
                    ++(*this);
                } while (index >_index);
            }
                
        }
        
        /**
         * \brief Iterator dereference.
         * \returns A no-null index on the billboard or past-the-end index.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        const std::size_t operator*() const { return index; }
        
        /**
         * \brief Prefix incrementation. Finds the next not-null index.
         * \returns The reference to the current iterator.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        reverse_iterator& operator++() {
            // The reverse iterator rely on bit manipulation that are quite fast on the CPU.
            //
            // The main ideas come from the chess programming community that deals with 8x8 sparse matrices.
            // The MSB1 (Most Significant Bit On) is computed using De Bruijn Multiplication.
            // See https://www.chessprogramming.org/BitScan#Index_of_LS1B_by_Popcount
            //
            // int bitScanReverse(U64 bb) {
            // const U64 debruijn64 = C64(0x03f79d71b4cb0a89);
            // assert (bb != 0);
            // bb |= bb >> 1;
            // bb |= bb >> 2;
            // bb |= bb >> 4;
            // bb |= bb >> 8;
            // bb |= bb >> 16;
            // bb |= bb >> 32;
            // return tab64[(bb * debruijn64) >> 58];
            // }
            
            // We can use this to iterate through all ones in reverse order.
            // The following code scales up this principle for multiple contiguous 64bit bitboards.
            // MSB are stored in the rightmost 64bits bitboard.
            
            std::size_t innerIndex ; //= bitIndex(iteratedBoards[arrayIndex] & -iteratedBoards[arrayIndex]);
            if (arrayIndex != size_t_maxvalue)
            {
                do {
                    std::size_t bb(iteratedBoards[arrayIndex]) ;
                    if (bb == 0)
                        innerIndex == 64 ;
                    else
                    {
                        // Compute MSB1
                        const std::size_t debruijn64 = std::size_t(0x03f79d71b4cb0a89);
                        bb |= bb >> 1;
                        bb |= bb >> 2;
                        bb |= bb >> 4;
                        bb |= bb >> 8;
                        bb |= bb >> 16;
                        bb |= bb >> 32;
                        innerIndex = index64[(bb * debruijn64) >> 58] ;
                        // Compute corresponding index
                        index = arrayIndex * BITBOARD_INT_SIZE + innerIndex;
                        // Reset bit
                        iteratedBoards[arrayIndex] ^= 1<<innerIndex ;
                    }
                    if (innerIndex == 0)
                        --arrayIndex ;
                    
                } while (innerIndex == 64 && arrayIndex != size_t_maxvalue);
            }
            else // set iterator for end()
                index = arrayIndex ;
            
            return *this;
        }
        
        /**
         * \brief Postfix incrementation. Finds the next not-null index.
         * \returns The pre-incremented iterator.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        reverse_iterator operator++(int) { reverse_iterator tmp = *this; ++(*this); return tmp; }
        
        /**
         * \brief Equality check.
         * \returns True if the indices are equal.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        friend bool operator==(const reverse_iterator &a, const reverse_iterator &b) { return a.index == b.index; }
        
        /**
         * \brief Inequality check.
         * \returns True if the indices are different.
         *
         * \author Bac A.
         * \version 0.2.0
         * \date 23/05/2024
         */
        friend bool operator!=(const reverse_iterator &a, const reverse_iterator &b) { return a.index != b.index; }
        
    private:
        std::size_t index;
        std::size_t arrayIndex;
        std::size_t size;
        std::vector<std::uint64_t> iteratedBoards;
    };

    /**
     * \brief Bitboard default constructor. Initialize bitboard with all zeros and size 64.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard(): boards({0UL}), _size(64) {}

    /**
     * \brief Bitboard value intializer. Initialize bitboard with given vector of bitboards.
     * 
     * \param[in] _bitboards The vector of 64bit ints.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard(const std::vector<std::uint64_t> &_bitboard):
        boards(_bitboard),
        _size(_bitboard.size() * 64)
    {}
    
    /**
     * \brief Bitboard size initializer. Initialize bitboard with given size.
     * 
     * \param[in] _size The size of the bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard(const std::size_t _size, bool empty = true):
        boards(std::vector<std::uint64_t>(_size / BITBOARD_INT_SIZE + (_size % BITBOARD_INT_SIZE != 0))),
        _size(_size) 
    { if (!empty) bit_not() ; }

    /**
     * \brief Bitboard copy constructor. Initialize bitboard with given bitboard.
     * 
     * \param[in] _bitboards The bitboard to copy.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard(const Bitboard &_bitboard):
        boards(_bitboard.boards),
        _size(_bitboard._size)
    {}
    

    /**
     * \brief Bitboard assign operator. Initialize bitboard with given bitboard.
     * 
     * \param[in] _bitboard The bitboard to copy.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator=(const Bitboard &_bitboard) {
        this->boards = _bitboard.boards;
        this->_size = _bitboard._size;

        return *this;
    }

    /**
     * \brief Bitwise NOT on a bitboard.
     *
     * \returns The reference to the NOTed bitboard.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& bit_not() {
        for (std::size_t i = 0 ; i < boards.size() ; ++i)
        {
            boards[i] = ~boards[i];
        }

        return *this;
    }
    
    /**
     * \brief Bitwise NOT of a bitboard.
     *
     * \param[in] _bitboard The argument bitboard.
     * \returns A new bitboard which is the result of the NOT.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator~(const Bitboard& _bitboard) {
        Bitboard res(_bitboard);
        return res.bit_not();
    }
    
    /**
     * \brief Bitboard begin iterator that allows to loop through all non-null indices.
     * \returns The begin iterator.
     * 
     * \note The iterator is constant.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    iterator begin() const { return iterator(0, _size, boards); }

    /**
     * \brief Bitboard past-the-end iterator.
     * \returns The past-the-end iterator.
     * 
     * \note The iterator is constant.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    iterator end() const { return iterator(_size, _size, boards); }
    
    /**
     * \brief Bitboard reverse_begin reverse_iterator that allows to loop through all non-null indices in decreasing order.
     * \returns The reverse_begin iterator.
     *
     * \note The reverse_iterator is constant.
     *
     * \author Bac A.
     * \version 0.2.0
     * \date 23/05/2024
     */
    reverse_iterator reverse_begin() const { return reverse_iterator(_size-1, _size, boards); }
    reverse_iterator reverse_begin(size_t index) const { std::cout << "reverse with index " << index << std::endl ; return reverse_iterator(index, _size, boards); }

    /**
     * \brief Bitboard past-the-end reverse_iterator.
     * \returns The past-the-end reverse_iterator.
     *
     * \note The reverse_iterator is constant.
     *
     * \author Bac A.
     * \version 0.2.0
     * \date 23/05/2024
     */
    reverse_iterator reverse_end() const { return reverse_iterator(size_t_maxvalue, _size, boards); }

    /**
     * \brief Stream operator that displays all bits.
     * 
     * \param[in,out] _stream The stream to edit.
     * \param[in] _bitboard The bitboard to display.
     * \returns The edited stream.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend std::ostream& operator<<(std::ostream &_stream, const Bitboard &_bitboard) {
        // Cheat for speed, iterate through all non-null indices and fill the in-between with zeros.
        std::size_t last = -1;
        for (auto index: _bitboard) {
            _stream << std::string(index - last - 1, '0') << "1"; 
            last = index;
        }
        // Add the leading zeros.
        _stream << std::string(_bitboard._size - last - 1, '0');

        return _stream;
    }

    /**
     * \brief Bitwise OR between two bitboards.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the OR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator|(const Bitboard &_left, const Bitboard &_right) {
        Bitboard res = _left;
        res |= _right;
        return res;
    }

    /**
     * \brief Bitwise OR between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the OR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator|(const Bitboard &_left, const std::size_t &_right) {
        Bitboard res = _left;
        res |= _right;
        return res;
    }

    /**
     * \brief Bitwise OR between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the OR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator|(const std::size_t &_left, const Bitboard &_right) {
        Bitboard res = _right;
        res |= _left;
        return res;
    }

    /**
     * \brief Bitwise AND between two bitboards.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the AND.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator&(const Bitboard &_left, const Bitboard &_right) {
        Bitboard res = _left;
        res &= _right;
        return res;
    }

    /**
     * \brief Bitwise AND between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the AND.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator&(const Bitboard &_left, const std::size_t &_right) {
        Bitboard res = _left;
        res &= _right;
        return res;
    }

    /**
     * \brief Bitwise AND between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the AND.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator&(const std::size_t &_left, const Bitboard &_right) {
        Bitboard res = _right;
        res &= _left;
        return res;
    }

    /**
     * \brief Bitwise XOR between two bitboards.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the XOR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator^(const Bitboard &_left, const Bitboard &_right) {
        Bitboard res = _left;
        res ^= _right;
        return res;
    }

    /**
     * \brief Bitwise XOR between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the XOR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator^(const Bitboard &_left, const std::size_t &_right) {
        Bitboard res = _left;
        res ^= _right;
        return res;
    }

    /**
     * \brief Bitwise XOR between a bitboard and a single bit at given position.
     * 
     * \param[in] _left The left hand side.
     * \param[in] _right The right hand side.
     * \returns A new bitboard which is the result of the XOR.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    friend Bitboard operator^(const std::size_t &_left, const Bitboard &_right) {
        Bitboard res = _right;
        res ^= _left;
        return res;
    }

    /**
     * \brief Bitwise OR between two bitboards.
     * 
     * \param[in] _other The other bitboard.
     * \returns The reference to the ORed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator|=(const Bitboard &_other) {
        for (std::size_t i = 0 ; i < boards.size() ; i++) {
            boards[i] |= _other.boards[i];
        }

        return *this;
    }
    
    /**
     * \brief Bitwise OR between a bitboard and a single bit at given position.
     * 
     * \param[in] _other The bit position.
     * \returns The reference to the ORed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator|=(const std::size_t &_other) {
        std::size_t boardIndex = _other / BITBOARD_INT_SIZE;
        std::uint64_t bit = 1ULL << (_other % BITBOARD_INT_SIZE);

        boards[boardIndex] |= bit;
        return *this;
    }

    /**
     * \brief Bitwise AND between two bitboards.
     * 
     * \param[in] _other The other bitboard.
     * \returns The reference to the ANDed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator&=(const Bitboard &_other) {
        for (std::size_t i = 0 ; i < boards.size() ; i++) {
            boards[i] &= _other.boards[i];
        }

        return *this;
    }

    /**
     * \brief Bitwise AND between a bitboard and a single bit at given position.
     * 
     * \param[in] _other The bit position.
     * \returns The reference to the ANDed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator&=(const std::size_t &_other) {
        std::size_t boardIndex = _other / BITBOARD_INT_SIZE;
        std::uint64_t bit = 1ULL << (_other % BITBOARD_INT_SIZE);

        boards[boardIndex] &= bit;
        return *this;
    }

    /**
     * \brief Bitwise XOR between two bitboards.
     * 
     * \param[in] _other The other bitboard.
     * \returns The reference to the XORed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator^=(const Bitboard &_other) {
        for (std::size_t i = 0 ; i < boards.size() ; i++) {
            boards[i] ^= _other.boards[i];
        }

        return *this;
    }

    /**
     * \brief Bitwise XOR between a bitboard and a single bit at given position.
     * 
     * \param[in] _other The bit position.
     * \returns The reference to the XORed bitboard.
     * 
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    Bitboard& operator^=(const std::size_t &_other) {
        std::size_t boardIndex = _other / BITBOARD_INT_SIZE;
        std::uint64_t bit = 1ULL << (_other % BITBOARD_INT_SIZE);

        boards[boardIndex] ^= bit;
        return *this;
    }

    /**
     * \brief Toggle on and off a given bit.
     * 
     * \param[in] _index The bit position.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    void toggle(const std::size_t &_index) {
        *this ^= _index;
    }

    /**
     * \brief Toggle on a given bit.
     * 
     * \param[in] _index The bit position.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    void setOn(const std::size_t &_index) {
        *this |= _index;
    }

    /**
     * \brief Toggle off a given bit.
     * 
     * \param[in] _index The bit position.
     *
     * \author Fedyna K.
     * \version 0.2.0
     * \date 23/05/2024
     */
    void setOff(const std::size_t &_index) {
        std::size_t boardIndex = _index / BITBOARD_INT_SIZE;
        std::uint64_t bit = 1ULL << (_index % BITBOARD_INT_SIZE);

        boards[boardIndex] &= ~bit;
    }

    bool isOn(const std::size_t &_index) const {
        std::size_t boardIndex = _index / BITBOARD_INT_SIZE;
        std::uint64_t bit = 1ULL << (_index % BITBOARD_INT_SIZE);

        return (boards[boardIndex] & bit) != 0;
    }
    
    size_t size() const { return _size ; }
};

}

#endif
