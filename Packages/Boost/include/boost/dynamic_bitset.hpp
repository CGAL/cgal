// (C) Copyright Chuck Allison and Jeremy Siek 2001, 2002.
//
// Permission to copy, use, modify, sell and distribute this software
// is granted provided this copyright notice appears in all
// copies. This software is provided "as is" without express or
// implied warranty, and with no claim as to its suitability for any
// purpose.

// With optimizations, bug fixes, and improvements by Gennaro Prota.

// See http://www.boost.org/libs/dynamic_bitset for documentation.

// -------------------------------------
// CHANGE LOG:
//
// - corrected workaround for Dinkum lib's allocate() [GP]
// - changed macro test for old iostreams [GP]
// - removed #include <vector> for now. [JGS]
// - Added __GNUC__ to compilers that cannot handle the constructor from basic_string. [JGS]
// - corrected to_block_range [GP]
// - corrected from_block_range [GP]
// - Removed __GNUC__ from compilers that cannot handle the constructor
//     from basic_string and added the workaround suggested by GP. [JGS]
// - Removed __BORLANDC__ from the #if around the basic_string
//     constructor. Luckily the fix by GP for g++ also fixes Borland. [JGS]

#ifndef BOOST_DYNAMIC_BITSET_HPP
#define BOOST_DYNAMIC_BITSET_HPP

#include <boost/config.hpp>
#include <cassert>
#include <string>
#include <iosfwd>
#include <cstring>             // for memset, memcpy, memcmp, etc.
#include <stdexcept>           // for std::overflow_error
#include <algorithm>           // for std::swap, std::min, std::copy, std::fill

#if defined(__GNUC__) && !defined(__SGI_STL_PORT)
#include <ctype.h>
#else
#include <cctype>               // for isspace
#endif

#include "boost/dynamic_bitset_fwd.hpp" //G.P.S.
#include "boost/detail/dynamic_bitset.hpp"

#if defined (__STL_CONFIG_H) && !defined (__STL_USE_NEW_IOSTREAMS)
#define BOOST_OLD_IOSTREAMS
#endif


#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
//  in certain situations VC++ requires a redefinition of
//  default template arguments, in contrast with 14.1/12
//
# define BOOST_WORKAROUND_REPEAT_DEFAULT_TEMPLATE_ARGUMENTS // macro 'local' to this file
#endif



namespace boost {

#ifdef BOOST_WORKAROUND_REPEAT_DEFAULT_TEMPLATE_ARGUMENTS
 template <typename Block = unsigned long, typename Allocator = std::allocator<Block> >
#else
 template <typename Block, typename Allocator>
#endif

class dynamic_bitset :
#ifdef BOOST_DYN_BITSET_USE_FRIENDS
    private
#else
    public
#endif
    detail::dynamic_bitset_base<Block, Allocator>
{
  // Portability note: member function templates are defined inside
  // this class definition to avoid problems with VC++. Similarly,
  // with the member functions of the nested class.
public:
    typedef Block block_type;
    typedef std::size_t size_type;
    enum { bits_per_block = CHAR_BIT * sizeof(Block) };

    // reference to a bit
    class reference
    {
        friend class dynamic_bitset<Block, Allocator>;
        dynamic_bitset* bs;
        size_type bit;
        reference(); // intentionally not implemented
        reference(dynamic_bitset& bs_, size_type bit_) : bs(&bs_), bit(bit_){ }
    public:
        reference& operator=(bool value)          // for b[i] = x
        {
            if (value)
                bs->set(bit);
            else
                bs->reset(bit);
            return *this;
        }
        reference& operator|=(bool value)         // for b[i] |= x
        {
            if (value)
                bs->set(bit);
            return *this;
        }
        reference& operator&=(bool value)         // for b[i] &= x
        {
            if (! (value && bs->test(bit)))
                bs->reset(bit);
            return *this;
        }
        reference& operator^=(bool value)         // for b[i] ^= x
        {
            bs->set(bit, bs->test(bit) ^ value);
            return *this;
        }
        reference& operator-=(bool value)         // for b[i] -= x
        {
            if (!value)
                bs->reset(bit);
            return *this;
        }
        reference& operator=(const reference& j)  // for b[i] = b[j]
        {
            if (j.bs->test(j.bit))
                bs->set(bit);
            else
                bs->reset(bit);
            return *this;
        }
        reference& operator|=(const reference& j) // for b[i] |= b[j]
        {
            if (j.bs->test(j.bit))
                bs->set(bit);
            return *this;
        }
        reference& operator&=(const reference& j) // for b[i] &= b[j]
        {
            if (! (j.bs->test(j.bit) && bs->test(bit)))
                bs->reset(bit);
            return *this;
        }
        reference& operator^=(const reference& j) // for b[i] ^= b[j]
        {
            bs->set(bit, bs->test(bit) ^ j.bs->test(j.bit));
            return *this;
        }
        reference& operator-=(const reference& j) // for b[i] -= b[j]
        {
            if (!j.bs->test(j.bit))
                bs->reset(bit);
            return *this;
        }
        bool operator~() const                    // flips the bit
        {
            return ! bs->test(bit);
        }
        operator bool() const                     // for x = b[i]
        {
            return bs->test(bit);
        }
        reference& flip()                         // for b[i].flip();
        {
            bs->flip(bit);
            return *this;
        }
    };
    typedef bool const_reference;

    // constructors, etc.
    explicit
    dynamic_bitset(const Allocator& alloc = Allocator());

    explicit
    dynamic_bitset(size_type num_bits, unsigned long value = 0,
               const Allocator& alloc = Allocator());

    // from string
#if defined(BOOST_OLD_IOSTREAMS)
    explicit
    dynamic_bitset(const std::string& s,
               std::string::size_type pos = 0,
               std::string::size_type n = std::string::npos,
               const Allocator& alloc = Allocator())
        : detail::dynamic_bitset_base<Block, Allocator>
            (std::min(n, s.size() - pos), alloc)
#else
    // The parenthesis around std::basic_string<CharT, Traits, Alloc>::npos
    // in the code below are to avoid a g++ 3.2 bug and a Borland bug. -JGS
    template <typename CharT, typename Traits, typename Alloc>
    explicit
    dynamic_bitset(const std::basic_string<CharT, Traits, Alloc>& s,
        typename std::basic_string<CharT, Traits, Alloc>::size_type pos = 0,
        typename std::basic_string<CharT, Traits, Alloc>::size_type n
            = (std::basic_string<CharT, Traits, Alloc>::npos),
        const Allocator& alloc = Allocator())
        : detail::dynamic_bitset_base<Block, Allocator>
            (std::min(n, s.size() - pos), alloc)
#endif
    {
        // Locate sub string
        assert(pos <= s.length());
        from_string(s, pos, std::min(n, s.size() - pos));
    }

    // The first bit in *first is the least significant bit, and the
    // last bit in the block just before *last is the most significant bit.
    template <typename BlockInputIterator>
    dynamic_bitset(BlockInputIterator first, BlockInputIterator last,
               const Allocator& alloc = Allocator())
        : detail::dynamic_bitset_base<Block, Allocator>
            (detail::initial_num_blocks(first, last)
            * bits_per_block, alloc)
    {
        if (first != last) {
            if (this->m_num_bits == 0) { // dealing with input iterators
                this->append(first, last);
            } else {
                // dealing with forward iterators, memory has been allocated
                for (std::size_t i = 0; first != last; ++first, ++i)
                    set_block_(i, *first);
            }
        }
    }


    // copy constructor
    dynamic_bitset(const dynamic_bitset& b);

    void swap(dynamic_bitset& b);

    dynamic_bitset& operator=(const dynamic_bitset& b);

    // size changing operations
    void resize(size_type num_bits, bool value = false);
    void clear();
    void push_back(bool bit);
    void append(Block block);

    // This is declared inside the class to avoid compiler bugs.
    template <typename BlockInputIterator>
    void append(BlockInputIterator first, BlockInputIterator last)
    {
        if (first != last) {
            std::size_t nblocks = detail::initial_num_blocks(first, last);
            if (nblocks == 0) { // dealing with input iterators
                for (; first != last; ++first)
                    append(*first);
            } else { // dealing with forward iterators
                if (size() % bits_per_block == 0) {
                    std::size_t old_nblocks = this->m_num_blocks;
                    resize(size() + nblocks * bits_per_block);
                    for (std::size_t i = old_nblocks; first != last; ++first)
                        set_block_(i++, *first);
                } else {
                    // probably should optimize this,
                    // but I'm sick of bit twiddling
                    for (; first != last; ++first)
                        append(*first);
                }
            }
        }
    }


    // bitset operations
    dynamic_bitset& operator&=(const dynamic_bitset& b);
    dynamic_bitset& operator|=(const dynamic_bitset& b);
    dynamic_bitset& operator^=(const dynamic_bitset& b);
    dynamic_bitset& operator-=(const dynamic_bitset& b);
    dynamic_bitset& operator<<=(size_type n);
    dynamic_bitset& operator>>=(size_type n);
    dynamic_bitset operator<<(size_type n) const;
    dynamic_bitset operator>>(size_type n) const;

    // basic bit operations
    dynamic_bitset& set(size_type n, bool val = true);
    dynamic_bitset& set();
    dynamic_bitset& reset(size_type n);
    dynamic_bitset& reset();
    dynamic_bitset& flip(size_type n);
    dynamic_bitset& flip();
    bool test(size_type n) const;
    bool any() const;
    bool none() const;
    dynamic_bitset operator~() const;
    size_type count() const;

    // subscript
    reference operator[](size_type pos) { return reference(*this, pos); }
    bool operator[](size_type pos) const { return test_(pos); } //[gps]

    unsigned long to_ulong() const;

    size_type size() const;
    size_type num_blocks() const;

    bool is_subset_of(const dynamic_bitset& a) const;
    bool is_proper_subset_of(const dynamic_bitset& a) const;


#ifdef BOOST_DYN_BITSET_USE_FRIENDS
    // lexicographical comparison
    template <typename B, typename A>
    friend bool operator==(const dynamic_bitset<B, A>& a,
                           const dynamic_bitset<B, A>& b);
    template <typename B, typename A>
    friend bool operator<(const dynamic_bitset<B, A>& a,
                          const dynamic_bitset<B, A>& b);
    template <typename B, typename A>
    friend bool operator>(const dynamic_bitset<B, A>& a,
                          const dynamic_bitset<B, A>& b);

    template <typename B, typename A, typename BlockOutputIterator>
    friend void to_block_range(const dynamic_bitset<B, A>& b,
                               BlockOutputIterator result);

    template <typename BlockIterator, typename B, typename A>
    friend void from_block_range(BlockIterator first, BlockIterator last,
                                 dynamic_bitset<B, A>& result);

    template <typename B, typename A, typename CharT, typename Alloc>
    friend void dump_to_string(const dynamic_bitset<B, A>& b,
                               std::basic_string<CharT, Alloc>& s);
#endif

private:
    void m_zero_unused_bits();
    void set_(size_type bit);
    bool set_(size_type bit, bool val);
    void reset_(size_type bit);
    bool test_(size_type bit) const;
    void set_block_(size_type blocknum, Block b);

public:

    // This is templated on the whole String instead of just CharT,
    // Traits, Alloc to avoid compiler bugs.
    template <typename String>
    void from_string(const String& s, typename String::size_type pos,
                     typename String::size_type rlen)
    {
        reset(); // bugfix [gps]
        size_type const tot = std::min (rlen, s.length()); // bugfix [gps]

        // Assumes string contains only 0's and 1's
        for (size_type i = 0; i < tot; ++i) {
        if (s[pos + tot - i - 1] == '1') {
                set_(i);
            } else {
            assert(s[pos + tot - i - 1] == '0');
            }
        }
    }

};

// Global Functions:

// comparison
template <typename Block, typename Allocator>
bool operator!=(const dynamic_bitset<Block, Allocator>& a,
                const dynamic_bitset<Block, Allocator>& b);

template <typename Block, typename Allocator>
bool operator<=(const dynamic_bitset<Block, Allocator>& a,
                const dynamic_bitset<Block, Allocator>& b);

template <typename Block, typename Allocator>
bool operator>(const dynamic_bitset<Block, Allocator>& a,
               const dynamic_bitset<Block, Allocator>& b);

template <typename Block, typename Allocator>
bool operator>=(const dynamic_bitset<Block, Allocator>& a,
                const dynamic_bitset<Block, Allocator>& b);

// stream operators
#ifdef BOOST_OLD_IOSTREAMS
template <typename Block, typename Allocator>
std::ostream& operator<<(std::ostream& os,
                         const dynamic_bitset<Block, Allocator>& b);

template <typename Block, typename Allocator>
std::istream& operator>>(std::istream& is, dynamic_bitset<Block,Allocator>& b);
#else
template <typename CharT, typename Traits, typename Block, typename Allocator>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& os,
           const dynamic_bitset<Block, Allocator>& b);

template <typename CharT, typename Traits, typename Block, typename Allocator>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& is,
           dynamic_bitset<Block, Allocator>& b);
#endif

// bitset operations
template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator&(const dynamic_bitset<Block, Allocator>& b1,
          const dynamic_bitset<Block, Allocator>& b2);

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator|(const dynamic_bitset<Block, Allocator>& b1,
          const dynamic_bitset<Block, Allocator>& b2);

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator^(const dynamic_bitset<Block, Allocator>& b1,
          const dynamic_bitset<Block, Allocator>& b2);

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator-(const dynamic_bitset<Block, Allocator>& b1,
          const dynamic_bitset<Block, Allocator>& b2);


template <typename Block, typename Allocator, typename CharT, typename Alloc>
void
to_string(const dynamic_bitset<Block, Allocator>& b,
          std::basic_string<CharT, Alloc>& s);

template <typename Block, typename Allocator, typename BlockOutputIterator>
void
to_block_range(const dynamic_bitset<Block, Allocator>& b,
               BlockOutputIterator result);

template <typename BlockIterator, typename B, typename A>
inline void
from_block_range(BlockIterator first, BlockIterator last,
                 dynamic_bitset<B, A>& result);

//=============================================================================
// dynamic_bitset implementation

#ifdef BOOST_OLD_IOSTREAMS
template <typename Block, typename Allocator>
inline std::ostream&
operator<<(std::ostream& os,
           const typename dynamic_bitset<Block, Allocator>::reference& br)
{
    return os << (bool)br;
}
#else
template <typename CharT, typename Traits, typename Block, typename Allocator>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& os,
           const typename dynamic_bitset<Block, Allocator>::reference& br)
{
    return os << (bool)br;
}
#endif

//-----------------------------------------------------------------------------
// constructors, etc.

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>::dynamic_bitset(const Allocator& alloc)
  : detail::dynamic_bitset_base<Block, Allocator>(0, alloc) { }

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>::
dynamic_bitset(size_type num_bits, unsigned long value, const Allocator& alloc)
  : detail::dynamic_bitset_base<Block, Allocator>(num_bits, alloc)
{
  const size_type M = std::min(sizeof(unsigned long) * CHAR_BIT, num_bits);
  for(size_type i = 0; i < M; ++i, value >>= 1) // [G.P.S.] to be optimized
    if ( value & 0x1 )
      set_(i);
}

// copy constructor
template <typename Block, typename Allocator>
inline dynamic_bitset<Block, Allocator>::
dynamic_bitset(const dynamic_bitset& b)
  : detail::dynamic_bitset_base<Block, Allocator>(b.size(), b.m_alloc)
{
    using namespace std;
    memcpy(this->m_bits, b.m_bits, this->m_num_blocks * sizeof(Block));
}

template <typename Block, typename Allocator>
inline void dynamic_bitset<Block, Allocator>::
swap(dynamic_bitset<Block, Allocator>& b)
{
    std::swap(this->m_bits, b.m_bits);
    std::swap(this->m_num_bits, b.m_num_bits);
    std::swap(this->m_num_blocks, b.m_num_blocks);
    std::swap(this->m_alloc, b.m_alloc);
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>& dynamic_bitset<Block, Allocator>::
operator=(const dynamic_bitset<Block, Allocator>& b)
{
    dynamic_bitset<Block, Allocator> tmp(b);
    this->swap(tmp);
    return *this;
}

//-----------------------------------------------------------------------------
// size changing operations

template <typename Block, typename Allocator>
void dynamic_bitset<Block, Allocator>::
resize(size_type num_bits, bool value)
{
  if (num_bits == size())
    return;
  size_type new_nblocks = this->calc_num_blocks(num_bits);
  Block* d = this->m_alloc.allocate(new_nblocks, static_cast<void const *>(0));
  if (num_bits < size()) { // shrink
    std::copy(this->m_bits, this->m_bits + new_nblocks, d);
    std::swap(d, this->m_bits);
    this->m_alloc.deallocate(d, this->m_num_blocks);
  } else { // grow
    std::copy(this->m_bits, this->m_bits + this->m_num_blocks, d);
    Block val = value? ~static_cast<Block>(0) : static_cast<Block>(0);
    std::fill(d + this->m_num_blocks, d + new_nblocks, val);
    std::swap(d, this->m_bits);
    for (std::size_t i = this->m_num_bits;
         i < this->m_num_blocks * bits_per_block; ++i)
      set_(i, value);
    if (d != 0)
      this->m_alloc.deallocate(d, this->m_num_blocks);
  }
  this->m_num_bits = num_bits;
  this->m_num_blocks = this->calc_num_blocks(num_bits);
  m_zero_unused_bits();
}

template <typename Block, typename Allocator>
void dynamic_bitset<Block, Allocator>::
clear()
{
  if (this->m_bits != 0) {
    this->m_alloc.deallocate(this->m_bits, this->m_num_blocks);
    this->m_bits = 0;
    this->m_num_bits = 0;
    this->m_num_blocks = 0;
  }
}


template <typename Block, typename Allocator>
void dynamic_bitset<Block, Allocator>::
push_back(bool bit)
{
  this->resize(this->size() + 1);
  set_(this->size() - 1, bit);
}

template <typename Block, typename Allocator>
void dynamic_bitset<Block, Allocator>::
append(Block value)
{
  std::size_t old_size = size();
  resize(old_size + bits_per_block);
  if (size() % bits_per_block == 0)
    set_block_(this->m_num_blocks - 1, value);
  else {
      // G.P.S. to be optimized
    for (std::size_t i = old_size; i < size(); ++i, value >>= 1)
      set_(i, value & 1);
  }
}


//-----------------------------------------------------------------------------
// bitset operations
template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator&=(const dynamic_bitset& rhs)
{
    assert(size() == rhs.size());
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        this->m_bits[i] &= rhs.m_bits[i];
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator|=(const dynamic_bitset& rhs)
{
    assert(size() == rhs.size());
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        this->m_bits[i] |= rhs.m_bits[i];
    m_zero_unused_bits();
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator^=(const dynamic_bitset& rhs)
{
    assert(size() == rhs.size());
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        this->m_bits[i] ^= rhs.m_bits[i];
    m_zero_unused_bits();
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator-=(const dynamic_bitset& rhs)
{
    assert(size() == rhs.size());
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        this->m_bits[i] = this->m_bits[i] & ~rhs.m_bits[i];
    m_zero_unused_bits();
    return *this;
}

/* [gps] - snipped
template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator<<=(size_type n)
{
    if (n >= this->m_num_bits)
        reset();
    else
    {
        size_type i;
        for (i = this->m_num_bits - 1; i > n; --i)
            set_(i,test_(i-n));
        if (i == n) // careful, unsigned can't go negative!
            set_(i,test_(i-n));
        for (i = 0; i < n; ++i)
            reset_(i);
    }
    return *this;
}*/


template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator<<=(size_type n)
{
    if (n >= this->m_num_bits)
        return reset();
    //else
    if (n > 0)
    {
        size_type  const last  = this->m_num_blocks - 1; // m_num_blocks is >= 1
        size_type  const div   = n / bits_per_block; // div is <= last
        size_type  const r     = n % bits_per_block;

        // PRE: div != 0  or  r != 0

        if (r != 0) {

            block_type const rs = bits_per_block - r;

            for (size_type i = last-div; i>0; --i) {
                this->m_bits[i+div] = (this->m_bits[i] << r) | (this->m_bits[i-1] >> rs);
            }
            this->m_bits[div] = this->m_bits[0] << r;

        }
        else {
            for (size_type i = last-div; i>0; --i) {
                this->m_bits[i+div] = this->m_bits[i];
            }
            this->m_bits[div] = this->m_bits[0];
        }


        // div blocks are zero filled at the less significant end
        std::fill(this->m_bits, this->m_bits+div, static_cast<block_type>(0));


    }

    return *this;


}




/* [gps] - snipped
template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::operator>>=(size_type n)
{
    if (n >= this->m_num_bits)
        reset();
    else
    {
        size_type i;
        for (i = 0; i < this->m_num_bits - n; ++i)
            set_(i,test_(i+n));

        for (i = this->m_num_bits - n; i < this->m_num_bits; ++i)
            reset_(i);
    }
    return *this;
}*/




// NOTE: this assumes that within a single block bits are
//       numbered from right to left. G.P.S.
//
//      static Block offset(size_type bit)
//        { return  bit % bits_per_block; }
//
//
// In the implementation below the 'if (r != 0)' is logically
// unnecessary. It's there as an optimization only: in fact
// for r==0 the first branch becomes the second one with the
// b[last-div] = b[last] >> r; statement that does the work of
// the last iteration.
//
template <typename B, typename A>
dynamic_bitset<B, A> & dynamic_bitset<B, A>::operator>>=(size_type n) {
    if (n >= this->m_num_bits) {
        return reset();
    }
    //else
    if (n>0){

        size_type  const last  = this->m_num_blocks - 1; // m_num_blocks is >= 1
        size_type  const div   = n / bits_per_block; // div is <= last
        size_type  const r     = n % bits_per_block;

        // PRE: div != 0  or  r != 0

        if (r != 0) {

            block_type const ls = bits_per_block - r;

            for (size_type i = div; i < last; ++i) {
                this->m_bits[i-div] = (this->m_bits[i] >> r) | (this->m_bits[i+1]  << ls);
            }
            // r bits go to zero
            this->m_bits[last-div] = this->m_bits[last] >> r;
        }

        else {
            for (size_type i = div; i <= last; ++i) {
                this->m_bits[i-div] = this->m_bits[i];
            }
            // note the '<=': the last iteration 'absorbs'
            // this->m_bits[last-div] = this->m_bits[last] >> 0;
        }



        // div blocks are zero filled at the most significant end
        std::fill(this->m_bits+(this->m_num_blocks-div), this->m_bits+this->m_num_blocks, static_cast<block_type>(0));
    }

    return *this;
}







template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
dynamic_bitset<Block, Allocator>::operator<<(size_type n) const
{
    dynamic_bitset r(*this);
    return r <<= n;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
dynamic_bitset<Block, Allocator>::operator>>(size_type n) const
{
    dynamic_bitset r(*this);
    return r >>= n;
}


//-----------------------------------------------------------------------------
// basic bit operations

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::set(size_type pos, bool val)
{
    assert(pos < this->m_num_bits);
    set_(pos, val);
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::set()
{
  if (this->m_num_bits > 0) {
    using namespace std;
    memset(this->m_bits, ~0u, this->m_num_blocks * sizeof(this->m_bits[0]));
    m_zero_unused_bits();
  }
  return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::reset(size_type pos)
{
    assert(pos < this->m_num_bits);
    reset_(pos);
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::reset()
{
  if (this->m_num_bits > 0) {
    using namespace std;
    memset(this->m_bits, 0, this->m_num_blocks * sizeof(this->m_bits[0]));
  }
  return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::flip(size_type pos)
{
    assert(pos < this->m_num_bits);
    this->m_bits[this->word(pos)] ^= this->mask1(pos);
    return *this;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>&
dynamic_bitset<Block, Allocator>::flip()
{
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        this->m_bits[i] = ~this->m_bits[i];
    m_zero_unused_bits();
    return *this;
}

template <typename Block, typename Allocator>
bool dynamic_bitset<Block, Allocator>::test(size_type pos) const
{
    assert(pos < this->m_num_bits);
    return test_(pos);
}

template <typename Block, typename Allocator>
bool dynamic_bitset<Block, Allocator>::any() const
{
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        if (this->m_bits[i])
            return 1;
    return 0;
}

template <typename Block, typename Allocator>
inline bool dynamic_bitset<Block, Allocator>::none() const
{
    return !any();
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
dynamic_bitset<Block, Allocator>::operator~() const
{
    dynamic_bitset b(*this);
    b.flip();
    return b;
}


/* snipped: [gps]

The following is the straightforward implementation of count(), which
we leave here in a comment for documentation purposes.

template <typename Block, typename Allocator>
typename dynamic_bitset<Block, Allocator>::size_type
dynamic_bitset<Block, Allocator>::count() const
{
    size_type sum = 0;
    for (size_type i = 0; i != this->m_num_bits; ++i)
        if (test_(i))
            ++sum;
    return sum;
}

The actual algorithm used is based on using a lookup
table.


  The basic idea of the method is to pick up X bits at a time
  from the internal array of blocks and consider those bits as
  the binary representation of a number N. Then, to use a table
  of 1<<X elements where table[N] is the number of '1' digits
  in the binary representation of N (i.e. in our X bits).

  Note that the table can be oversized (i.e. can even have more
  than 1<<X elements; in that case only the first 1<<X will be
  actually used) but it cannot be undersized.
  In this implementation X is 8 (but can be easily changed: you
  just have to change the definition of count<>::max_bits) and
  the internal array of blocks is seen as an array of bytes: if
  a byte has exactly 8 bits then it's enough to sum the value
  of table[B] for each byte B. Otherwise 8 bits at a time are
  'extracted' from each byte by using another loop. As a further
  efficiency consideration note that even if you have, let's say,
  32-bit chars the inner loop will not do 4 (i.e. 32/8) iterations,
  unless you have at least one bit set in the highest 8 bits of the
  byte.

  Note also that the outmost if/else is not necessary but is there
  to help the optimizer (and one of the two branches is always dead
  code).

*/


template <typename Block, typename Allocator>
typename dynamic_bitset<Block, Allocator>::size_type
dynamic_bitset<Block, Allocator>::count() const
{
    using detail::byte_t;

    const byte_t * p = reinterpret_cast<const byte_t*>(this->m_bits);
    const byte_t * past_end = p + this->m_num_blocks * sizeof(Block);

    size_type num = 0;
    unsigned int const max_bit = detail::count<>::max_bit;

    if (CHAR_BIT <= max_bit) { // table is large enough
        while (p < past_end) {
            num += detail::count<>::table[*p];
            ++p;
        }
    }
    else {
        while (p < past_end) {
            // this inner loop 'extracts' max_bit bits at a time
            // from the byte at address p, and thus allows to use
            // our (small) table for any (high) CHAR_BIT value
            byte_t value = *p;
            do {
                num += detail::count<>::table[value & ((1<<max_bit)-1)];
            } while (value >>= max_bit);

            ++p;
        }
    }


    return num;
}


//-----------------------------------------------------------------------------
// conversions

// take as ref param instead?
template <typename Block, typename Allocator, typename CharT, typename Alloc>
void
to_string(const dynamic_bitset<Block, Allocator>& b,
          std::basic_string<CharT, Alloc>& s)
{
    s.assign(b.size(), '0');
    for (std::size_t i = 0; i < b.size(); ++i)
        if (b.test(i)) // [G.P.S.]
            s[b.size() - 1 - i] = '1';
}


// Differently from to_string this function dumps out
// every bit of the internal representation (useful
// for debugging purposes)
//
template <typename B, typename A, typename CharT, typename Alloc>
void
dump_to_string(const dynamic_bitset<B, A>& b,
               std::basic_string<CharT, Alloc>& s)
{
    std::size_t const len = b.m_num_blocks * (dynamic_bitset<B, A>::bits_per_block);
    s.assign(len, '0');
    for (std::size_t i = 0; i != len; ++i)
      if (b[i])// could use test_ here, but we have friend issues.-JGS
            s[len - 1 - i] = '1';
}

template <typename Block, typename Allocator, typename BlockOutputIterator>
void
to_block_range(const dynamic_bitset<Block, Allocator>& b,
               BlockOutputIterator result)
{
    assert(b.size() != 0 || b.num_blocks() == 0);
    std::copy (b.m_bits, b.m_bits + b.m_num_blocks, result);
}

template <typename BlockIterator, typename B, typename A>
inline void
from_block_range(BlockIterator first, BlockIterator last,
                 dynamic_bitset<B, A>& result)
{
    // PRE: distance(first, last) == numblocks()
    std::copy (first, last, result.m_bits);
}

template <typename Block, typename Allocator>
unsigned long dynamic_bitset<Block, Allocator>::
to_ulong() const
{
  const std::overflow_error
    overflow("boost::bit_set::operator unsigned long()");

  if (this->m_num_blocks == 0)
    return 0;

  if (sizeof(Block) >= sizeof(unsigned long)) {
    for (size_type i = 1; i < this->m_num_blocks; ++i)
      if (this->m_bits[i])
        throw overflow;
    const Block mask = static_cast<Block>(static_cast<unsigned long>(-1));
    if (this->m_bits[0] & ~mask)
      throw overflow;
    size_type N = std::min(sizeof(unsigned long) * CHAR_BIT, this->size());
    unsigned long num = 0;
    for (size_type j = 0; j < N; ++j)
      if (this->test(j))
        num |= (1 << j);
    return num;
  }
  else { // sizeof(Block) < sizeof(unsigned long).
    const size_type nwords =
      (sizeof(unsigned long) + sizeof(Block) - 1) / sizeof(Block);

    if (this->m_num_blocks > nwords) {
      for (size_type i = nwords; i < this->m_num_blocks; ++i)
        if (this->m_bits[i])
          throw overflow;
    }

    unsigned long result = 0;
    size_type N = std::min(sizeof(unsigned long) * CHAR_BIT, this->size());
    for (size_type i = 0; i < N; ++i)
      if (this->test(i))
        result |= (1 << i);
    return result;
  }
}


template <typename Block, typename Allocator>
inline typename dynamic_bitset<Block, Allocator>::size_type
dynamic_bitset<Block, Allocator>::size() const
{
    return this->m_num_bits;
}

template <typename Block, typename Allocator>
inline typename dynamic_bitset<Block, Allocator>::size_type
dynamic_bitset<Block, Allocator>::num_blocks() const
{
    return this->m_num_blocks;
}


template <typename Block, typename Allocator>
bool dynamic_bitset<Block, Allocator>::
is_subset_of(const dynamic_bitset<Block, Allocator>& a) const
{
    assert(this->size() == a.size());
    for (size_type i = 0; i < this->m_num_blocks; ++i)
        if (this->m_bits[i] & ~a.m_bits[i])
            return false;
    return true;
}

template <typename Block, typename Allocator>
bool dynamic_bitset<Block, Allocator>::
is_proper_subset_of(const dynamic_bitset<Block, Allocator>& a) const
{
    assert(this->size() == a.size());
    bool proper = false;
    for (size_type i = 0; i < this->m_num_blocks; ++i) {
        Block bt = this->m_bits[i], ba = a.m_bits[i];
        if (ba & ~bt)
            proper = true;
        if (bt & ~ba)
            return false;
    }
    return proper;
}

//-----------------------------------------------------------------------------
// comparison

template <typename Block, typename Allocator>
bool operator==(const dynamic_bitset<Block, Allocator>& a,
                const dynamic_bitset<Block, Allocator>& b)
{
    using namespace std;
    return (a.m_num_bits == b.m_num_bits) &&
      ((a.m_num_bits == 0) ||
        !memcmp(a.m_bits, b.m_bits, a.m_num_blocks * sizeof(a.m_bits[0])));
}

template <typename Block, typename Allocator>
inline bool operator!=(const dynamic_bitset<Block, Allocator>& a,
                       const dynamic_bitset<Block, Allocator>& b)
{
    return !(a == b);
}

template <typename Block, typename Allocator>
bool operator<(const dynamic_bitset<Block, Allocator>& a,
               const dynamic_bitset<Block, Allocator>& b)
{
    assert(a.size() == b.size());
    typedef typename dynamic_bitset<Block, Allocator>::size_type size_type;

    if (a.size() == 0)
      return false;

    // Since we are storing the most significant bit
    // at pos == size() - 1, we need to do the memcmp in reverse.

    // Compare a block at a time
    for (size_type i = a.m_num_blocks - 1; i > 0; --i)
      if (a.m_bits[i] < b.m_bits[i])
        return true;
      else if (a.m_bits[i] > b.m_bits[i])
        return false;

    if (a.m_bits[0] < b.m_bits[0])
      return true;
    else
      return false;
}

template <typename Block, typename Allocator>
inline bool operator<=(const dynamic_bitset<Block, Allocator>& a,
                       const dynamic_bitset<Block, Allocator>& b)
{
    return !(a > b);
}

template <typename Block, typename Allocator>
inline bool operator>(const dynamic_bitset<Block, Allocator>& a,
                      const dynamic_bitset<Block, Allocator>& b)
{
    assert(a.size() == b.size());
    typedef typename dynamic_bitset<Block, Allocator>::size_type size_type;

    if (a.size() == 0)
      return false;

    // Since we are storing the most significant bit
    // at pos == size() - 1, we need to do the memcmp in reverse.

    // Compare a block at a time
    for (size_type i = a.m_num_blocks - 1; i > 0; --i)
      if (a.m_bits[i] < b.m_bits[i])
        return false;
      else if (a.m_bits[i] > b.m_bits[i])
        return true;

    if (a.m_bits[0] > b.m_bits[0])
      return true;
    else
      return false;
}

template <typename Block, typename Allocator>
inline bool operator>=(const dynamic_bitset<Block, Allocator>& a,
                       const dynamic_bitset<Block, Allocator>& b)
{
    return !(a < b);
}

//-----------------------------------------------------------------------------
// stream operations

#ifdef BOOST_OLD_IOSTREAMS
template < typename Block, typename Allocator>
std::ostream&
operator<<(std::ostream& os, const dynamic_bitset<Block, Allocator>& b)
{
    std::string s;
    to_string(b, s);
    os << s.c_str();
    return os;
}
#else
template <typename CharT, typename Traits, typename Block, typename Allocator>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& os,
           const dynamic_bitset<Block, Allocator>& b)
{
    std::basic_string<CharT, Traits> s;
    to_string(b, s);
    os << s;
    return os;
}
#endif

#ifdef BOOST_OLD_IOSTREAMS
template <typename Block, typename Allocator>
std::istream&
operator>>(std::istream& is, dynamic_bitset<Block, Allocator>& b)
{
    typedef char CharT;
    std::string buf;
    typedef typename std::string::size_type size_type;
    const size_type N = b.size();
    buf.reserve(N);

    // skip whitespace
    if (is.flags() & std::ios::skipws) {
        char c;
        do
            is.get(c);
#if defined(__GNUC__) && !defined(__SGI_STL_PORT)
        while (is && isspace(c));
#else
        while (is && std::isspace(c));
#endif
        if (is)
            is.putback(c);
    }

    size_type i;
    for (i = 0; i < N; ++i)
    {
        CharT c;
        is.get(c);
        if (c == '0' || c == '1')
            buf += c;
        else
        {
            is.putback(c);
            break;
        }
    }

    if (i == 0)
        is.clear(is.rdstate() | std::ios::failbit);
    else
    {
        dynamic_bitset<Block, Allocator> tmp(buf);
        b.swap(tmp);
    }
    return is;
}
#else
template <typename CharT, typename Traits, typename Block, typename Allocator>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& in_stream,
           dynamic_bitset<Block, Allocator>& b)
{
    std::basic_string<CharT,Traits> tmp;
    typedef typename std::basic_string<CharT,Traits>::size_type size_type;
    const size_type N = b.size();
    tmp.reserve(N);

    // skip whitespace
    typename std::basic_istream<CharT, Traits>::sentry sentry(in_stream);
    if (sentry) {
      std::basic_streambuf<CharT, Traits>* read_buf = in_stream.rdbuf();

      for (size_type i = 0; i < N; ++i) {
        typename Traits::int_type c1 = read_buf->sbumpc();
        if (Traits::eq_int_type(c1, Traits::eof())) {
          in_stream.setstate(std::ios_base::eofbit);
          break;
        } else {
          char c2 = Traits::to_char_type(c1);
          char c  = in_stream.narrow(c2, '*');

          if (c == '0' || c == '1')
            tmp += c; // old dinkumware basic_string missing push_back
          else if (Traits::eq_int_type(read_buf->sputbackc(c2), Traits::eof()))
          {
            in_stream.setstate(std::ios_base::failbit);
            break;
          }
        }
      } // for

      if (tmp.empty()) // did not read in enough bits
        in_stream.setstate(std::ios_base::failbit);
      else
          b.from_string(tmp, static_cast<size_type>(0), N);
    } // if (sentry)
    return in_stream;
}
#endif


//-----------------------------------------------------------------------------
// bitset operations

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator&(const dynamic_bitset<Block, Allocator>& x,
          const dynamic_bitset<Block, Allocator>& y)
{
    dynamic_bitset<Block, Allocator> b(x);
    return b &= y;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator|(const dynamic_bitset<Block, Allocator>& x,
          const dynamic_bitset<Block, Allocator>& y)
{
    dynamic_bitset<Block, Allocator> b(x);
    return b |= y;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator^(const dynamic_bitset<Block, Allocator>& x,
          const dynamic_bitset<Block, Allocator>& y)
{
    dynamic_bitset<Block, Allocator> b(x);
    return b ^= y;
}

template <typename Block, typename Allocator>
dynamic_bitset<Block, Allocator>
operator-(const dynamic_bitset<Block, Allocator>& x,
          const dynamic_bitset<Block, Allocator>& y)
{
    dynamic_bitset<Block, Allocator> b(x);
    return b -= y;
}


//-----------------------------------------------------------------------------
// private member functions

template <typename Block, typename Allocator>
inline void dynamic_bitset<Block, Allocator>::
set_(size_type bit)
{
    this->m_bits[this->word(bit)] |= this->mask1(bit);
}

template <typename Block, typename Allocator>
inline void dynamic_bitset<Block, Allocator>::
set_block_(size_type blocknum, Block value)
{
    this->m_bits[blocknum] = value;
}

template <typename Block, typename Allocator>
inline void dynamic_bitset<Block, Allocator>::
reset_(size_type b)
{
    this->m_bits[this->word(b)] &= this->mask0(b);
}

template <typename Block, typename Allocator>
inline bool dynamic_bitset<Block, Allocator>::test_(size_type b) const
{
    return (this->m_bits[this->word(b)] & this->mask1(b)) != static_cast<Block>(0);
}

template <typename Block, typename Allocator>
bool dynamic_bitset<Block, Allocator>::set_(size_type n, bool value)
{
    if (value)
        set_(n);
    else
        reset_(n);
    return value != static_cast<Block>(0);
}


// If size() is not a multiple of bits_per_block
// then not all the bits in the last block are used.
// This function resets the unused bits (convenient
// for the implementation of many member functions)
//
template <typename Block, typename Allocator>
inline void dynamic_bitset<Block, Allocator>::m_zero_unused_bits()
{
    assert (this->m_num_blocks == this->calc_num_blocks(this->m_num_bits));

    // if != 0 this is the number of bits used in the last block
    const size_type used_bits = this->m_num_bits % bits_per_block;

    if (used_bits != 0)
        this->m_bits[this->m_num_blocks - 1] &= ~(~static_cast<Block>(0) << used_bits);

}


} // namespace boost


#undef BOOST_WORKAROUND_REPEAT_DEFAULT_TEMPLATE_ARGUMENTS

#endif // BOOST_DYNAMIC_BITSET_HPP


