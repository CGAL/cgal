// -*- C++ -*-
//  Boost general library format ---------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2001
//                  krempp@crans.ens-cachan.fr
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ideas taken from Rüdiger Loos's format class
// and Karl Nelson's ofstream

// ----------------------------------------------------------------------------
// format_implementation.hpp  Implementation of the basic_format class
// ----------------------------------------------------------------------------


#ifndef BOOST_FORMAT_IMPLEMENTATION_HPP
#define BOOST_FORMAT_IMPLEMENTATION_HPP

#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>
#include <boost/format/format_class.hpp>
#include <algorithm> // std::swap

namespace boost {

// --------  format:: -------------------------------------------
    template< class Ch, class Tr>
    basic_format<Ch, Tr>:: basic_format(const Ch* str)
        : style_(0), cur_arg_(0), num_args_(0), dumped_(false),
          exceptions_(io::all_error_bits)
    {
        if( str)
            parse( str );
    }

#ifndef BOOST_NO_STD_LOCALE
    template< class Ch, class Tr>
    basic_format<Ch, Tr>:: basic_format(const Ch* str, const std::locale & loc)
        : style_(0), cur_arg_(0), num_args_(0), dumped_(false),
          exceptions_(io::all_error_bits)
    {
        oss_.imbue( loc );
        if(str) parse( str );
    }

    template< class Ch, class Tr>
    basic_format<Ch, Tr>:: basic_format(const string_t& s, const std::locale & loc)
        : style_(0), cur_arg_(0), num_args_(0), dumped_(false),
          exceptions_(io::all_error_bits)
    {
        oss_.imbue( loc );
        parse(s);  
    }
#endif //BOOST_NO_STD_LOCALE

    template< class Ch, class Tr>
    basic_format<Ch, Tr>:: basic_format(const string_t& s)
        : style_(0), cur_arg_(0), num_args_(0), dumped_(false),
          exceptions_(io::all_error_bits)
    {
        parse(s);  
    }

    template< class Ch, class Tr>
    basic_format<Ch, Tr>:: basic_format(const basic_format& x)
        : items_(x.items_), bound_(x.bound_), style_(x.style_), 
          cur_arg_(x.cur_arg_), num_args_(x.num_args_), dumped_(false), 
          prefix_(x.prefix_), exceptions_(x.exceptions_) 
    { 
    } 

    template< class Ch, class Tr>
    void  basic_format<Ch, Tr>:: swap (basic_format & x) {
        std::swap(exceptions_, x.exceptions_);
        std::swap(style_, x.style_); 
        std::swap(cur_arg_, x.cur_arg_); 
        std::swap(num_args_, x.num_args_);
        std::swap(dumped_, x.dumped_);

        items_.swap(x.items_);
        prefix_.swap(x.prefix_);
        bound_.swap(x.bound_);
    }

    template< class Ch, class Tr>
    basic_format<Ch, Tr>& basic_format<Ch, Tr>:: operator= (const basic_format& x) {
        if(this == &x)
            return *this;
        (basic_format<Ch, Tr>(x)).swap(*this);
        return *this;
    }

    template< class Ch, class Tr>
    unsigned char basic_format<Ch,Tr>:: exceptions() const {
        return exceptions_; 
    }

    template< class Ch, class Tr>
    unsigned char basic_format<Ch,Tr>:: exceptions(unsigned char newexcept) { 
        unsigned char swp = exceptions_; 
        exceptions_ = newexcept; 
        return swp; 
    }

    template<class Ch, class Tr>
      void basic_format<Ch, Tr>:: make_or_reuse_data(std::size_t nbitems) {
        Ch fill = oss_.widen(' ');
        if(items_.size() == 0)
            items_.assign( nbitems, format_item_t(fill) );
        else {
            bound_.resize(0);
            items_.resize(nbitems, format_item_t(fill));
            for(std::size_t i=0; i < nbitems; ++i)
                items_[i].reset(fill); //  strings are resized to "", instead of reallocated
        }
    }

    template< class Ch, class Tr>
    basic_format<Ch,Tr>& basic_format<Ch,Tr>:: clear() {
        // empty the string buffers (except bound arguments)
        // and make the format object ready for formatting a new set of arguments

        BOOST_ASSERT( bound_.size()==0 || num_args_ == static_cast<int>(bound_.size()) );

        for(unsigned long i=0; i<items_.size(); ++i) {
            // clear converted strings only if the corresponding argument is not  bound :
            if( bound_.size()==0 || !bound_[ items_[i].argN_ ] )  items_[i].res_.resize(0);
        }
        cur_arg_=0; dumped_=false;
        // maybe first arg is bound:
        if(bound_.size() != 0) {
            while(cur_arg_ < num_args_ && bound_[cur_arg_] )
                  ++cur_arg_;
        }
        return *this;
    }

    template< class Ch, class Tr>
    basic_format<Ch,Tr>& basic_format<Ch,Tr>:: clear_binds() {
        // remove all binds, then clear()
        bound_.resize(0);
        clear();
        return *this;
    }

    template< class Ch, class Tr>
    basic_format<Ch,Tr>& basic_format<Ch,Tr>:: clear_bind(int argN) {
        // remove the bind of ONE argument then clear()

        if(argN<1 || argN > num_args_ || bound_.size()==0 || !bound_[argN-1] ) {
            if( exceptions() & io::out_of_range_bit )
                boost::throw_exception(io::out_of_range()); // arg not in range.
            else return *this;
        }
        bound_[argN-1]=false;
        clear();
        return *this;
    }



    template< class Ch, class Tr>
    std::basic_string<Ch,Tr> basic_format<Ch,Tr>:: str() const {
        dumped_=true;
        if(items_.size()==0)
            return prefix_;
        if( cur_arg_ < num_args_)
            if( exceptions() & io::too_few_args_bit )
                boost::throw_exception(io::too_few_args()); // not enough variables supplied
        unsigned long i;
        string_t res;
        res.reserve(size());
        res += prefix_;
        for(i=0; i < items_.size(); ++i) {
            const format_item_t& item = items_[i];
            res += item.res_;
            if( item.argN_ == format_item_t::argN_tabulation) { 
                BOOST_ASSERT( item.pad_scheme_ & format_item_t::tabulation);
                std::streamsize  n = item.fmtstate_.width_ - res.size();
                if( n > 0 )
                    res.append( n, item.fmtstate_.fill_ );
            }
            res += item.appendix_;
        }
        return res;
    }
    template< class Ch, class Tr>
    typename basic_format<Ch, Tr>::size_type  basic_format<Ch,Tr>:: 
    size () const {
        std::streamsize sz = prefix_.size();
        unsigned long i;
        for(i=0; i < items_.size(); ++i) {
            const format_item_t& item = items_[i];
            sz += item.res_.size();
            if( item.argN_ == format_item_t::argN_tabulation)
                sz = std::max(sz, item.fmtstate_.width_);
            sz +=  + item.appendix_.size();
        }
        return static_cast<size_type> (sz);
    }

namespace io {
namespace detail {

    template<class Ch, class Tr, class T> 
    basic_format<Ch, Tr>&  bind_arg_body( basic_format<Ch, Tr>& self, 
                                          int argN, const T& val) {
        // bind one argument to a fixed value
        // this is persistent over clear() calls, thus also over str() and <<
        if(self.dumped_) 
            self.clear(); // needed because we will modify cur_arg_
        if(argN<1 || argN > self.num_args_) {
            if( self.exceptions() & io::out_of_range_bit )
                boost::throw_exception(io::out_of_range()); // arg not in range.
            else return self;
        }
        if(self.bound_.size()==0) 
            self.bound_.assign(self.num_args_,false);
        else 
            BOOST_ASSERT( self.num_args_ == static_cast<signed int>(self.bound_.size()) );
        int o_cur_arg = self.cur_arg_;
        self.cur_arg_ = argN-1; // arrays begin at 0

        self.bound_[self.cur_arg_]=false; // if already set, we unset and re-sets..
        self.operator%(val); // put val at the right place, because cur_arg is set
    

        // Now re-position cur_arg before leaving :
        self.cur_arg_ = o_cur_arg; 
        self.bound_[argN-1]=true;
        if(self.cur_arg_ == argN-1 ) {
            // hum, now this arg is bound, so move to next free arg
            while(self.cur_arg_ < self.num_args_ && self.bound_[self.cur_arg_])   
                ++self.cur_arg_;
        }
        // In any case, we either have all args, or are on a non-binded arg :
        BOOST_ASSERT( self.cur_arg_ >= self.num_args_ || ! self.bound_[self.cur_arg_]);
        return self;
    }

    template<class Ch, class Tr, class T> 
    basic_format<Ch, Tr>&  modify_item_body( basic_format<Ch, Tr>& self,
                                             int itemN, T manipulator) {
        // applies a manipulator to the format_item describing a given directive.
        // this is a permanent change, clear or reset won't cancel that.
        if(itemN<1 || itemN > static_cast<signed int>(self.items_.size() )) {
            if( self.exceptions() & io::out_of_range_bit ) 
                boost::throw_exception(io::out_of_range()); // item not in range.
            else return self;
        }
        self.items_[itemN-1].fmtstate_. template apply_manip<T> ( manipulator );
        return self;
    }

} // namespace detail
} // namespace io
} // namespace boost



#endif  // BOOST_FORMAT_IMPLEMENTATION_HPP
