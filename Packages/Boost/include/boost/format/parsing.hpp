// -*- C++ -*-
//  Boost general library 'format'   ---------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2001
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ideas taken from Rudiger Loos's format class
// and Karl Nelson's ofstream (also took its parsing code as basis for printf parsing)

// ------------------------------------------------------------------------------
// parsing.hpp :  implementation of the parsing member functions
//                      ( parse, parse_printf_directive)
// ------------------------------------------------------------------------------


#ifndef BOOST_FORMAT_PARSING_HPP
#define BOOST_FORMAT_PARSING_HPP


#include <boost/format/format_class.hpp>
#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>


namespace boost {
namespace io {
namespace detail {

    template<class Ch, class Stream> inline
    bool wrap_isdigit(Ch c, Stream &os) {
#if ! defined( BOOST_NO_LOCALE_ISIDIGIT )
        return std::isdigit(c, os.rdbuf()->getloc() );
# else
        using namespace std;
        return isdigit(c); 
#endif 
    } //end- wrap_isdigit(..)
 
    template<class Iter, class Stream> 
    Iter wrap_scan_notdigit(Iter beg, Iter end, const Stream & os) {
        using namespace std;
        for( ; beg!=end && wrap_isdigit(*beg,os); ++beg) ;
        return beg;
    }


    template<class Res, class Iter, class Stream>
    Iter str2int(const Iter & start, const Iter & last, Res & res, Stream &os) 
        // Input : [start, last) iterators range and a
        //          a basic_ios& merely to use its widen/narrow member function
        // Effects : reads sequence and converts digits into an integral n, of type Res
        // Returns : n
    {
        using namespace std;
        Iter it;
        res=0;
        for(it=start; it != last && wrap_isdigit(*it, os); ++it ) {
            char cur_ch = os.narrow( *it, 0); // cant fail.
            res *= 10;
            res += cur_ch - '0'; // 22.2.1.1.2.13 of the C++ standard
        }
        return it;
    }

    template<class Iter, class Stream>
    Iter skip_asterisk(Iter start, Iter last, Stream &os) 
        // skip printf's "asterisk-fields" directives in the format-string buf
        // Input : char string, with starting index *pos_p
        //         a basic_ios& merely to use its widen/narrow member function
        // Effects : advance *pos_p by skipping printf's asterisk fields.
        // Returns : nothing
    {
        using namespace std;
        ++ start;
        start = wrap_scan_notdigit(start, last, os);
        if(start!=last && *start== os.widen('$') )
            ++start;
        return start;
    }


    inline void maybe_throw_exception( unsigned char exceptions)
        // auxiliary func called by parse_printf_directive
        // for centralising error handling
        // it either throws if user sets the corresponding flag, or does nothing.
    {
        if(exceptions & io::bad_format_string_bit)
            boost::throw_exception(io::bad_format_string());
    }
    


    template<class Ch, class Tr, class Iter, class Stream>
    bool parse_printf_directive(Iter & start, const Iter& last, 
                                detail::format_item<Ch, Tr> * fpar,
                                Stream &os,
                                unsigned char exceptions)
        // Input: a 'printf-directive' in the format-string, starting at buf[ *pos_p ]
        //        a basic_ios& merely to use its widen/narrow member function
        //        a bitset'excpetions' telling whether to throw exceptions on errors.
        // Returns: true if parse somehow succeeded (ignore some errors if exceptions disabled)
        //          false if it failed so bad that the directive should be printed verbatim
        // Effects:  *pos_p is incremented so that buf[*pos_p] is the first char after the directive
        //           *fpar is set with the parameters read in the directive
    {
        typedef format_item<Ch, Tr>  format_item_t;
        //BOOST_ASSERT( pos_p != 0);

        fpar->argN_ = format_item_t::argN_no_posit;  // if no positional-directive
        bool precision_set = false;
        bool in_brackets=false;
        if(*start== os.widen('|')) {
            in_brackets=true;
            if( ++start >= last ) {
                maybe_throw_exception(exceptions);
                return false;
            }
        }

        // the flag '0' would be picked as a digit for argument order, but here it's a flag :
        if(*start== os.widen('0')) 
            goto parse_flags;

        // handle argument order (%2$d)  or possibly width specification: %2d
        if(wrap_isdigit(*start, os)) {
            int n;
            start = str2int(start, last, n, os);
            if( start >= last ) {
                maybe_throw_exception(exceptions);
                return false;
            }
            
            // %N% case : this is already the end of the directive
            if( *start ==  os.widen('%') ) {
                fpar->argN_ = n-1;
                ++start;
                if( in_brackets) 
                    maybe_throw_exception(exceptions); 
                // but don't return.  maybe "%" was used in lieu of '$', so we go on.
                else
                    return true;
            }

            if ( *start== os.widen('$') ) {
                fpar->argN_ = n-1;
                ++start;
            } 
            else {
                // non-positionnal directive
                fpar->fmtstate_.width_ = n;
                fpar->argN_  = format_item_t::argN_no_posit;
                goto parse_precision;
            }
        }
    
      parse_flags: 
        // handle flags
        while ( start != last) { // as long as char is one of + - = _ # 0 l h   or ' '
            // misc switches
            switch ( os.narrow(*start, 0)) {
            case '\'' : break; // no effect yet. (painful to implement)
            case 'l':
            case 'h':  // short/long modifier : for printf-comaptibility (no action needed)
                break;
            case '-':
                fpar->fmtstate_.flags_ |= std::ios_base::left;
                break;
            case '=':
                fpar->pad_scheme_ |= format_item_t::centered;
                break;
            case '_':
                fpar->fmtstate_.flags_ |= std::ios_base::internal;
                break;
            case ' ':
                fpar->pad_scheme_ |= format_item_t::spacepad;
                break;
            case '+':
                fpar->fmtstate_.flags_ |= std::ios_base::showpos;
                break;
            case '0':
                fpar->pad_scheme_ |= format_item_t::zeropad;
                // need to know alignment before really setting flags,
                // so just add 'zeropad' flag for now, it will be processed later.
                break;
            case '#':
                fpar->fmtstate_.flags_ |= std::ios_base::showpoint | std::ios_base::showbase;
                break;
            default:
                goto parse_width;
            }
            ++start;
        } // loop on flag.

        if( start>=last) {
            maybe_throw_exception(exceptions);
            return true; 
        }
      parse_width:
        // handle width spec
        // first skip 'asterisk fields' :  *, or *N$
        if(*start == os.widen('*') )
            start = skip_asterisk(start, last, os); 
        if(start!=last && wrap_isdigit(*start, os))
            start = str2int(start, last, fpar->fmtstate_.width_, os);

      parse_precision:
        if( start>= last) { 
            maybe_throw_exception(exceptions);
            return true;
        }
        // handle precision spec
        if (*start== os.widen('.')) {
            ++start;
            if(start != last && *start == os.widen('*') )
                start = skip_asterisk(start, last, os); 
            if(start != last && wrap_isdigit(*start, os)) {
                start = str2int(start, last, fpar->fmtstate_.precision_, os);
                precision_set = true;
            }
            else
                fpar->fmtstate_.precision_ =0;
        }
    
        // handle  formatting-type flags :
        while( start != last && 
               ( *start== os.widen('l') || *start== os.widen('L') || *start== os.widen('h')) )
            ++start;
        if( start>=last) {
            maybe_throw_exception(exceptions);
            return true;
        }

        if( in_brackets && *start== os.widen('|') ) {
            ++start;
            return true;
        }
        switch ( os.narrow(*start, 0) ) {
        case 'X':
            fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
        case 'p': // pointer => set hex.
        case 'x':
            fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
            fpar->fmtstate_.flags_ |= std::ios_base::hex;
            break;

        case 'o':
            fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
            fpar->fmtstate_.flags_ |=  std::ios_base::oct;
            break;

        case 'E':
            fpar->fmtstate_.flags_ |=  std::ios_base::uppercase;
        case 'e':
            fpar->fmtstate_.flags_ &= ~std::ios_base::floatfield;
            fpar->fmtstate_.flags_ |=  std::ios_base::scientific;

            fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
            fpar->fmtstate_.flags_ |=  std::ios_base::dec;
            break;
      
        case 'f':
            fpar->fmtstate_.flags_ &= ~std::ios_base::floatfield;
            fpar->fmtstate_.flags_ |=  std::ios_base::fixed;
        case 'u':
        case 'd':
        case 'i':
            fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
            fpar->fmtstate_.flags_ |=  std::ios_base::dec;
            break;

        case 'T':
            ++start;
            if( start >= last)
                maybe_throw_exception(exceptions);
            else
                fpar->fmtstate_.fill_ = *start;
            fpar->pad_scheme_ |= format_item_t::tabulation;
            fpar->argN_ = format_item_t::argN_tabulation; 
            break;
        case 't': 
            fpar->fmtstate_.fill_ = os.widen(' ');
            fpar->pad_scheme_ |= format_item_t::tabulation;
            fpar->argN_ = format_item_t::argN_tabulation; 
            break;

        case 'G':
            fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
            break;
        case 'g': // 'g' conversion is default for floats.
            fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
            fpar->fmtstate_.flags_ |=  std::ios_base::dec;

            // CLEAR all floatield flags, so stream will CHOOSE
            fpar->fmtstate_.flags_ &= ~std::ios_base::floatfield; 
            break;

        case 'C':
        case 'c': 
            fpar->truncate_ = 1;
            break;
        case 'S':
        case 's': 
            if(precision_set) // handle truncation manually, with own parameter.
                fpar->truncate_ = fpar->fmtstate_.precision_;
            fpar->fmtstate_.precision_ = 6; // default stream precision.
            break;
        case 'n' :  
            fpar->argN_ = format_item_t::argN_ignored;
            break;
        default: 
            maybe_throw_exception(exceptions);
        }
        ++start;

        if( in_brackets ) {
            if( start != last && *start== os.widen('|') ) {
                ++start;
                return true;
            }
            else  maybe_throw_exception(exceptions);
        }
        return true;
    }


    template<class string_t, class Stream>
    int upper_bound_from_fstring(const string_t& buf, 
                                 const typename string_t::value_type arg_mark,
                                 Stream& os,  // just to carry the locale
                                 unsigned char exceptions) {
        // quick-parsing of the format-string to count arguments mark (arg_mark, '%')
        // returns : upper bound on the number of format items in the format strings
        typename string_t::size_type i1=0;
        int num_items=0;
        while( (i1=buf.find(arg_mark,i1)) != string_t::npos ) {
            if( i1+1 >= buf.size() ) {
                if(exceptions & io::bad_format_string_bit)
                    boost::throw_exception(io::bad_format_string()); // must not end in ".. %"
                else break; // stop there, ignore last '%'
            }
            if(buf[i1+1] == buf[i1] ) {// escaped "%%"
                i1+=2; continue; 
            }

            ++i1;
            // in case of %N% directives, dont count it double (wastes allocations..) :
            i1 = wrap_scan_notdigit(buf.begin()+i1, buf.end(), os) - buf.begin();
            if( i1 < buf.size() && buf[i1] == arg_mark )
                ++i1;
            ++num_items;
        }
        return num_items;
    }
} // detail namespace
} // io namespace



// -----------------------------------------------
//  format :: parse(..)

    template<class Ch, class Tr>
    basic_format<Ch, Tr>& basic_format<Ch, Tr>:: parse(const string_t& buf) {
        // parse the format-string 
        using namespace std;


        const Ch arg_mark = oss_.widen('%');
        bool ordered_args=true; 
        int max_argN=-1;

        // A: find upper_bound on num_items and allocates arrays
        int num_items = io::detail::upper_bound_from_fstring(buf, arg_mark, oss_, exceptions());
        make_or_reuse_data(num_items);

        // B: Now the real parsing of the format string :
        num_items=0;
        typename string_t::size_type i0=0, i1=0;
        typename string_t::const_iterator it;
        bool special_things=false;
        int cur_item=0;
        while( (i1=buf.find(arg_mark,i1)) != string_t::npos ) {
            string_t & piece = (cur_item==0) ? prefix_ : items_[cur_item-1].appendix_;
            if( buf[i1+1] == buf[i1] ) { // escaped mark, '%%' 
                piece += buf.substr(i0, i1+1-i0);
                i1+=2; i0=i1;
                continue; 
            }
            BOOST_ASSERT(  static_cast<unsigned int>(cur_item) < items_.size() || cur_item==0);

            if(i1!=i0)
                piece += buf.substr(i0, i1-i0);
            ++i1;
            it = buf.begin()+i1;
            bool parse_ok = io::detail::parse_printf_directive(
                it, buf.end(), &items_[cur_item], oss_, exceptions());
            i1 = it - buf.begin();
            if( ! parse_ok ) // the directive will be printed verbatim
                continue; 
            i0=i1;
            items_[cur_item].compute_states(); // process complex options, like zeropad, into params

            int argN=items_[cur_item].argN_;
            if(argN == format_item_t::argN_ignored)
                continue;
            if(argN ==format_item_t::argN_no_posit)
                ordered_args=false;
            else if(argN == format_item_t::argN_tabulation) special_things=true;
            else if(argN > max_argN) max_argN = argN;
            ++num_items;
            ++cur_item;
        } // loop on %'s
        BOOST_ASSERT(cur_item == num_items);
    
        // store the final piece of string
        {
            string_t & piece = (cur_item==0) ? prefix_ : items_[cur_item-1].appendix_;
            piece += buf.substr(i0);
        }
    
        if( !ordered_args) {
            if(max_argN >= 0 ) {  // dont mix positional with non-positionnal directives
                if(exceptions() & io::bad_format_string_bit)
                    boost::throw_exception(io::bad_format_string());
                // else do nothing. => positionnal arguments are processed as non-positionnal
            }
            // set things like it would have been with positional directives :
            int non_ordered_items = 0;
            for(int i=0; i< num_items; ++i)
                if(items_[i].argN_ == format_item_t::argN_no_posit) {
                    items_[i].argN_ = non_ordered_items;
                    ++non_ordered_items;
                }
            max_argN = non_ordered_items-1;
        }
    
        // C: set some member data :
        items_.resize(num_items, format_item_t(oss_.fill()) );

        if(special_things) style_ |= special_needs;
        num_args_ = max_argN + 1;
        if(ordered_args) style_ |=  ordered;
        else style_ &= ~ordered;
        return *this;
    }

} // namespace boost


#endif //  BOOST_FORMAT_PARSING_HPP
