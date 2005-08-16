/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_EXPRESSION_VALUE_HPP_452FE66D_8754_4107_AF1E_E42255A0C18A_INCLUDED)
#define CPP_EXPRESSION_VALUE_HPP_452FE66D_8754_4107_AF1E_E42255A0C18A_INCLUDED

#if defined (BOOST_SPIRIT_DEBUG)
#include <iostream>
#endif // defined(BOOST_SPIRIT_DEBUG)

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace grammars {
namespace closures {

///////////////////////////////////////////////////////////////////////////////
//
//  The closure_value class represents the closure type, which is used for the 
//  expression grammar. 
//
//      This class was introduced to allow the expression grammar to respect 
//      the numeric type of a numeric literal or expression result.
//
///////////////////////////////////////////////////////////////////////////////
class closure_value {
public:

    enum value_type {
        is_int = 1,
        is_uint = 2,
        is_bool = 3
    };
    
    enum value_error {
        error_noerror = 0x0,
        error_division_by_zero = 0x1,
        error_overflow = 0x2
    };
    
    closure_value(value_error valid_ = error_noerror) 
    : type(is_int), valid(valid_) 
    { value.i = 0; }
    explicit closure_value(int i, value_error valid_ = error_noerror) 
    : type(is_int), valid(valid_) 
    { value.i = i; }
    explicit closure_value(unsigned int ui, value_error valid_ = error_noerror) 
    : type(is_uint), valid(valid_) 
    { value.ui = ui; }
    explicit closure_value(long i, value_error valid_ = error_noerror) 
    : type(is_int), valid(valid_) 
    { value.i = i; }
    explicit closure_value(unsigned long ui, value_error valid_ = error_noerror) 
    : type(is_uint), valid(valid_) 
    { value.ui = ui; }
    explicit closure_value(bool b, value_error valid_ = error_noerror) 
    : type(is_bool), valid(valid_) 
    { value.b = b; }

    value_type get_type() const { return type; }
    value_error is_valid() const { return valid; }
    
// implicit conversion
    operator int() const 
    {
        switch (type) {
        case is_uint:   return value.ui;
        case is_bool:   return value.b ? 1 : 0;
        case is_int:    break;
        }
        return value.i;
    }
    operator unsigned int() const 
    {
        switch (type) {
        case is_uint:   return value.ui;
        case is_bool:   return value.b ? 1 : 0;
        case is_int:    break;
        }
        return value.i;
    }
    operator long() const 
    {
        switch (type) {
        case is_uint:   return value.ui;
        case is_bool:   return value.b ? 1 : 0;
        case is_int:    break;
        }
        return value.i;
    }
    operator unsigned long() const 
    {
        switch (type) {
        case is_uint:   return value.ui;
        case is_bool:   return value.b ? 1 : 0;
        case is_int:    break;
        }
        return value.i;
    }
    operator bool() const 
    {
        switch (type) {
        case is_uint:   return value.ui != 0;
        case is_bool:   return value.b;
        case is_int:    break;
        }
        return value.i != 0.0;
    }

// assignment    
    closure_value &operator= (closure_value const &rhs)
    {
        switch (rhs.get_type()) {
        case is_int:    
            value.i = long(rhs); 
            type = is_int;
            break;
        
        case is_uint:   
            value.ui = (unsigned long)(rhs); 
            type = is_uint;
            break;
        
        case is_bool:   
            value.b = bool(rhs);
            type = is_bool;
            break;
        }
        valid = rhs.valid;
        return *this;
    }
    closure_value &operator= (int rhs)
    {
        type = is_int;
        value.i = rhs;
        valid = error_noerror;
        return *this;
    }
    closure_value &operator= (unsigned int rhs)
    {
        type = is_uint;
        value.ui = rhs;
        valid = error_noerror;
        return *this;
    }
    closure_value &operator= (long rhs)
    {
        type = is_int;
        value.i = rhs;
        valid = error_noerror;
        return *this;
    }
    closure_value &operator= (unsigned long rhs)
    {
        type = is_uint;
        value.ui = rhs;
        valid = error_noerror;
        return *this;
    }
    closure_value &operator= (bool rhs)
    {
        type = is_bool;
        value.b = rhs;
        valid = error_noerror;
        return *this;
    }

// arithmetics
    closure_value &operator+= (closure_value const &rhs)
    {
        switch (type) {
        case is_int:    
            switch(rhs.type) {
            case is_bool:
                {
                    long result = value.i + long(rhs); 
                    if (rhs.value.i > 0L && value.i > result || 
                        rhs.value.i < 0L && value.i < result)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.i = result;
                    }
                }
                break;
                
            case is_int:
                {
                    long result = value.i + rhs.value.i;
                    if (rhs.value.i > 0L && value.i > result || 
                        rhs.value.i < 0L && value.i < result)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.i = result;
                    }
                }
                break;
                
            case is_uint:
                {
                    unsigned long result = value.ui + rhs.value.ui; 
                    if (result < value.ui) {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                        type = is_uint; 
                    }
                }
                break;
            }
            break;
            
        case is_uint:
            {
                unsigned long result = value.ui + (unsigned long)(rhs); 
                if (result < value.ui) {
                    valid = error_overflow;
                }
                else {
                    value.ui = result;
                }
            }
            break;
            
        case is_bool:   
            value.i = value.b + bool(rhs);
            type = is_int;
        }
        valid = (value_error)(valid | rhs.valid);
        return *this;
    }
    closure_value &operator-= (closure_value const &rhs)
    {
        switch (type) {
        case is_int:
            switch(rhs.type) {
            case is_bool:
                {
                    long result = value.i - long(rhs); 
                    if (rhs.value.i > 0L && result > value.i || 
                        rhs.value.i < 0L && result < value.i)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.i = result;
                    }
                }
                break;

            case is_int:
                {
                    long result = value.i - rhs.value.i;
                    if (rhs.value.i > 0L && result > value.i || 
                        rhs.value.i < 0L && result < value.i)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.i = result;
                    }
                }
                break;
                
            case is_uint:
                {
                    unsigned long result = value.ui - rhs.value.ui; 
                    if (result > value.ui) {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                        type = is_uint; 
                    }
                }
                break;
            }
            break;
            
        case is_uint:
            switch(rhs.type) {
            case is_bool:
                {
                    unsigned long result = value.ui - (unsigned long)(rhs); 
                    if (result > value.ui)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                    }
                }
                break;

            case is_int:
                {
                    unsigned long result = value.ui - rhs.value.i;
                    if (rhs.value.i > 0L && result > value.ui || 
                        rhs.value.i < 0L && result < value.ui)
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                    }
                }
                break;
                
            case is_uint:
                {
                    unsigned long result = value.ui - rhs.value.ui; 
                    if (result > value.ui) {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                    }
                }
                break;
            }
            break;

        case is_bool:   
            value.i = value.b - bool(rhs);
            type = is_int;
        }
        valid = (value_error)(valid | rhs.valid);
        return *this;
    }
    closure_value &operator*= (closure_value const &rhs)
    {
        switch (type) {
        case is_int:    
            switch(rhs.type) {
            case is_bool:   value.i *= long(rhs); break;
            case is_int:
                {
                    long result = value.i * rhs.value.i; 
                    if (0 != value.i && 0 != rhs.value.i &&
                        (result / value.i != rhs.value.i ||
                         result / rhs.value.i != value.i)
                       )
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.i = result;
                    }
                }
                break;
                
            case is_uint:
                {
                    unsigned long result = value.ui * rhs.value.ui; 
                    if (0 != value.ui && 0 != rhs.value.ui &&
                        (result / value.ui != rhs.value.ui ||
                         result / rhs.value.ui != value.ui)
                       )
                    {
                        valid = error_overflow;
                    }
                    else {
                        value.ui = result;
                        type = is_uint; 
                    }
                }
                break;
            }
            break;
            
        case is_uint:
            {
                unsigned long rhs_val = (unsigned long)(rhs);
                unsigned long result = value.ui * rhs_val; 
                if (0 != value.ui && 0 != rhs_val &&
                    (result / value.ui != rhs_val ||
                      result / rhs_val != value.ui)
                    )
                {
                    valid = error_overflow;
                }
                else {
                    value.ui = result;
                    type = is_uint; 
                }
            }
            break;
            
        case is_bool:
            switch (rhs.type) {
            case is_int:
                value.i = (value.b ? 1 : 0) * rhs.value.i; 
                type = is_int; 
                break;
                
            case is_uint:
                value.ui = (value.b ? 1 : 0) * rhs.value.ui; 
                type = is_uint; 
                break;
                
            case is_bool:
                value.b = 0 != ((value.b ? 1 : 0) * (rhs.value.b ? 1 : 0));
                break;
            }
        }
        valid = (value_error)(valid | rhs.valid);
        return *this;
    }
    closure_value &operator/= (closure_value const &rhs)
    {
        switch (type) {
        case is_int:    
            switch(rhs.type) {
            case is_bool:   
            case is_int:
                if (long(rhs) != 0) {
                    if (value.i == -value.i && -1 == rhs.value.i) {
                    // LONG_MIN / -1 on two's complement
                        valid = error_overflow;
                    }
                    else {
                        value.i /= long(rhs); 
                    }
                }
                else {
                    valid = error_division_by_zero;   // division by zero
                }
                break;
                
            case is_uint:
                if (rhs.value.ui != 0) {
                    value.ui /= rhs.value.ui; 
                    type = is_uint; 
                }
                else {
                    valid = error_division_by_zero;      // division by zero
                }
                break;
            }
            break;
            
        case is_uint: 
            if ((unsigned long)(rhs) != 0) 
                value.ui /= (unsigned long)(rhs); 
            else
                valid = error_division_by_zero;         // division by zero
            break;

        case is_bool:  
            if (bool(rhs)) {
                switch(rhs.type) {
                case is_int:
                    value.i = (value.b ? 1 : 0) / rhs.value.i;
                    type = is_int;
                    break;
                    
                case is_uint:
                    value.i = (value.b ? 1 : 0) / rhs.value.ui;
                    type = is_int;
                    break;
                    
                case is_bool:
                    break;
                }
            }
            else {
                valid = error_division_by_zero;         // division by zero
            }
        }
        return *this;
    }
    closure_value &operator%= (closure_value const &rhs)
    {
        switch (type) {
        case is_int:    
            switch(rhs.type) {
            case is_bool:   
            case is_int:
                if (long(rhs) != 0) {
                    if (value.i == -value.i && -1 == rhs.value.i) {
                    // LONG_MIN % -1 on two's complement
                        valid = error_overflow;
                    }
                    else {
                        value.i %= long(rhs); 
                    }
                }
                else {
                    valid = error_division_by_zero;      // division by zero
                }
                break;
                
            case is_uint:
                if (rhs.value.ui != 0) {
                    value.ui %= rhs.value.ui; 
                    type = is_uint; 
                }
                else {
                    valid = error_division_by_zero;      // division by zero
                }
                break;
            }
            break;
            
        case is_uint: 
            if ((unsigned long)(rhs) != 0) 
                value.ui %= (unsigned long)(rhs); 
            else
                valid = error_division_by_zero;      // division by zero
            break;

        case is_bool:  
            if (bool(rhs)) {
                switch(rhs.type) {
                case is_int:
                    value.i = (value.b ? 1 : 0) % rhs.value.i;
                    type = is_int;
                    break;
                    
                case is_uint:
                    value.i = (value.b ? 1 : 0) % rhs.value.ui;
                    type = is_int;
                    break;
                    
                case is_bool:
                    break;
                }                    
            }
            else {
                valid = error_division_by_zero;      // division by zero
            }
        }
        return *this;
    }

    friend closure_value 
    operator- (closure_value const &rhs)
    {
        switch (rhs.type) {
        case is_int:
            {
                long value = long(rhs);
                if (value != 0 && value == -value)
                    return closure_value(-value, error_overflow);
                return closure_value(-value, rhs.valid);
            }
            
        case is_bool:   return closure_value(-long(rhs), rhs.valid); 
        case is_uint:   break;
        }

        long value = (unsigned long)(rhs);
        if (value != 0 && value == -value)
            return closure_value(-value, error_overflow);
        return closure_value(-value, rhs.valid);
    }
    friend closure_value 
    operator~ (closure_value const &rhs)
    {
        return closure_value(~(unsigned long)(rhs), rhs.valid);
    }
    friend closure_value 
    operator! (closure_value const &rhs)
    {
        switch (rhs.type) {
        case is_int:    return closure_value(!long(rhs), rhs.valid);
        case is_bool:   return closure_value(!bool(rhs), rhs.valid); 
        case is_uint:   break;
        }
        return closure_value(!(unsigned long)(rhs), rhs.valid);
    }
    
// comparison
    friend closure_value 
    operator== (closure_value const &lhs, closure_value const &rhs)
    {
        bool cmp = false;
        switch (lhs.type) {
        case is_int:
            switch(rhs.type) {
            case is_bool:   cmp = bool(lhs) == rhs.value.b; break;
            case is_int:    cmp = lhs.value.i == rhs.value.i; break;
            case is_uint:   cmp = lhs.value.ui == rhs.value.ui; break;
            }
            break;
            
        case is_uint:   cmp = lhs.value.ui == (unsigned long)(rhs); break;
        case is_bool:   cmp = lhs.value.b == bool(rhs); break;
        }
        return closure_value(cmp, (value_error)(lhs.valid & rhs.valid));
    }
    friend closure_value 
    operator!= (closure_value const &lhs, closure_value const &rhs)
    {
        return closure_value(!bool(lhs == rhs), (value_error)(lhs.valid & rhs.valid));
    }
    friend closure_value 
    operator> (closure_value const &lhs, closure_value const &rhs)
    {
        bool cmp = false;
        switch (lhs.type) {
        case is_int:
            switch(rhs.type) {
            case is_bool:   cmp = lhs.value.i > long(rhs); break;
            case is_int:    cmp = lhs.value.i > rhs.value.i; break;
            case is_uint:   cmp = lhs.value.ui > rhs.value.ui; break;
            }
            break;
            
        case is_uint:   cmp = lhs.value.ui > (unsigned long)(rhs); break;
        case is_bool:   cmp = lhs.value.b > bool(rhs); break;
        }
        return closure_value(cmp, (value_error)(lhs.valid & rhs.valid));
    }
    friend closure_value 
    operator< (closure_value const &lhs, closure_value const &rhs)
    {
        bool cmp = false;
        switch (lhs.type) {
        case is_int:    cmp = long(lhs) < long(rhs); break;
            switch(rhs.type) {
            case is_bool:   cmp = lhs.value.i < long(rhs); break;
            case is_int:    cmp = lhs.value.i < rhs.value.i; break;
            case is_uint:   cmp = lhs.value.ui < rhs.value.ui; break;
            }
            break;
            
        case is_uint:   cmp = lhs.value.ui < (unsigned long)(rhs); break;
        case is_bool:   cmp = bool(lhs) < bool(rhs); break;
        }
        return closure_value(cmp, (value_error)(lhs.valid & rhs.valid));
    }
    friend closure_value 
    operator<= (closure_value const &lhs, closure_value const &rhs)
    {
        return closure_value(!bool(lhs > rhs), (value_error)(lhs.valid & rhs.valid));
    }
    friend closure_value 
    operator>= (closure_value const &lhs, closure_value const &rhs)
    {
        return closure_value(!bool(lhs < rhs), (value_error)(lhs.valid & rhs.valid));
    }

    closure_value &
    operator<<= (closure_value const &rhs)
    {
        switch (type) {
        case is_bool:
        case is_int:
            switch (rhs.type) {
            case is_bool:
            case is_int:
                {
                long shift_by = long(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    else if (shift_by < -64)
                        shift_by = -64;
                    value.i <<= shift_by; 
                }
                break;
                
            case is_uint:
                {
                unsigned long shift_by = (unsigned long)(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    value.ui <<= shift_by; 
                
                // Note: The usual arithmetic conversions are not performed on 
                //       bit shift operations.
                }
                break;
            }
            break;

        case is_uint:
            switch (rhs.type) {
            case is_bool:
            case is_int:
                {
                long shift_by = long(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    else if (shift_by < -64)
                        shift_by = -64;
                    value.ui <<= shift_by; 
                }
                break;
                
            case is_uint:
                {
                unsigned long shift_by = (unsigned long)(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    value.ui <<= shift_by; 
                }
                break;
            }
        }
        valid = (value_error)(valid | rhs.valid);
        return *this;
    }

    closure_value &
    operator>>= (closure_value const &rhs)
    {
        switch (type) {
        case is_bool:
        case is_int:
            switch (rhs.type) {
            case is_bool:
            case is_int:
                {
                long shift_by = long(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    else if (shift_by < -64)
                        shift_by = -64;
                    value.i >>= shift_by; 
                }
                break;
                
            case is_uint:
                {
                unsigned long shift_by = (unsigned long)(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    value.ui >>= shift_by; 
                
                // Note: The usual arithmetic conversions are not performed on 
                //       bit shift operations.
                }
                break;
            }
            break;
            
        case is_uint:
            switch (rhs.type) {
            case is_bool:
            case is_int:
                {
                long shift_by = long(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    else if (shift_by < -64)
                        shift_by = -64;
                    value.ui >>= shift_by; 
                }
                break;
                
            case is_uint:
                {
                unsigned long shift_by = (unsigned long)(rhs);
                    
                    if (shift_by > 64) 
                        shift_by = 64;
                    value.ui >>= shift_by; 
                }
                break;
            }
            break;
        }
        valid = (value_error)(valid | rhs.valid);
        return *this;
    }

    friend closure_value 
    operator|| (closure_value const &lhs, closure_value const &rhs)
    {
        bool result = bool(lhs) || bool(rhs);
        return closure_value(result, (value_error)(lhs.valid & rhs.valid));
    }
    
    friend closure_value 
    operator&& (closure_value const &lhs, closure_value const &rhs)
    {
        bool result = bool(lhs) && bool(rhs);
        return closure_value(result, (value_error)(lhs.valid & rhs.valid));
    }

    // handle the ?: operator
    closure_value &
    handle_questionmark(closure_value const &cond, closure_value const &val2)
    {
        switch (type) {
        case is_int:
            switch (val2.type) {
            case is_bool: value.b = bool(cond) ? value.b : bool(val2); break;
            case is_int:  value.i = bool(cond) ? value.i : long(val2); break;
            case is_uint: 
                value.ui = bool(cond) ? value.ui : (unsigned long)(val2); 
                type = is_uint;   // changing type!
                break;
            }
            break;
            
        case is_uint:   value.ui = bool(cond) ? value.ui : (unsigned long)(val2); break;
        case is_bool:   value.b = bool(cond) ? value.b : bool(val2); break;
        }
        valid = bool(cond) ? valid : val2.valid;
        return *this;
    }
    
#if defined (BOOST_SPIRIT_DEBUG)
    friend std::ostream&
    operator<< (std::ostream &o, closure_value const &val)
    {
        switch (val.type) {
        case is_int:    o << "int(" << long(val) << ")"; break;
        case is_uint:   o << "unsigned int(" << (unsigned long)(val) << ")"; break;
        case is_bool:   o << "bool(" << bool(val) << ")"; break;
        }
        return o;
    }
#endif // defined(BOOST_SPIRIT_DEBUG)

private:
    value_type type;
    union {
        long i;
        unsigned long ui;
        bool b;
    } value;
    value_error valid;
};

///////////////////////////////////////////////////////////////////////////////
}   // namespace closures
}   // namespace grammars
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_EXPRESSION_VALUE_HPP_452FE66D_8754_4107_AF1E_E42255A0C18A_INCLUDED)
