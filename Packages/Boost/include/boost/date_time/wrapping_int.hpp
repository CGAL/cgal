#ifndef _DATE_TIME_WRAPPING_INT_HPP__
#define _DATE_TIME_WRAPPING_INT_HPP__

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */


namespace boost {
namespace date_time {

//! A wrapping integer used to support time durations 
/*! In composite date and time types this type is used to
 *  wrap at the day boundary.
 *
 */
template<typename int_type_, int_type_ wrap_val>
class wrapping_int {
public:
  typedef int_type_ int_type;
  //typedef overflow_type_ overflow_type;
  static int_type wrap_value() {return wrap_val;}
  //!Add, return true if wrapped
  wrapping_int(int_type v) : value_(v) {};
  //! Explicit converion method
  int_type as_int()   const   {return value_;}
  operator int_type() const   {return value_;}
  int_type add(int_type v) 
  {
    //take the mod here and assign it....
    int_type remainder = static_cast<int_type>(v % wrap_val);
    int_type overflow = static_cast<int_type>(v / wrap_val);
    value_ = static_cast<int_type>(value_ + remainder);
    if ((value_) >= wrap_val) {
      value_ -= wrap_val;
      overflow++;
    }
    return overflow;
  }
  int_type subtract(int_type v) 
  {
    //take the mod here and assign it....
    int_type remainder = v % wrap_val;
    int_type underflow = v / wrap_val;
    
//     std::cout << "wi :" << value_ << "|"
//               << v << "|"
//               << underflow << "|" 
//               << remainder << "|" 
//               << wrap_val  << std::endl;
    if (remainder > value_) {
      underflow++;
      //      value_ = remainder - value_;
      value_ = wrap_val - (remainder-value_);
    }
    else {
      value_ -= remainder;
      //value_ = wrap_val-(remainder-value_);
      //      value_ = wrap_val -(value_-remainder);
    }
//     std::cout << "wi final uf: " << underflow 
//               << " value: "  << value_ << std::endl;
    return underflow;
  }
            
private:
  int_type value_;
                  
};


//! A wrapping integer used to wrap around at the top
/*! Bad name, quick impl to fix a bug -- fix later!!
 *  This allows the wrap to restart at a value other than 0.
 *  Currently this only works if wrap_min == 1 
 */
template<typename int_type_, int_type_ wrap_min, int_type_ wrap_max>
class wrapping_int2 {
public:
  typedef int_type_ int_type;
  static unsigned long wrap_value() {return wrap_max;}
  static unsigned long min_value()  {return wrap_min;}
  /*! If initializing value is out of range of [wrap_min, wrap_max],
   * value will be initialized to closest of min or max */
  wrapping_int2(int_type v) : value_(v) {
    if(value_ < wrap_min)
    {
      value_ = wrap_min;
    }
    if(value_ > wrap_max)
    {
      value_ = wrap_max;
    }
  }
  //! Explicit converion method
  int_type as_int()   const   {return value_;}
  operator int_type() const {return value_;}
  //!Add, return number of wraps performed 
  int_type add(int_type v) 
  {
    int_type remainder = static_cast<int_type>(v % (wrap_max - wrap_min + 1));
    int_type overflow = static_cast<int_type>(v / (wrap_max - wrap_min + 1));
    value_ = static_cast<int_type>(value_ + remainder);
    if ((value_) > wrap_max) 
    {
      overflow++;
      value_ -= (wrap_max - wrap_min + 1);
    }
    return overflow;
  }
  //! Subtract will return '-d' if wrapping took place ('d' is the number of wraps)
  int_type subtract(int_type v) 
  {
    int_type remainder = static_cast<int_type>(v % (wrap_max - wrap_min + 1));
    int_type underflow = static_cast<int_type>(-(v / (wrap_max - wrap_min + 1)));
    value_ = static_cast<int_type>(value_ - remainder);
    if ((value_) < wrap_min) 
    {
      underflow--;
      value_ += (wrap_max - wrap_min + 1);
    }
    return underflow;
  }
            
private:
  int_type value_;
                  
};



} } //namespace date_time



#endif

