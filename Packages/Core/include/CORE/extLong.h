/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: extLong.h
 *
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>,
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>,
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_extLong_H
#define CORE_extLong_H

#include <CORE/CoreImpl.h>
#include <CORE/CoreAux.h>

CORE_BEGIN_NAMESPACE

#ifndef LONG_MAX 
#error "haven't define LONG_MAX"
#endif

#ifndef LONG_MIN 
#error "haven't define LONG_MIN"
#endif

// LONG_MAX is assumed in this class:

/// \class extLong
/// \brief extended long integer
class extLong
{
private:  
  long val;  ///< internal representation
  int  flag; ///< flags 
             /**<  0 -- Normal;
                   1 -- Overflow (positive);
                  -1 -- Overflow (negative);
                   2 -- NaN (sign can not be determined) */
  
  /// comparison  
  int compare(const extLong &) const;
  static long add4Long(long x, long y);
  static long sub4Long(long x, long y);

public:
                 
  /// \name Constructors
  //@{
  /// default constructor
  extLong(bool isNaN = false);
  /// constructor for \c int
  extLong(int);
  /// constructor for \c unsigned int
  extLong(unsigned int);
  /// constructor for \c long
  extLong(long);
  /// constructor for \c unsigned long
  extLong(unsigned long);
  //@}
  
  /// \name Comparison operators
  //@{
  friend bool operator==(const extLong&, const extLong&);
  friend bool operator!=(const extLong&, const extLong&);
  friend bool operator< (const extLong&, const extLong&);
  friend bool operator<=(const extLong&, const extLong&);
  friend bool operator> (const extLong&, const extLong&);
  friend bool operator>=(const extLong&, const extLong&);
  //@}

  /// \name Arithmetic operators
  //@{
  friend extLong operator+(const extLong&, const extLong&);
  friend extLong operator-(const extLong&, const extLong&);
  friend extLong operator*(const extLong&, const extLong&);
  friend extLong operator/(const extLong&, const extLong&);
  //@}
  
  /// \name Arithmetic and assignment operators  
  //@{
  extLong& operator +=(const extLong&);
  extLong& operator -=(const extLong&);
  extLong& operator *=(const extLong&);
  extLong& operator /=(const extLong&);
  //@}
  
  /// \name Incremental, Decremental, Unary minus operators
  //@{
  extLong& operator++();
  extLong  operator++(int);
  extLong& operator--();
  extLong  operator--(int);
  extLong  operator-() const;
  //@}
  
  /// \name Conversion Function
  //@{
  CORE_INLINE long toLong() const;
  //@}
 
  /// \name Builtin functions  
  //@{
  CORE_INLINE long asLong() const;
  CORE_INLINE bool isInfty() const;
  CORE_INLINE bool isTiny() const;
  CORE_INLINE bool isNaN() const;
  int  sign() const;
  //@}

  /// \name I/O Stream
  ///@{
  friend std::ostream& operator <<(std::ostream&, const extLong&);
  //@}

  static const extLong& getNaNLong();
  static const extLong& getPosInfty();
  static const extLong& getNegInfty();
};

// constants (Globally)
#define CORE_NaNLong extLong::getNaNLong()
#define CORE_posInfty extLong::getPosInfty()
#define CORE_negInfty extLong::getNegInfty()

#ifdef CORE_ENABLE_INLINES
#include <CORE/extLong.inl>
#endif

CORE_END_NAMESPACE
#endif

