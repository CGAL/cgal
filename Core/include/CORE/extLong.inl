/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: extLong.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

//  private comparison function
CORE_INLINE int extLong::compare(const extLong& x) const {
  if (isNaN() || x.isNaN()) {
    core_error("Two extLong NaN's cannot be compared!", 
                __FILE__, __LINE__, false);
  }

  if (val > x.val)
    return 1;
  else if (val == x.val)
    return 0;
  else
    return - 1;
}

//  constructors
//CORE_INLINE extLong::extLong() : val(0), flag(0) { }

CORE_INLINE extLong::extLong(int i) : val(i), flag(0) {
  if (val == LONG_MAX)
    flag = 1;
  else if (val <= LONG_MIN + 1)
    flag = - 1;
}

CORE_INLINE extLong::extLong(unsigned int ui) : val(ui), flag(0) {
  if (val >= LONG_MAX)
    flag = 1;
}

CORE_INLINE extLong::extLong(long l) : val(l), flag(0) {
  if (val >= LONG_MAX)
    flag = 1;
  else if (val <= LONG_MIN + 1)
    flag = - 1;
}

CORE_INLINE extLong::extLong(unsigned long u) {
  if (u >= (unsigned long)(LONG_MAX)) {
    val  = LONG_MAX;
    flag = 1;
  } else {
    val = static_cast<long>(u);
    flag = 0;
  }
}

// isNaN defaults to false
CORE_INLINE extLong::extLong(bool isNaN) : val(0), flag(0) {
  if (isNaN) {
    val = LONG_MIN;
    flag = 2;
  }
}

// comparison operators
CORE_INLINE bool operator== (const extLong& x, const extLong& y) {
  return x.compare(y) == 0;
}

CORE_INLINE bool operator!= (const extLong& x, const extLong& y) {
  return x.compare(y) != 0;
}

CORE_INLINE bool operator< (const extLong& x, const extLong& y) {
  return x.compare(y) < 0;
}

CORE_INLINE bool operator<= (const extLong& x, const extLong& y) {
  return x.compare(y) <= 0;
}

CORE_INLINE bool operator> (const extLong& x, const extLong& y) {
  return x.compare(y) > 0;
}

CORE_INLINE bool operator>= (const extLong& x, const extLong& y) {
  return x.compare(y) >= 0;
}

//  arithmetic operators
CORE_INLINE extLong operator+ (const extLong& x, const extLong& y) {
  extLong r = x; 
  r += y; 
  return r;
}

CORE_INLINE extLong operator- (const extLong& x, const extLong& y) {
  extLong r = x;
  r -= y;
  return r;
}

CORE_INLINE extLong operator* (const extLong& x, const extLong& y) {
  extLong r = x; 
  r *= y; 
  return r;
}

CORE_INLINE extLong operator/ (const extLong& x, const extLong& y) {
  extLong r = x;
  r /= y;
  return r;
}

CORE_INLINE extLong& extLong::operator++ () {
  *this += 1;
  return *this;
}
CORE_INLINE extLong extLong::operator++ (int) {
  extLong t = *this;
  *this += 1;
  return t;
}

CORE_INLINE extLong& extLong::operator-- () {
  *this -= 1;
  return *this;
}

CORE_INLINE extLong extLong::operator-- (int) {
  extLong t = *this;
  *this -= 1;
  return t;
} 

//  conversion to long
CORE_INLINE long extLong::toLong() const {
  return val;
}

// builtin functions
CORE_INLINE long extLong::asLong() const {
  return val; 
}

CORE_INLINE bool extLong::isInfty() const {
  return (flag == 1); 
}

CORE_INLINE bool extLong::isTiny() const {
  return (flag == -1); 
}

CORE_INLINE bool extLong::isNaN() const { 
  return (flag == 2); 
}

#endif

