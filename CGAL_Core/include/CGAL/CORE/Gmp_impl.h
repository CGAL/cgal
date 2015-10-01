/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * file: GmpIO.cpp
 * 		Adapted from multi-files under /cxx in GMP's source distribution
 *
 * Zilin Du, 2003
 *
 * $URL$
 * $Id$
 ***************************************************************************/

/* Auxiliary functions for C++-style input of GMP types. 

Copyright 2001 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/CORE/Gmp.h>
#include <cctype>
#include <iostream>
#include <string>
#include <cstdio>

using namespace std;

namespace CORE { 

CGAL_INLINE_FUNCTION
int
__gmp_istream_set_base (istream &i, char &c, bool &zero, bool &showbase)
{
  int base;

  zero = showbase = false;
  switch (i.flags() & ios::basefield)
    {
    case ios::dec:
      base = 10;
      break;
    case ios::hex:
      base = 16;
      break;
    case ios::oct:
      base = 8;
      break;
    default:
      showbase = true; // look for initial "0" or "0x" or "0X"
      if (c == '0')
	{
	  if (! i.get(c))
	    c = 0; // reset or we might loop indefinitely

	  if (c == 'x' || c == 'X')
	    {
	      base = 16;
	      i.get(c);
	    }
	  else
	    {
	      base = 8;
	      zero = true; // if no other digit is read, the "0" counts
	    }
	}
      else
	base = 10;
      break;
    }

  return base;
}

CGAL_INLINE_FUNCTION
void
__gmp_istream_set_digits (string &s, istream &i, char &c, bool &ok, int base)
{
  switch (base)
    {
    case 10:
      while (isdigit(c))
	{
	  ok = true; // at least a valid digit was read
	  s += c;
	  if (! i.get(c))
	    break;
	}
      break;
    case 8:
      while (isdigit(c) && c != '8' && c != '9')
	{
	  ok = true; // at least a valid digit was read
	  s += c;
	  if (! i.get(c))
	    break;
	}
      break;
    case 16:
      while (isxdigit(c))
	{
	  ok = true; // at least a valid digit was read
	  s += c;
	  if (! i.get(c))
	    break;
	}
      break;
    }
}

CGAL_INLINE_FUNCTION
istream &
//operator>> (istream &i, mpz_ptr z)
io_read (istream &i, mpz_ptr z)
{
  int base;
  char c = 0;
  string s;
  bool ok = false, zero, showbase;

  i.get(c); // start reading

  if (i.flags() & ios::skipws) // skip initial whitespace
    while (isspace(c) && i.get(c))
      ;

  if (c == '-' || c == '+') // sign
    {
      if (c == '-') // mpz_set_str doesn't accept '+'
	s = "-";
      i.get(c);
    }

  while (isspace(c) && i.get(c)) // skip whitespace
    ;

  base = __gmp_istream_set_base(i, c, zero, showbase); // select the base
  __gmp_istream_set_digits(s, i, c, ok, base);         // read the number

  if (i.good()) // last character read was non-numeric
    i.putback(c);
  else if (i.eof() && (ok || zero)) // stopped just before eof
    i.clear();

  if (ok)
    mpz_set_str(z, s.c_str(), base); // extract the number
  else if (zero)
    mpz_set_ui(z, 0);
  else
    i.setstate(ios::failbit); // read failed

  return i;
}

CGAL_INLINE_FUNCTION
istream &
//operator>> (istream &i, mpq_ptr q)
io_read (istream &i, mpq_ptr q)
{
  int base;
  char c = 0;
  string s;
  bool ok = false, zero, showbase;

  i.get(c); // start reading

  if (i.flags() & ios::skipws) // skip initial whitespace
    while (isspace(c) && i.get(c))
      ;

  if (c == '-' || c == '+') // sign
    {
      if (c == '-')
	s = "-";
      i.get(c);
    }

  while (isspace(c) && i.get(c)) // skip whitespace
    ;

  base = __gmp_istream_set_base(i, c, zero, showbase); // select the base
  __gmp_istream_set_digits(s, i, c, ok, base);         // read the numerator

  if (! ok && zero) // the only digit read was "0"
    {
      base = 10;
      s += '0';
      ok = true;
    }

  if (i.flags() & ios::skipws)
    while (isspace(c) && i.get(c)) // skip whitespace
      ;

  if (c == '/') // there's a denominator
    {
      bool zero2 = false;
      int base2 = base;

      s += '/';
      ok = false; // denominator is mandatory
      i.get(c);

      while (isspace(c) && i.get(c)) // skip whitespace
	;

      if (showbase) // check base of denominator
	base2 = __gmp_istream_set_base(i, c, zero2, showbase);

      if (base2 == base || base2 == 10) // read the denominator
	__gmp_istream_set_digits(s, i, c, ok, base);

      if (! ok && zero2) // the only digit read was "0"
	{                // denominator is 0, but that's your business
	  s += '0';
	  ok = true;
	}
    }

  if (i.good()) // last character read was non-numeric
    i.putback(c);
  else if (i.eof() && ok) // stopped just before eof
    i.clear();

  if (ok)
    mpq_set_str(q, s.c_str(), base); // extract the number
  else
    i.setstate(ios::failbit); // read failed

  return i;
}

CGAL_INLINE_FUNCTION
ostream&
//operator<< (ostream &o, mpz_srcptr z)
io_write (ostream &o, mpz_srcptr z)
{ 
  char *str = new char [mpz_sizeinbase(z,10) + 2];
  str = mpz_get_str(str, 10, z);
  o << str ;
  delete[] str;
  return o; 
}

CGAL_INLINE_FUNCTION
ostream&
//operator<< (ostream &o, mpq_srcptr q)
io_write (ostream &o, mpq_srcptr q)
{ 
  // size according to GMP documentation
  char *str = new char [mpz_sizeinbase(mpq_numref(q), 10) +
                        mpz_sizeinbase (mpq_denref(q), 10) + 3];
  str = mpq_get_str(str, 10, q);
  o << str ;
  delete[] str;
  return o;
}

} //namespace CORE
