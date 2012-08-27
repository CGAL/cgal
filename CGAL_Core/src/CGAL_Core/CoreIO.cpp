/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: CoreIO.cpp
 *
 * Written by 
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#include <CGAL/CORE/BigFloatRep.h>
#include <CGAL/CORE/BigFloat.h>

namespace CORE { 

void core_io_error_handler(const char *f, const char *m) {
  std::cout << "\n error_handler";
  std::cout << "::" << f << "::" << m << "\n";
  std::cout.flush();
  std::abort();
}

void core_io_memory_handler(char *t, const char *f, const char *m) {
  if (t == NULL) {
    std::cout << "\n memory_handler";
    std::cout << "::" << f << "::" << m;
    std::cout << "memory exhausted\n";
    std::cout.flush();
    std::abort();
  }
}

// s has size old_size and will be resized to new_size.
void allocate (char * &s, int old_size, int new_size) {
  if (old_size > new_size)
    old_size = new_size;

  if (s == NULL)
    old_size = 0;

  char *t = new char[new_size];
  core_io_memory_handler(t, "CoreIO", "allocate::out of memory error");

  int i;
  for (i = 0; i < old_size; i++)
    t[i] = s[i];

  delete[] s;
  s = t;
}

// appends c to s at position pos.
// sz is the size of s
void append_char (char * &s, int & sz, int pos, char c) {
  if (pos > sz)
    core_io_error_handler("CoreIO", "append_char::invalid argument");

  if (pos == sz) {
    allocate(s, sz, 2*sz);
    sz *= 2;
  }

  s[pos] = c;
}

// skip blanks, tabs, line breaks and comment lines
int skip_comment_line (std::istream & in) {
  int c;

  do {
    c = in.get();
    while ( c == '#' ) {
      do {
        c = in.get();
      } while ( c != '\n' );
      c = in.get();
    }
  } while (c == ' ' || c == '\t' || c == '\n');

  if (c == EOF)
    core_io_error_handler("CoreIO::read_from_file()","unexpected end of file.");

  in.putback(c);
  return c;
}

// skips '\\' followed by '\n'
int skip_backslash_new_line (std::istream & in) {
  int c = in.get();

  while (c == '\\') {
    c = in.get();

    if (c == '\n')
      c = in.get();
    else
      core_io_error_handler("CoreIO::operator>>", "\\ must be immediately followed by new line.");
  }

  return c;
}

void read_string(std::istream& in, char* &buffer, int sz) {
  int c, pos=0;
  skip_comment_line(in);

  while ( (c = in.get()) != EOF ) {
    if ( c == ' ' || c == '\t' || c == '\n' || c == '#')
      break;
    else
      append_char(buffer, sz, pos++, c);
  }
  append_char(buffer, sz, pos, '\0');
}

void read_base_number(std::istream& in, BigInt& m, long length, long maxBits) {
  char *buffer;
  int size, offset;
  int base;
  bool is_negate;

  int c, pos = 0;
  skip_comment_line(in);

  // read sign
  c = in.get();
  if (c == '-') {
    is_negate = true;
    c = in.get();
  } else
    is_negate = false;

  // read base and compute digits
  if (c == '0') {
    c = in.get();
    if (c == 'b') {
      base = 2;
      size = (maxBits == 0 || maxBits > length) ? length : maxBits;
      offset = length - size;
    } else if (c == 'x') {
      base = 16;
      size = (maxBits == 0) ? length : (maxBits+3) >> 2;
      size = (size > length) ? length : size;
      offset = (length - size) << 2;
    } else {
      base = 8;
      size = (maxBits == 0) ? length : (maxBits+2) / 3;
      size = (size > length) ? length : size;
      offset = (length - size) * 3;
      in.putback(c);
    }
  } else {
    base = 10;
    size = (maxBits == 0) ? length : (int)std::ceil(maxBits*std::log(2.0)/std::log(10.0));
    size = (size > length) ? length : size;
    offset = length - size;
    in.putback(c);
  }

  buffer = new char[size+2];
  // read digits
  for (int i=0; (i<size)&&((c=skip_backslash_new_line(in)) != EOF ); i++) {
    if (c != ' ' && c != '\t' && c != '\n')
      append_char(buffer, size, pos++, c);
  }
  if (base == 10) {
    for(int j=0; j<offset; j++)
      append_char(buffer, size, pos++, '0');
  }
  append_char(buffer, size, pos, '\0');

  // convert string to bigint.
  if (m.set_str(buffer, base) < 0)
    core_io_error_handler("CoreIO::read_from_file()","bad big number format.");
  delete[] buffer;

  // shift left if neccessary
  if (offset > 0 && base != 10) {
    m <<= offset;
  }

  if (is_negate)
    negate(m);
}


void write_base_number(std::ostream& out, char* buffer, int length, int base, int charsPerLine) {
  // write big number in a format that gmp's mpz_set_str() can
  // automatically recognize with argument base = 0.
  if (base == 2)
    out << "0b";
  else if (base == 16)
    out << "0x";
  else if (base == 8)
    out << '0';

  // write big number in charsPerLine.
  char* start, *end, c;
  for (int i=0; i<length; i += charsPerLine) {
    start = buffer + i;
    if (i + charsPerLine >= length)
      out << start;
    else {
      end = start + charsPerLine;
      c = *end;
      *end = '\0';

      out << start << "\\\n";
      *end = c;
    }
  }
}

void readFromFile(BigInt& z, std::istream& in, long maxLength) {
  char *buffer;
  long length;

  // check type name whether it is Integer or not.
  buffer = new char[8];
  read_string(in, buffer, sizeof(buffer));
  if ( std::strcmp(buffer, "Integer") != 0)
    core_io_error_handler("BigInt::read_from_file()","type name expected.");
  delete[] buffer;

  // read the bit length field.
  buffer = new char[100];
  read_string(in, buffer, sizeof(buffer));
  length = std::atol(buffer);
  delete[] buffer;

  // read bigint
  read_base_number(in, z, length, maxLength);
}

void writeToFile(const BigInt& z, std::ostream& out, int base, int charsPerLine) {
  BigInt c = abs(z);

  // get the absoulte value string
  char* buffer = new char[mpz_sizeinbase(c.get_mp(), base) + 2];
  mpz_get_str(buffer, base, c.get_mp());
  int length = std::strlen(buffer);

  // write type name of big number and length
  //out << "# This is an experimental big number format.\n";
  out << "Integer " << length << "\n";

  // if bigint is negative, then write an sign '-'.
  if ( sign(z) < 0  )
    out << '-';

  write_base_number(out, buffer, length, base, charsPerLine);
  out << "\n";
  delete[] buffer;
}

void readFromFile(BigFloat& bf, std::istream& in, long maxLength) {
  char *buffer;
  long length;
  long exponent;
  BigInt mantissa;

  // check type name whether it is Float
  buffer = new char[6];
  read_string(in, buffer, sizeof(buffer));
  if (std::strcmp(buffer, "Float") != 0)
    core_io_error_handler("BigFloat::read_from_file()", "type name expected");
  delete[] buffer;

  // read base (default is 16384)
  buffer = new char[8];
  read_string(in, buffer, sizeof(buffer));
  if (std::strcmp(buffer, "(16384)") != 0)
    core_io_error_handler("BigFloat::read_from_file()", "base expected");
  delete[] buffer;

  // read the bit length field.
  buffer = new char[100];
  read_string(in, buffer, sizeof(buffer));
  length = std::atol(buffer);
  delete[] buffer;

  // read exponent
  buffer = new char[100];
  read_string(in, buffer, sizeof(buffer));
  exponent = std::atol(buffer);
  delete[] buffer;

  // read mantissa
  read_base_number(in, mantissa, length, maxLength);

  // construct BigFloat
  bf = BigFloat(mantissa, 0, exponent);
}

void writeToFile(const BigFloat& bf, std::ostream& out, int base, int charsPerLine) {
  BigInt c(CORE::abs(bf.m()));

  // get the absoulte value string
  char* buffer = new char[mpz_sizeinbase(c.get_mp(), base) + 2];
  mpz_get_str(buffer, base, c.get_mp());
  int length = std::strlen(buffer);


  // write type name, base, length
  //out << "# This is an experimental Big Float format." << std::endl;
  out << "Float (16384) " << length << std::endl;
  // write exponent
  out << bf.exp() << std::endl;

  // write mantissa
  if ( CORE::sign(bf.m()) < 0 )
    out << '-';

  write_base_number(out, buffer, length, base, charsPerLine);
  out << '\n';
  delete[] buffer;
}

/* Underconstruction ----
void BigFloat::read_from_file2(std::istream& in, long maxLength) {
  long length = 1024;
  char *buffer;
  
  // check type name whether it is Float
  buffer = new char[7];
  BigInt::read_string(in, buffer, sizeof(buffer));
  if (strcmp(buffer, "NFloat") != 0)
    core_io_error_handler("BigFloat::read_from_file2()", "type name expected");
  delete[] buffer;
 
  // read base (default is 16) 
  buffer = new char[5];
  BigInt::read_string(in, buffer, sizeof(buffer));
  if (strcmp(buffer, "(16)") != 0)
    core_io_error_handler("BigFloat::read_from_file2()", "base expected");
  delete[] buffer;
  
  // read length field
  buffer = new char[100];
  BigInt::read_string(in, buffer, sizeof(buffer));
  
  // get the length field if it is not null.
  if (buffer[0] != '\0') {
    length = atol(buffer);
    if (maxLength > 0 && length >= maxLength)
      length = maxLength;
  }
  delete[] buffer;
 
  // read exponent
  buffer = new char[100];
  BigInt::read_string(in, buffer, sizeof(buffer));
  long exp16 = atol(buffer);
  delete[] buffer;
 
  // read mantissa
  buffer = new char[length+2];
  //BigInt::read_base_number(in, buffer, length);
 
  BigInt m16(buffer);
  delete[] buffer;
  
  // convert to base CHUNK_BIT
  exp16 = exp16 - length + 1; 
  if ( m16.is_negative() )
    exp16 ++;
 
  long tmp_exp = exp16 * 4;
  long q = tmp_exp / CHUNK_BIT;
  long r = tmp_exp % CHUNK_BIT;
  if ( r < 0 ) {
    r += CHUNK_BIT;
    q --;
  }
  
  BigInt mantissa = m16 << r;
  long exponent = q;
 
  // construct BigFloat
  if (--rep->refCount == 0)
    delete rep;
  
  rep = new BigFloatRep(mantissa, 0, exponent);
  rep->refCount++;
  
}
 
// write normal float
// now it assumed to write in hex base, i.e. B=2^4=16
// (note: our default base B=2^(CHUNK_BIT)=2^14=16384
void BigFloat::write_to_file2(std::ostream& out, int base, int charsPerLine) {
  // convert to base 16.
  long new_base = 4; // 2^4 = 16
  long tmp_exp = (rep->exp) * CHUNK_BIT;
  long q = tmp_exp / new_base;
  long r = tmp_exp % new_base;
  std::cout << "CORE_DEBUG: q=" << q << ", r=" << r << std::endl;
  if ( r < 0 ) {
    r += new_base;
    q--;
  }
  std::cout << "CORE_DEBUG: q=" << q << ", r=" << r << std::endl;
  
  BigInt m16 = (rep->m) << r;
 
  int size = mpz_sizeinbase(m16.I, base) + 2;
  std::cout << "size=" << size << std::endl;
  char* buffer = new char[size];
 
  int length = bigint_to_string(m16, buffer, base);
  std::cout << "length=" << length << std::endl;
 
  long exp16 = q + length - 1; 
  if ( m16.is_negative() )
    exp16 --;
 
  // write type name, base, length
  out << "# This is an experimental Big Float format." << std::endl;
  out << "NFloat (16) " << length << std::endl;
  
  // write exponent
  out << exp16 << std::endl;
  
  // write mantissa
  if ( m16.is_negative() ) {
    out << '-';
    buffer ++;
  }
  
  BigInt::write_base_number(out, buffer, length, base, charsPerLine);
  out << '\n';
  delete[] buffer;
}
*/

} //namespace CORE

