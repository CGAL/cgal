/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigInt.cpp
 * Synopsis:
 *      Provides the basic big integer capabilities of Core Library
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/
// LiDIA - a library for computational number theory
//   Copyright (c) 1994, 1995, 2000 by the LiDIA Group
//
// File        : bigint.c (implementation of the gmp interface) 
// Author      : Thomas Papanikolaou (TP)
// Last change : TP, Feb  7 1995, initial version
//               MM, Sep 17 1996, 


#include "BigInt.h"

CORE_BEGIN_NAMESPACE

#ifndef CORE_ENABLE_INLINES
#include "BigInt.inl"
#endif

const bigint& bigint::getAnonymous() {
  static bigint anonymous(0);
  return anonymous;
}

void bigint::getKaryExpo(bigint& m, int& e, unsigned long k) const {
  bigint q, r;
  m = *this; e = 0;
  mpz_tdiv_qr_ui(q.I, r.I, m.I, k);
  while (mpz_sgn(r.I) == 0) {
    e ++; m = q; 
    mpz_tdiv_qr_ui(q.I, r.I, m.I, k);
  }
}

/**
 **  The maximal number of characters of which the
 **  operating system is able to put in one line.
 **/

int bigint::chars_per_line = 80;

void remainder(bigint & c, const bigint & a, const bigint & b)
{ 
  // can be optimized by using static tmp
  //
  // Problem: see first MM-entry in Last Changes, mpz_mod

  if (a < 0) // rest must be negative
   {
     if (&b == &c)  // b will be overwritten
      {
        mpz_t tmp;   // set tmp = |b|
        mpz_init_set(tmp, b.I);
        if (b < 0) // b is negative
          mpz_neg(tmp, tmp);

        mpz_mod(c.I, a.I, b.I); // compute non-negative rest

        if (c != 0) // rest is not equal to zero
           mpz_sub(c.I, c.I, tmp); // create negative rest
        
        mpz_clear(tmp);
      }
     else // b won't be overwritten
      {
        mpz_mod(c.I, a.I, b.I); // compute non-negative rest

        if (c != 0) // rest is not equal to zero
         {
           if (b < 0)  // create negative rest
              mpz_add(c.I, c.I, b.I);
           else
              mpz_sub(c.I, c.I, b.I);
         }
      }
   }
  else
     mpz_mod(c.I, a.I, b.I);
}


void power(bigint & c, const bigint & a, const bigint & b)
{
  bigint exponent, multiplier;
  if (mpz_cmpabs_ui(a.I, 1) == 0) // a = +-1 -> a^(|b|)
    {
      if (b.is_odd())
        c = a;
      else
        c.assign_one();
    }
  else if (b.is_negative()) // b < 0, a != +-1 -> 0
    c.assign_zero();
  else if (b.is_zero()) // b == 0 -> 1
    c.assign_one();
  else if (b.is_one())  // b == 1 -> a
    c.assign(a);
  else
  {
    exponent.assign(b);

    multiplier = a;
    c = a;

    // Note: The bitlength is at least 2.
    // Zilin: change lidia_std::size_t to int   
    int length = exponent.bit_length()-2;
    bigint tmp;

    shift_left(tmp,1,length);

    while (tmp.is_gt_zero())
     {
       square(c, c);
       if (!((exponent&tmp).is_zero()))
         multiply(c, c, multiplier);
       tmp.divide_by_2();
     }
  }
}

void power(bigint & c, const bigint & a, long i)
{
  if (mpz_cmpabs_ui(a.I, 1) == 0) // a = +-1 -> a^(|b|)
    {
      if (i&1)
        c = a;
      else
        c.assign_one();
    }
  else if (i < 0)  // i < 0, a != +-1 -> 0
    c.assign_zero();
  else if (i == 0) // i == 0 -> 1
    c.assign_one();
  else 
    mpz_pow_ui(c.I, a.I, i); 
} 

// Note: Remove below 4 functions since it seems we never use it in Core Library
//      Zilin Du, June 14, 2001

/**
** using fread/fwrite
**/

/*  
void bigint::read_from_file(FILE * fp)
{ mpz_inp_raw(&I, fp); }

void bigint::write_to_file(FILE * fp)
{ mpz_out_raw(fp, &I); }
*/

/**
** using fscanf/fprintf
**/

/*
void bigint::scan_from_file(FILE * fp)
{ mpz_inp_str(&I, fp, 10); }

void bigint::print_to_file(FILE * fp)
{ mpz_out_str(fp, 10, &I); }
*/

// Note: Add 2 function to implement read/write BigInt in our own format
//       Zilin Du, June 14, 2001

// skip blanks, tabs, line breaks and comment lines
int bigint::skip_comment_line (std::istream & in)
{
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
    lidia_error_handler("bigint::read_from_file()","unexpected end of file.");
  
  in.putback(c);
  return c;
}

void bigint::read_string(std::istream& in, char* &buffer, int sz) {  
  int c, pos=0;  
  skip_comment_line(in);  
  
  while ( (c = in.get()) != EOF ) {
    if ( c == ' ' || c == '\t' || c == '\n' || c == '#')
      break;
    else
      bigint::append_char(buffer, sz, pos++, c);
  }
  bigint::append_char(buffer, sz, pos, '\0');
}

void bigint::read_base_number(std::istream& in, bigint& m, long length, long maxBits) {
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
    size = (maxBits == 0) ? length : (int)::ceil(maxBits*log(2.0)/log(10.0));
    size = (size > length) ? length : size;
    offset = length - size;
    in.putback(c);
  }

  buffer = new char[size+2];  
  // read digits    
  for (int i=0; (i<size)&&((c=skip_backslash_new_line(in)) != EOF ); i++) {
    if (c != ' ' && c != '\t' && c != '\n')
      bigint::append_char(buffer, size, pos++, c); 
  }
  if (base == 10) {
    for(int j=0; j<offset; j++)
      bigint::append_char(buffer, size, pos++, '0'); 
  } 
  bigint::append_char(buffer, size, pos, '\0');

  // convert string to bigint. 
  if (string_to_bigint(buffer, m, base) < 0) 
    lidia_error_handler("bigint::read_from_file()","bad big number format.");  
  delete[] buffer;

  // shift left if neccessary
  if (offset > 0 && base != 10) {
    m <<= offset;
  }
    
  if (is_negate)
    m.negate();
}  
  
void bigint::read_from_file(std::istream& in, long maxLength)
{
  char *buffer;
  long length;
  
  // check type name whether it is Integer or not.
  buffer = new char[8];
  read_string(in, buffer, sizeof(buffer));
  if ( strcmp(buffer, "Integer") != 0)
    lidia_error_handler("bigint::read_from_file()","type name expected.");  
  delete[] buffer;
  
  // read the bit length field.
  buffer = new char[100];
  read_string(in, buffer, sizeof(buffer));
  length = atol(buffer);
  delete[] buffer;  
  
  // read bigint
  read_base_number(in, *this, length, maxLength);
}

void bigint::write_base_number(std::ostream& out, char* buffer, int length,
        int base, int charsPerLine) {
  
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

void bigint::write_to_file(std::ostream& out, int base, int charsPerLine)
{
   // Bug fixed: allocate exact space to hold result string according
   //            to GNU MP manual.
   // Zilin Du, 07/11/01
   //
   // original code:
   //   char* buffer = new char[MAX_BIGINT_LENGTH];
   bigint c = *this; 
   c.abs();
   
   char* buffer = new char[mpz_sizeinbase(c.I, base)];
   int length = bigint_to_string(c, buffer, base);
   
   // write type name of big number and length  
   out << "# This is an experimental big number format.\n";
   out << "Integer " << length << "\n";
   
   // if bigint is negative, then write an sign '-'.
   if ( is_negative() ) {
     out << '-';
   }
   
   write_base_number(out, buffer, length, base, charsPerLine);
   out << "\n";

   delete[] buffer;  
}

//
// s has size old_size and will be resized to new_size.
//

void bigint
::allocate (char * &s, int old_size, int new_size)
 {
   if (old_size > new_size)
     old_size = new_size;

   if (s == NULL)
       old_size = 0;
   
   char *t = new char[new_size];
   memory_handler(t, "bigint", "allocate::out of memory error");

   int i;
   for (i = 0; i < old_size; i++)
       t[i] = s[i];

   delete[] s;
   s = t;
 }
 

//
// appends c to s at position pos.
// sz is the size of s
//

void bigint
::append_char (char * &s, int & sz, int pos, char c)
 {
   if (pos > sz)
       lidia_error_handler("bigint", "append_char::invalid argument");

   if (pos == sz)
    {
       bigint::allocate(s, sz, 2*sz);
       sz *= 2;
    }

   s[pos] = c;
 }



//
// skips '\\' followed by '\n'
//

int bigint
::skip_backslash_new_line (std::istream & in)
 {
   int c = in.get();

   while (c == '\\')
    {
      c = in.get();
      
      if (c == '\n')
        c = in.get();
      else
        {
          lidia_error_handler("bigint::operator>>",
                              "\\ must be immediately followed by new line.");
        }
    }

   return c;
 }

//
// reads in a bigint in ASCII-format from in
//

void bigint
::scan (std::istream & in)
{
                int  sz = 10000;
                char *buffer = new char[sz];
                buffer[0] = '\0';
                int pos = 0;

                bool stop;
                int c;

                
                // skip blanks and line breaks
                do
                  {
                    c = in.get();
                  }
                while (c == ' ' || c == '\n');

                in.putback(c);


                // skip blanks, line breaks, and backslashs 
                // followed by line break
                do 
                  {
                    c = bigint::skip_backslash_new_line(in);

                    // end of stream, complain

                    if (c == EOF)
                      {
                        delete[] buffer;
                        lidia_error_handler("bigint::operator>>",
                                            "bigint expected.");
                      }
                  } while (c == ' ' || c == '\n');


                // handle sign
                if (c == '+' || c == '-')
                  {
                    buffer[pos++] = c;
                    c = bigint::skip_backslash_new_line(in);
                  }

                // require digit now
                if (!isdigit(c))
                  { 
                    delete[] buffer;
                    stop = true;
                    lidia_error_handler ("bigint::operator>>",
                                         "digit expected");
                  }
                else
                  {
                    bigint::append_char(buffer, sz, pos++, c);
                    stop = false;
                  }


                while (!stop)
                  {
                    c = in.get();

                    // store digit and continue
                    if (isdigit(c))
                        bigint::append_char(buffer, sz, pos++, c);

                    // if line break, stop
                    else if (c == '\n')
                      {
                        stop = true;
                        in.putback(c);
                      }

                    // skip "\\\n".
                    // stop if "\\c" and putback "\\c"
                    else if (c == '\\')
                      {
                        c = in.get();
                       
                        if (c != '\n')
                         {
                           stop = true;
                           in.putback(c);
                           in.putback('\\');
                         }
                      }

                    // stop if EOF
                    else if (c == EOF)
                        stop = true;

                    // stop otherwise and putback c
                    else
                      {
                        stop = true;
                        in.putback(c);
                      }
                  }

                bigint::append_char(buffer, sz, pos, '\0');
                string_to_bigint(buffer, *this);
                delete [] buffer;
                return;         
}


//
// Prints a bigint (represented by the string s)
// to out, where at most bigint::chars_per_line
// characters are printed in one line.
//

void bigint
::print (std::ostream & out, char *s) const
{
  if (bigint::chars_per_line <= 1)
    out << s;
  else
    {
      int l = strlen(s);
      char *start, *end, c;
      int  i;

      for (i=0; i < l; i += bigint::chars_per_line)
        {
          start = s + i;

          if (i + bigint::chars_per_line >= l)
            out << start;
          else
            {
              end = start + bigint::chars_per_line;

              c = *end;
              *end = '\0';
      
              out << start;
              out << '\\';
              out << '\n';
              *end = c;
            }
        }
    }
}

void newton_root(bigint & b, const bigint & a, int n)
{
  bigint c, d;

  if (n < 0)
    lidia_error_handler("bigint", "newton_root::negative n.");

  if (a.is_negative())
    lidia_error_handler("bigint", "newton_root::negative a.");

  switch (n)
  {
  case 0:
    b.assign_one();
    break;
  case 1:
    b.assign(a);
    break;
  case 2:
    sqrt(b, a);
    break;
  default:
    b.assign_one();
    shift_left(b, b, (unsigned int)((a.bit_length() + n - 1) / n));
    do{
      power(c, b, n - 1);
      div_rem(c, d, a, c);
      subtract(c, b, c);
      divide(c, c, n);
      subtract(b, b, c);
    } while (c.sign() > 0);
    power(c, b, n);
    if (c.compare(a) > 0)
      dec(b);
    if (b.compare(3) == 0)
      {
        power(c, b, n);
        if (c.compare(a) > 0)
          dec(b);
      }
    break;
  }
}

CORE_END_NAMESPACE
