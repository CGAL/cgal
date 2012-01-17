// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Radu Ursu



#include <CGAL/Quotient.h>

template <class NT>
CGAL::Quotient<NT>
double_to_quotient(double x)
{
    NT num = 0;
    NT den = 1;

    if (x != 0.0)
    { int neg = (x < 0);
      if (neg) x = -x;

      const unsigned shift = 15;   // a safe shift per step
      const unsigned int shift_pow = 32768; // = 2^shift
      const double width = 32768;  // = 2^shift
      const int maxiter = 20;      // ought not be necessary, but just in case,
                                   // max 300 bits of precision
      int expt;
      double mantissa = std::frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;

      while (mantissa != 0.0 && k++ < maxiter)

      { mantissa *= width; // shift double mantissa
        mantissa = std::modf(mantissa, &intpart);
        num *= shift_pow;
        num += (long)intpart;
        exponent -= shift;
      }
      int expsign = (exponent>0 ? +1 : (exponent<0 ? -1 : 0));
      exponent *= expsign;
      NT twopot(2);
      NT exppot(1);
      while (exponent!=0) {
        if (exponent & 1)
          exppot *= twopot;
        exponent >>= 1;
        twopot *= twopot;
      }

      if (expsign > 0)
        num *= exppot;
      else if (expsign < 0)
        den *= exppot;
      if (neg)
        num = -num;
    }
    CGAL::Quotient<NT> q(num,den);
    q.normalize();
    return q;
}
/* XPM */
static const char *intersection_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c opaque",
". c #808080",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXX  XXXXX    XXXX   XXXXXXX",
"XXXXXX XXXXXXX . XXXXXXXX XXXXXX",
"XXXXX XXXXXXX ... XXXXXXXX XXXXX",
"XXXX XXXXXXX ..... XXXXXXXX XXXX",
"XXXX XXXXXXX ..... XXXXXXXX XXXX",
"XXX XXXXXXX ....... XXXXXXXX XXX",
"XXX XXXXXXX ....... XXXXXXXX XXX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XX XXXXXXX ......... XXXXXXXX XX",
"XXX XXXXXXX ....... XXXXXXXX XXX",
"XXX XXXXXXX ....... XXXXXXXX XXX",
"XXXX XXXXXXX ..... XXXXXXXX XXXX",
"XXXX XXXXXXX ..... XXXXXXXX XXXX",
"XXXXX XXXXXXX ... XXXXXXXX XXXXX",
"XXXXXX XXXXXXX . XXXXXXXX XXXXXX",
"XXXXXXX  XXXXX    XXXX   XXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static const char *union_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c opaque",
". c #808080",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXX  .....    ....   XXXXXXX",
"XXXXXX ....... . ........ XXXXXX",
"XXXXX ....... ... ........ XXXXX",
"XXXX ....... ..... ........ XXXX",
"XXXX ....... ..... ........ XXXX",
"XXX ....... ....... ........ XXX",
"XXX ....... ....... ........ XXX",
"XX ....... ......... ........ XX",
"XX ....... ......... ........ XX",
"XX ....... ......... ........ XX",
"XX ....... ......... ........ XX",
"XX ....... ......... ........ XX",
"XX ....... ......... ........ XX",
"XXX ....... ....... ........ XXX",
"XXX ....... ....... ........ XXX",
"XXXX ....... ..... ........ XXXX",
"XXXX ....... ..... ........ XXXX",
"XXXXX ....... ... ........ XXXXX",
"XXXXXX ....... . ........ XXXXXX",
"XXXXXXX  .....    ....   XXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static const char *difference_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c opaque",
". c #808080",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXX  .....    XXXX   XXXXXXX",
"XXXXXX ....... X XXXXXXXX XXXXXX",
"XXXXX ....... XXX XXXXXXXX XXXXX",
"XXXX ....... XXXXX XXXXXXXX XXXX",
"XXXX ....... XXXXX XXXXXXXX XXXX",
"XXX ....... XXXXXXX XXXXXXXX XXX",
"XXX ....... XXXXXXX XXXXXXXX XXX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XX ....... XXXXXXXXX XXXXXXXX XX",
"XXX ....... XXXXXXX XXXXXXXX XXX",
"XXX ....... XXXXXXX XXXXXXXX XXX",
"XXXX ....... XXXXX XXXXXXXX XXXX",
"XXXX ....... XXXXX XXXXXXXX XXXX",
"XXXXX ....... XXX XXXXXXXX XXXXX",
"XXXXXX ....... X XXXXXXXX XXXXXX",
"XXXXXXX  .....    XXXX   XXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static const char *symmetric_difference_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c opaque",
". c #808080",
"X c None",
/* pixels */
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXX  .....    ....   XXXXXXX",
"XXXXXX ....... X ........ XXXXXX",
"XXXXX ....... XXX ........ XXXXX",
"XXXX ....... XXXXX ........ XXXX",
"XXXX ....... XXXXX ........ XXXX",
"XXX ....... XXXXXXX ........ XXX",
"XXX ....... XXXXXXX ........ XXX",
"XX ....... XXXXXXXXX ........ XX",
"XX ....... XXXXXXXXX ........ XX",
"XX ....... XXXXXXXXX ........ XX",
"XX ....... XXXXXXXXX ........ XX",
"XX ....... XXXXXXXXX ........ XX",
"XX ....... XXXXXXXXX ........ XX",
"XXX ....... XXXXXXX ........ XXX",
"XXX ....... XXXXXXX ........ XXX",
"XXXX ....... XXXXX ........ XXXX",
"XXXX ....... XXXXX ........ XXXX",
"XXXXX ....... XXX ........ XXXXX",
"XXXXXX ....... X ........ XXXXXX",
"XXXXXXX  .....    ....   XXXXXXX",
"XXXXXXXXX     XXXX    XXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
};

/* XPM */
static const char *complement_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 2 1",
"  c opaque",
". c #808080",
/* pixels */
"................................",
"................................",
"................................",
"................................",
"................................",
".....  ..................  .....",
".....                      .....",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
"......                    ......",
".....                      .....",
".....  ..................  .....",
"................................",
"................................",
"................................",
"................................",
"................................"
};

/* XPM */
static const char *interior_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 2 1",
"  c opaque",
". c #808080",
/* pixels */
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"       ..................       ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                "
};

/* XPM */
static const char *closure_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 3 1",
"  c opaque",
". c #808080",
"X c white",
/* pixels */
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"     XXXXXXXXXXXXXXXXXXXXXX     ",
"     XXXXXXXXXXXXXXXXXXXXXX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XX..................XX     ",
"     XXXXXXXXXXXXXXXXXXXXXX     ",
"     XXXXXXXXXXXXXXXXXXXXXX     ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                "
};

/* XPM */
static const char *boundary_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 2 1",
"  c opaque",
". c white",
/* pixels */
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"     ......................     ",
"     ......................     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ..                  ..     ",
"     ......................     ",
"     ......................     ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                "
};

/* XPM */
static const char *regularization_xpm[] = {
/* columns rows colors chars-per-pixel */
"32 32 4 1",
"  c opaque",
". c yellow",
"X c #808080",
"o c white",
/* pixels */
"                                ",
"                   ..           ",
"                    ..      ..  ",
"                     ..    ..   ",
"                     ..   ...   ",
"                      ......    ",
"       ooooooooooooooo.....     ",
"       ooooooooooooooo.....     ",
"       ooXXXXXXXXXXXXX.....     ",
"       ooXXXXXXXXXXXX.......    ",
"       ooXXXXXXXXXXX..Xoo ...   ",
"       ooXXXXXXXXXX..XXoo  ...  ",
"       ooXXXXXXXXXXXXXXoo   .   ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       ooXXXXXXXXXXXXXXoo       ",
"       oooooooooooooooooo       ",
"       oooooooooooooooooo       ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                ",
"                                "
};
