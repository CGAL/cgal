// Copyright (c) 2000  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int getsubopt(char ** optionsp,  const char * tokens[], char ** valuep)
{
  register char *s = *optionsp, *p;
  register int i;
  size_t optlen;

  *valuep = NULL;
  if (*s == '\0')
    return (-1);
  p = strchr(s, ',');		/* find next option */
  if (p == NULL) {
    p = s + strlen(s);
  } else {
    *p++ = '\0';		/* mark end and point to next */
  }
  *optionsp = p;		/* point to next option */
  p = strchr(s, '=');		/* find value */
  if (p == NULL) {
    optlen =  strlen(s);
    *valuep = NULL;
  } else {
    optlen = p - s;
    *valuep = ++p;
  }
  for (i = 0; tokens[i] != NULL; i++) {
    if ((optlen == strlen(tokens[i])) &&
        (strncmp(s, tokens[i], optlen) == 0))
      return (i);
  }
  /* no match, point value at option and return error */
  *valuep = s;
  return (-1);
}
