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

#ifndef CGAL_MY_GETOPT_H
#define CGAL_MY_GETOPT_H

#ifdef __cplusplus
extern "C" {
#endif

extern int getopt(int argc, char * argv[], const char * opts);

extern char * optarg;
extern int opterr;
extern int optind;
extern int optopt;

#define	GETOPTDONE	(-1)

/*
 * This value is returned when an option argument is found that is not
 * in the option list
 */
#define	GETOPTHUH	'?'

#if (defined _MSC_VER)
extern int getsubopt(char ** optionsp, char * const * tokens, char ** valuep);
extern void getoptreset(void);
#endif
  
#ifdef __cplusplus
}
#endif

#endif
