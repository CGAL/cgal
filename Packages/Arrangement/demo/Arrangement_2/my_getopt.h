#ifndef __GETOPT_H__
#define __GETOPT_H__
#ident "$Revision$"
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Declarations for getopt(3C).
 *
 * Copyright 1990, Silicon Graphics, Inc.
 * All Rights Reserved.
 *
 * This is UNPUBLISHED PROPRIETARY SOURCE CODE of Silicon Graphics, Inc.;
 * the contents of this file may not be disclosed to third parties, copied or
 * duplicated in any form, in whole or in part, without the prior written
 * permission of Silicon Graphics, Inc.
 *
 * RESTRICTED RIGHTS LEGEND:
 * Use, duplication or disclosure by the Government is subject to restrictions
 * as set forth in subdivision (c)(1)(ii) of the Rights in Technical Data
 * and Computer Software clause at DFARS 252.227-7013, and/or in similar or
 * successor clauses in the FAR, DOD or NASA FAR Supplement. Unpublished -
 * rights reserved under the Copyright Laws of the United States.
 */

// #include <standards.h>

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

#endif /* !__GETOPT_H__ */
