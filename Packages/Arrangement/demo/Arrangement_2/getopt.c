/*  Copyright (c) 1990, 1991 UNIX System Laboratories, Inc. */
/*	Copyright (c) 1988 AT&T	*/
/*	  All Rights Reserved  	*/

/*	THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF     	*/
/*	UNIX System Laboratories, Inc.                     	*/
/*	The copyright notice above does not evidence any   	*/
/*	actual or intended publication of such source code.	*/

#ifndef __LINUX__
#ident  "@(#)libc-port:gen/getopt.c     1.28"

#ifdef __STDC__
        #pragma weak getopt = _getopt
#if !_LIBC_ABI
        #pragma weak getoptreset = _getoptreset
#endif
#endif
#endif // linux

// #include "synonyms.h"
// #include <unistd.h>
#include <string.h>
#include <stdio.h>
// #include <pfmt.h>
// #include "../locale/pfmt_data.h"
#include <stdlib.h>

#include "my_getopt.h"

char * optarg = NULL;
int  opterr;
int  optind = 1;
int  optopt = 0;

int _sp = 1;
static void pr_error(const char *, const char *, int);

int
getopt(int argc, char * argv[], const char *opts)
{
	register char c;
	register char *cp;

	if(_sp == 1)
		if(optind >= argc ||
		   argv[optind][0] != '-' || argv[optind][1] == '\0')
			return(-1);
		else if(strcmp(argv[optind], "--") == 0) {
			optind++;
			return(-1);
		}
	optopt = c = (unsigned char)argv[optind][_sp];
	if(c == ':' || (cp=strchr(opts, c)) == NULL) {
		/*
		 *  X/Open XPG4 change:
		 *	If 'getopt' detects a missing option-argument, it
		 *	returns the colon character (:) if the first character
		 *	of 'opts' was a colon, or a question-mark (?) character
		 *	otherwise.  If the application has not set the variable
		 *	'opterr' to zero and the first character of 'opts' is
		 *	not a colon, 'getopt' also prints a diagnostic message
		 *	to 'stderr'.
		 */
		if(opterr && *opts != ':')
			pr_error(argv[0], "uxlibc:1:Illegal option -- %c\n", c);
		if(argv[optind][++_sp] == '\0') {
			optind++;
			_sp = 1;
		}
		return('?');
	}
	if(*++cp == ':') {
		if(argv[optind][_sp+1] != '\0')
			optarg = &argv[optind++][_sp+1];
		else if(++optind >= argc) {
		/*
		 *  X/Open XPG4 change:
		 *	If 'getopt' detects a missing option-argument, it
		 *	returns the colon character (:) if the first character
		 *	of 'opts' was a colon, or a question-mark (?) character
		 *	otherwise.  If the application has not set the variable
		 *	'opterr' to zero and the first character of 'opts' is
		 *	not a colon, 'getopt' also prints a diagnostic message
		 *	to 'stderr'.
		 */
		optarg = NULL;
		if (*opts == ':') {
			_sp = 1;
			return(':');
		}
		if (opterr)
			pr_error(argv[0],
				"uxlibc:2:Option requires an argument -- %c\n",c);
			_sp = 1;
			return('?');
		} else
			optarg = argv[optind++];
		_sp = 1;
	} else {
		if(argv[optind][++_sp] == '\0') {
			_sp = 1;
			optind++;
		}
		optarg = NULL;
	}
	return(c);
}

/*
 * reset everyting.
 * We opterr alone since thats already under user control
 */
void
getoptreset(void)
{
	optind = 1;
	optarg = NULL;
	optopt = 0;
	_sp = 1;
}

/* Print the error message.
 * uses argv[0] as label if setlabel() was not called before
 */
static void
pr_error(const char *cmdname, const char *msg, int opt)
{
#ifdef NOTDEF
	int tmp_label;

	if (!__pfmt_label[0]){
		const char *lab = cmdname;
		if (strlen(lab) > (size_t)25) {
			lab = (char *)strrchr(lab, '/');
			if (lab)
				lab++;
			else
				lab = cmdname;
		}
		tmp_label = 1;
		(void)setlabel(lab);
	}
	else
		tmp_label = 0;

	(void)pfmt(stderr, MM_ERROR, msg, opt);

	if (tmp_label)
		__pfmt_label[0] = '\0';
#else
    fprintf(stderr, msg, opt);
#endif
}
