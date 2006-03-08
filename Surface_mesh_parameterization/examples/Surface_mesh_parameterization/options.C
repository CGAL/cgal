// ****************************************************************************
// ^FILE: options.c - implement the functions defined in <options.h>
//
// ^HISTORY:
//    01/16/92	Brad Appleton	<bradapp@enteract.com>	Created
//
//    03/23/93	Brad Appleton	<bradapp@enteract.com>
//    - Added OptIstreamIter class
//
//    10/08/93	Brad Appleton	<bradapp@enteract.com>
//    - Added "hidden" options
//
//    02/08/94	Brad Appleton	<bradapp@enteract.com>
//    - Added "OptionSpec" class
//    - Permitted use of stdio instead of iostreams via #ifdef USE_STDIO
//
//    03/08/94	Brad Appleton	<bradapp@enteract.com>
//    - completed support for USE_STDIO
//    - added #ifdef NO_USAGE for people who always want to print their own
//    - Fixed stupid NULL pointer error in OptionsSpec class
//
//    07/31/97	Brad Appleton	<bradapp@enteract.com>
//    - Added PARSE_POS control flag and POSITIONAL return value.
// ^^**************************************************************************

//#ifdef USE_STDIO
//# include <stdio.h>
//#else
//# include <iostream>
//using namespace std;
//#endif

// #include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "options.h"

// LS 03/2006: include stdlib.h instead of defining directly exit()
#include <stdlib.h>
//extern "C" {
//   void  exit(int);
//}

// LS 08/2005: removed unused variable ident
//static const char ident[] = "@(#)Options  1.05" ;

   // I need a portable version of "tolower" that does NOT modify
   // non-uppercase characters.
   //
#define  TOLOWER(c)  (isupper(c) ? tolower(c) : c)

   // Use this to shut the compiler up about NULL strings
#define  NULLSTR  (char *)NULL

// ******************************************************** insertion operators

  // If you are using <stdio.h> then you need this stuff!
  // If you are using <iostream.h> then #ifdef this stuff out
  //


#ifdef  USE_STDIO

   // Implement just enough of ostream to get this file to compile
   //

static const char endl = '\n' ;

class  ostream {
public:
   ostream(FILE * fileptr) : fp(fileptr) {}

   ostream &
   operator<<(char ch);

   ostream &
   operator<<(const char * str);

   ostream &
   write(const char * buf, unsigned bufsize);

private:
   FILE * fp;
} ;

ostream &
ostream::operator<<(char ch) {
   fputc(ch, fp);
   return *this;
}

ostream &
ostream::operator<<(const char * str) {
   fputs(str, fp);
   return *this;
}

ostream &
ostream::write(const char * buf, unsigned ) {
   fputs(buf, fp);
   return *this;
}

static  ostream  cerr(stderr);
static  ostream  cout(stdout);

#endif  /* USE_STDIO */

// ************************************************************** OptIter

OptIter::~OptIter(void) {}

const char *
OptIter::operator()(void)  {
   const char * elt = curr();
   (void) next();
   return  elt;
}

// ************************************************************** OptIterRwd

OptIterRwd::OptIterRwd(void) {}

OptIterRwd::~OptIterRwd(void) {}

// ************************************************************** OptArgvIter

OptArgvIter::~OptArgvIter(void) {}

const char *
OptArgvIter::curr(void) {
   return ((ndx == ac) || (av[ndx] == NULL)) ? NULLSTR : av[ndx];
}

void
OptArgvIter::next(void) {
   if ((ndx != ac) && av[ndx]) ++ndx;
}

const char *
OptArgvIter::operator()(void) {
   return ((ndx == ac) || (av[ndx] == NULL)) ? NULLSTR : av[ndx++];
}

void
OptArgvIter::rewind(void) { ndx = 0; }

// ************************************************************** OptStrTokIter

static const char WHITESPACE[] = " \t\n\r\v\f" ;
const char * OptStrTokIter::default_delims = WHITESPACE ;

OptStrTokIter::OptStrTokIter(const char * tokens, const char * delimiters)
   : len(unsigned(strlen(tokens))), str(tokens), seps(delimiters),
     cur(NULLSTR), tokstr(NULLSTR)
{
   if (seps == NULL)  seps = default_delims;
   tokstr = new char[len + 1];
   (void) ::strcpy(tokstr, str);
   cur = ::strtok(tokstr, seps);
}


OptStrTokIter::~OptStrTokIter(void) { delete [] tokstr; }

const char *
OptStrTokIter::curr(void) { return cur; }

void
OptStrTokIter::next(void) { if (cur) cur = ::strtok(NULL, seps); }

const char *
OptStrTokIter::operator()(void) {
   const char * elt = cur;
   if (cur) cur = ::strtok(NULL, seps);
   return  elt;
}

void
OptStrTokIter::rewind(void) {
   (void) ::strcpy(tokstr, str);
   cur = ::strtok(tokstr, seps);
}

// ************************************************************* OptIstreamIter

#ifdef vms
   enum { c_COMMENT = '!' } ;
#else
   enum { c_COMMENT = '#' } ;
#endif

const unsigned  OptIstreamIter::MAX_LINE_LEN = 1024 ;

   // Constructor
OptIstreamIter::OptIstreamIter(istream & input) : is(input), tok_iter(NULL)
{
#ifdef  USE_STDIO
   fprintf(stderr, "%s: Can't use OptIstreamIter class:\n",
                   "OptIstreamIter::OptIstreamIter");
   fprintf(stderr, "\tOptions(3C++) was compiled with USE_STDIO #defined.\n");
   exit(-1);
#endif  /* USE_STDIO */
}

   // Destructor
OptIstreamIter::~OptIstreamIter(void) {
   delete  tok_iter;
}

const char *
OptIstreamIter::curr(void) {
#ifdef  USE_STDIO
   return  NULLSTR;
#else
   const char * result = NULLSTR;
   if (tok_iter)  result = tok_iter->curr();
   if (result)  return  result;
   fill();
   return (! is) ? NULLSTR : tok_iter->curr();
#endif  /* USE_STDIO */
}

void
OptIstreamIter::next(void) {
#ifdef  USE_STDIO
   return;
#else
   const char * result = NULLSTR;
   if (tok_iter)  result = tok_iter->operator()();
   if (result)  return;
   fill();
   if (! is) tok_iter->next();
#endif  /* USE_STDIO */
}

const char *
OptIstreamIter::operator()(void) {
#ifdef  USE_STDIO
   return  NULLSTR;
#else
   const char * result = NULLSTR;
   if (tok_iter)  result = tok_iter->operator()();
   if (result)  return  result;
   fill();
   return (! is) ? NULLSTR : tok_iter->operator()();
#endif  /* USE_STDIO */
}

   // What we do is this: for each line of text in the istream, we use
   // a OptStrTokIter to iterate over each token on the line.
   //
   // If the first non-white character on a line is c_COMMENT, then we
   // consider the line to be a comment and we ignore it.
   //
void
OptIstreamIter::fill(void) {
#ifdef USE_STDIO
   return;
#else
   char buf[OptIstreamIter::MAX_LINE_LEN];
   do {
      *buf = '\0';
      is.getline(buf, sizeof(buf));
      char * ptr = buf;
      while (isspace(*ptr)) ++ptr;
      if (*ptr && (*ptr != c_COMMENT)) {
         delete  tok_iter;
         tok_iter = new OptStrTokIter(ptr);
         return;
      }
   } while (is);
#endif  /* USE_STDIO */
}

// **************************************************** Options class utilities

   // Is this option-char null?
inline static int
isNullOpt(char optchar) {
   return  ((! optchar) || isspace(optchar) || (! isprint(optchar)));
}

   // Check for explicit "end-of-options"
inline static int
isEndOpts(const char * token) {
   return ((token == NULL) || (! ::strcmp(token, "--"))) ;
}

   // See if an argument is an option
inline static int
isOption(unsigned  flags, const char * arg) {
   return  (((*arg != '\0') || (arg[1] != '\0')) &&
            ((*arg == '-')  || ((flags & Options::PLUS) && (*arg == '+')))) ;
}

   // See if we should be parsing only options or if we also need to
   // parse positional arguments
inline static int
isOptsOnly(unsigned  flags) {
   return  (flags & Options::PARSE_POS) ? 0 : 1;
}

   // return values for a keyword matching function
enum kwdmatch_t { NO_MATCH, PARTIAL_MATCH, EXACT_MATCH } ;

// ---------------------------------------------------------------------------
// ^FUNCTION: kwdmatch - match a keyword
//
// ^SYNOPSIS:
//    static kwdmatch_t kwdmatch(src, attempt, len)
//
// ^PARAMETERS:
//    char * src -- the actual keyword to match
//    char * attempt -- the possible keyword to compare against "src"
//    int len -- number of character of "attempt" to consider
//               (if 0 then we should use all of "attempt")
//
// ^DESCRIPTION:
//    See if "attempt" matches some prefix of "src" (case insensitive).
//
// ^REQUIREMENTS:
//    - attempt should be non-NULL and non-empty
//
// ^SIDE-EFFECTS:
//    None.
//
// ^RETURN-VALUE:
//    An enumeration value of type kwdmatch_t corresponding to whether
//    We had an exact match, a partial match, or no match.
//
// ^ALGORITHM:
//    Trivial
// ^^-------------------------------------------------------------------------
static kwdmatch_t
kwdmatch(const char * src, const char * attempt, int len =0) {
   int  i;

   if (src == attempt)  return  EXACT_MATCH ;
   if ((src == NULL) || (attempt == NULL))  return  NO_MATCH ;
   if ((! *src) && (! *attempt))  return  EXACT_MATCH ;
   if ((! *src) || (! *attempt))  return  NO_MATCH ;

   for (i = 0 ; ((i < len) || (len == 0)) &&
                (attempt[i]) && (attempt[i] != ' ') ; i++) {
      if (TOLOWER(src[i]) != TOLOWER(attempt[i]))  return  NO_MATCH ;
   }

   return  (src[i]) ? PARTIAL_MATCH : EXACT_MATCH ;
}

// **************************************************************** OptionSpec

   // Class that represents an option-specification
   //    *NOTE*:: Assumes that the char-ptr given to the constructor points
   //             to storage that will NOT be modified and whose lifetime will
   //             be as least as long as the OptionSpec object we construct.
   //
class OptionSpec {
public:
   OptionSpec(const char * decl =NULLSTR)
      : hidden(0), spec(decl)
   {
      if (spec == NULL)  spec = NULL_spec;
      CheckHidden();
   }

   OptionSpec(const OptionSpec & cp) : hidden(cp.hidden), spec(cp.spec) {}

   // NOTE: use default destructor!

      // Assign to another OptionSpec
   OptionSpec &
   operator=(const OptionSpec & cp) {
      if (this != &cp) {
         spec = cp.spec;
         hidden = cp.hidden;
      }
      return *this;
   }

      // Assign to a string
   OptionSpec &
   operator=(const char * decl) {
      if (spec != decl) {
         spec = decl;
         hidden = 0;
         CheckHidden();
      }
      return *this;
   }

      // Convert to char-ptr by returning the original declaration-string
   operator const char*() { return  isHiddenOpt() ? (spec - 1) : spec; }

      // Is this option NULL?
   int
   isNULL(void) const { return ((spec == NULL) || (spec == NULL_spec)); }

      // Is this options incorrectly specified?
   int
   isSyntaxError(const char * name) const;

      // See if this is a Hidden option
   int
   isHiddenOpt(void) const { return  hidden; }

      // Get the corresponding option-character
   char
   OptChar(void) const { return  *spec; }

      // Get the corresponding long-option string
   const char *
   LongOpt(void) const {
       return  (spec[1] && spec[2] && (! isspace(spec[2]))) ? (spec + 2) : NULLSTR;
   }

      // Does this option require an argument?
   int
   isValRequired(void) const {
      return  ((spec[1] == ':') || (spec[1] == '+'));
   }

      // Does this option take an optional argument?
   int
   isValOptional(void) const {
      return  ((spec[1] == '?') || (spec[1] == '*'));
   }

      // Does this option take no arguments?
   int
   isNoArg(void) const {
      return  ((spec[1] == '|') || (! spec[1]));
   }

      // Can this option take more than one argument?
   int
   isList(void) const {
      return  ((spec[1] == '+') || (spec[1] == '*'));
   }

      // Does this option take any arguments?
   int
   isValTaken(void) const {
      return  (isValRequired() || isValOptional()) ;
   }

      // Format this option in the given buffer
   unsigned
   Format(char * buf, unsigned optctrls) const;

private:
   void
   CheckHidden(void) {
      if ((! hidden) && (*spec == '-')) {
         ++hidden;
         ++spec;
      }
   }

   unsigned     hidden : 1;  // hidden-flag
   const char * spec;        // string specification

   static const char NULL_spec[];
} ;

const char OptionSpec::NULL_spec[] = "\0\0\0" ;

int
OptionSpec::isSyntaxError(const char * name) const {
   int  error = 0;
   if ((! spec) || (! *spec)) {
      cerr << name << ": empty option specifier." << endl;
      cerr << "\tmust be at least 1 character long." << endl;
      ++error;
   } else if (spec[1] && (strchr("|?:*+", spec[1]) == NULL)) {
      cerr << name << ": bad option specifier \"" << spec << "\"." << endl;
      cerr << "\t2nd character must be in the set \"|?:*+\"." << endl;
      ++error;
   }
   return  error;
}

// ---------------------------------------------------------------------------
// ^FUNCTION: OptionSpec::Format - format an option-spec for a usage message
//
// ^SYNOPSIS:
//    unsigned OptionSpec::Format(buf, optctrls) const
//
// ^PARAMETERS:
//    char * buf -- where to print the formatted option
//    unsigned  optctrls -- option-parsing configuration flags
//
// ^DESCRIPTION:
//    Self-explanatory.
//
// ^REQUIREMENTS:
//    - buf must be large enough to hold the result
//
// ^SIDE-EFFECTS:
//    - writes to buf.
//
// ^RETURN-VALUE:
//    Number of characters written to buf.
//
// ^ALGORITHM:
//    Follow along in the source - it's not hard but it is tedious!
// ^^-------------------------------------------------------------------------
unsigned
OptionSpec::Format(char * buf, unsigned optctrls) const {
#ifdef NO_USAGE
   return  (*buf = '\0');
#else
   static  char default_value[] = "<value>";
   if (isHiddenOpt())  return (unsigned)(*buf = '\0');
   char optchar = OptChar();
   const char * longopt = LongOpt();
   char * p = buf ;

   const char * value = NULLSTR;
   int    longopt_len = 0;
   int    value_len = 0;

   if (longopt) {
      value = ::strchr(longopt, ' ');
      longopt_len = (value) ? (value - longopt) : ::strlen(longopt);
   } else {
      value = ::strchr(spec + 1, ' ');
   }
   while (value && (*value == ' '))  ++value;
   if (value && *value) {
      value_len = ::strlen(value);
   } else {
      value = default_value;
      value_len = sizeof(default_value) - 1;
   }

   if ((optctrls & Options::SHORT_ONLY) &&
       ((! isNullOpt(optchar)) || (optctrls & Options::NOGUESSING))) {
      longopt = NULLSTR;
   }
   if ((optctrls & Options::LONG_ONLY) &&
       (longopt || (optctrls & Options::NOGUESSING))) {
      optchar = '\0';
   }
   if (isNullOpt(optchar) && (longopt == NULL)) {
      *buf = '\0';
      return  0;
   }

   *(p++) = '[';

   // print the single character option
   if (! isNullOpt(optchar)) {
      *(p++) = '-';
      *(p++) = optchar;
   }

   if ((! isNullOpt(optchar)) && (longopt))  *(p++) = '|';

   // print the long option
   if (longopt) {
      *(p++) = '-';
      if (! (optctrls & (Options::LONG_ONLY | Options::SHORT_ONLY))) {
         *(p++) = '-';
      }
      strncpy(p, longopt, longopt_len);
      p += longopt_len;
   }

   // print any argument the option takes
   if (isValTaken()) {
      *(p++) = ' ' ;
      if (isValOptional())  *(p++) = '[' ;
      strcpy(p, value);
      p += value_len;
      if (isList()) {
         strcpy(p, " ...");
         p += 4;
      }
      if (isValOptional())  *(p++) = ']' ;
   }

   *(p++) = ']';
   *p = '\0';

   return  (unsigned) strlen(buf);
#endif  /* USE_STDIO */
}

// ******************************************************************* Options

#if (defined(MSWIN) || defined(OS2) || defined(MSDOS))
# define DIR_SEP_CHAR '\\'
#else
# define DIR_SEP_CHAR '/'
#endif

Options::Options(const char * name, const char * const optv[])
   : explicit_end(0), optctrls(DEFAULT), optvec(optv),
     nextchar(NULLSTR), listopt(NULLSTR), cmdname(name)
{
   const char * basename = ::strrchr(cmdname, DIR_SEP_CHAR);
   if (basename)  cmdname = basename + 1;
   check_syntax();
}

Options::~Options(void) {}

   // Make sure each option-specifier has correct syntax.
   //
   // If there is even one invalid specifier, then exit ungracefully!
   //
void
Options::check_syntax(void) const {
   int  errors = 0;
   if ((optvec == NULL) || (! *optvec))  return;

   for (const char * const * optv = optvec ; *optv ; optv++) {
      OptionSpec  optspec = *optv;
      errors += optspec.isSyntaxError(cmdname);
   }
   if (errors)  exit(127);
}

// ---------------------------------------------------------------------------
// ^FUNCTION: Options::match_opt - match an option
//
// ^SYNOPSIS:
//    const char * match_opt(opt, int  ignore_case) const
//
// ^PARAMETERS:
//    char opt -- the option-character to match
//    int  ignore_case -- should we ignore character-case?
//
// ^DESCRIPTION:
//    See if "opt" is found in "optvec"
//
// ^REQUIREMENTS:
//    - optvec should be non-NULL and terminated by a NULL pointer.
//
// ^SIDE-EFFECTS:
//    None.
//
// ^RETURN-VALUE:
//    NULL if no match is found,
//    otherwise a pointer to the matching option-spec.
//
// ^ALGORITHM:
//    foreach option-spec
//       - see if "opt" is a match, if so return option-spec
//    end-for
// ^^-------------------------------------------------------------------------
const char *
Options::match_opt(char opt, int ignore_case) const {
   if ((optvec == NULL) || (! *optvec))  return  NULLSTR;

   for (const char * const * optv = optvec ; *optv ; optv++) {
      OptionSpec  optspec = *optv;
      char optchar = optspec.OptChar();
      if (isNullOpt(optchar))  continue;
      if (opt == optchar) {
         return  optspec;
      } else if (ignore_case && (TOLOWER(opt) == TOLOWER(optchar))) {
         return  optspec;
      }
   }

   return  NULLSTR;  // not found
}

// ---------------------------------------------------------------------------
// ^FUNCTION: Options::match_longopt - match a long-option
//
// ^SYNOPSIS:
//   const char * Options::match_longopt(opt, len, ambiguous)
//
// ^PARAMETERS:
//    char * opt -- the long-option to match
//    int len -- the number of character of "opt" to match
//    int & ambiguous -- set by this routine before returning.
//
// ^DESCRIPTION:
//    Try to match "opt" against some unique prefix of a long-option
//    (case insensitive).
//
// ^REQUIREMENTS:
//    - optvec should be non-NULL and terminated by a NULL pointer.
//
// ^SIDE-EFFECTS:
//    - *ambiguous is set to '1' if "opt" matches >1 long-option
//      (otherwise it is set to 0).
//
// ^RETURN-VALUE:
//    NULL if no match is found,
//    otherwise a pointer to the matching option-spec.
//
// ^ALGORITHM:
//    ambiguous is FALSE
//    foreach option-spec
//       if we have an EXACT-MATCH, return the option-spec
//       if we have a partial-match then
//          if we already had a previous partial match then
//             set ambiguous = TRUE and return NULL
//          else
//             remember this options spec and continue matching
//          end-if
//       end-if
//    end-for
//    if we had exactly 1 partial match return it, else return NULL
// ^^-------------------------------------------------------------------------
const char *
Options::match_longopt(const char * opt, int  len, int & ambiguous) const {
   kwdmatch_t  result;
   const char * matched = NULLSTR ;

   ambiguous = 0;
   if ((optvec == NULL) || (! *optvec))  return  NULLSTR;

   for (const char * const * optv = optvec ; *optv ; optv++) {
      OptionSpec  optspec = *optv;
      const char * longopt = optspec.LongOpt();
      if (longopt == NULL)  continue;
      result = kwdmatch(longopt, opt, len);
      if (result == EXACT_MATCH) {
         return  optspec;
      } else if (result == PARTIAL_MATCH) {
         if (matched) {
            ++ambiguous;
            return  NULLSTR;
         } else {
            matched = optspec;
         }
      }
   }//for

   return  matched;
}

// ---------------------------------------------------------------------------
// ^FUNCTION: Options::parse_opt - parse an option
//
// ^SYNOPSIS:
//    int Options::parse_opt(iter, optarg)
//
// ^PARAMETERS:
//    OptIter & iter -- option iterator
//    const char * & optarg -- where to store any option-argument
//
// ^DESCRIPTION:
//    Parse the next option in iter (advancing as necessary).
//    Make sure we update the nextchar pointer along the way. Any option
//    we find should be returned and optarg should point to its argument.
//
// ^REQUIREMENTS:
//    - nextchar must point to the prospective option character
//
// ^SIDE-EFFECTS:
//    - iter is advanced when an argument completely parsed
//    - optarg is modified to point to any option argument
//    - if Options::QUIET is not set, error messages are printed on cerr
//
// ^RETURN-VALUE:
//    'c' if the -c option was matched (optarg points to its argument)
//    BADCHAR if the option is invalid (optarg points to the bad
//                                                   option-character).
//
// ^ALGORITHM:
//    It gets complicated -- follow the comments in the source.
// ^^-------------------------------------------------------------------------
int
Options::parse_opt(OptIter & iter, const char * & optarg) {
   listopt = NULLSTR;  // reset the list pointer

   if ((optvec == NULL) || (! *optvec))  return  Options::ENDOPTS;

      // Try to match a known option
   OptionSpec optspec = match_opt(*(nextchar++), (optctrls & Options::ANYCASE));

      // Check for an unknown option
   if (optspec.isNULL()) {
      // See if this was a long-option in disguise
      if (! (optctrls & Options::NOGUESSING)) {
         unsigned  save_ctrls = optctrls;
         const char * save_nextchar = nextchar;
         nextchar -= 1;
         optctrls |= (Options::QUIET | Options::NOGUESSING);
         int  optchar = parse_longopt(iter, optarg);
         optctrls = save_ctrls;
         if (optchar > 0) {
            return  optchar;
         } else {
            nextchar = save_nextchar;
         }
      }
      if (! (optctrls & Options::QUIET)) {
         cerr << cmdname << ": unknown option -"
              << *(nextchar - 1) << "." << endl ;
      }
      optarg = (nextchar - 1);  // record the bad option in optarg
      return  Options::BADCHAR;
   }

      // If no argument is taken, then leave now
   if (optspec.isNoArg()) {
      optarg = NULLSTR;
      return  optspec.OptChar();
   }

      // Check for argument in this arg
   if (*nextchar) {
      optarg = nextchar; // the argument is right here
      nextchar = NULLSTR;   // we've exhausted this arg
      if (optspec.isList())  listopt = optspec ;  // save the list-spec
      return  optspec.OptChar();
   }

      // Check for argument in next arg
   const char * nextarg = iter.curr();
   if ((nextarg != NULL)  &&
       (optspec.isValRequired() || (! isOption(optctrls, nextarg)))) {
      optarg = nextarg;    // the argument is here
      iter.next();         // end of arg - advance
      if (optspec.isList())  listopt = optspec ;  // save the list-spec
      return  optspec.OptChar();
   }

     // No argument given - if its required, thats an error
   optarg = NULLSTR;
   if (optspec.isValRequired() &&  !(optctrls & Options::QUIET)) {
      cerr << cmdname << ": argument required for -" << optspec.OptChar()
           << " option." << endl ;
   }
   return  optspec.OptChar();
}

// ---------------------------------------------------------------------------
// ^FUNCTION: Options::parse_longopt - parse a long-option
//
// ^SYNOPSIS:
//    int Options::parse_longopt(iter, optarg)
//
// ^PARAMETERS:
//    OptIter & iter -- option iterator
//    const char * & optarg -- where to store any option-argument
//
// ^DESCRIPTION:
//    Parse the next long-option in iter (advancing as necessary).
//    Make sure we update the nextchar pointer along the way. Any option
//    we find should be returned and optarg should point to its argument.
//
// ^REQUIREMENTS:
//    - nextchar must point to the prospective option character
//
// ^SIDE-EFFECTS:
//    - iter is advanced when an argument completely parsed
//    - optarg is modified to point to any option argument
//    - if Options::QUIET is not set, error messages are printed on cerr
//
// ^RETURN-VALUE:
//    'c' if the the long-option corresponding to the -c option was matched
//         (optarg points to its argument)
//    BADKWD if the option is invalid (optarg points to the bad long-option
//                                                                     name).
//
// ^ALGORITHM:
//    It gets complicated -- follow the comments in the source.
// ^^-------------------------------------------------------------------------
int
Options::parse_longopt(OptIter & iter, const char * & optarg) {
   int  len = 0, ambiguous = 0;

   listopt = NULLSTR ;  // reset the list-spec

   if ((optvec == NULL) || (! *optvec))  return  Options::ENDOPTS;

      // if a value is supplied in this argv element, get it now
   const char * val = strpbrk(nextchar, ":=") ;
   if (val) {
      len = val - nextchar ;
      ++val ;
   }

      // Try to match a known long-option
   OptionSpec  optspec = match_longopt(nextchar, len, ambiguous);

      // Check for an unknown long-option
   if (optspec.isNULL()) {
      // See if this was a short-option in disguise
      if ((! ambiguous) && (! (optctrls & Options::NOGUESSING))) {
         unsigned  save_ctrls = optctrls;
         const char * save_nextchar = nextchar;
         optctrls |= (Options::QUIET | Options::NOGUESSING);
         int  optchar = parse_opt(iter, optarg);
         optctrls = save_ctrls;
         if (optchar > 0) {
            return  optchar;
         } else {
            nextchar = save_nextchar;
         }
      }
      if (! (optctrls & Options::QUIET)) {
         cerr << cmdname << ": " << ((ambiguous) ? "ambiguous" : "unknown")
              << " option "
              << ((optctrls & Options::LONG_ONLY) ? "-" : "--")
              << nextchar << "." << endl ;
      }
      optarg = nextchar;  // record the bad option in optarg
      nextchar = NULLSTR;    // we've exhausted this argument
      return  (ambiguous) ? Options::AMBIGUOUS : Options::BADKWD;
   }

      // If no argument is taken, then leave now
   if (optspec.isNoArg()) {
      if ((val) && ! (optctrls & Options::QUIET)) {
         cerr << cmdname << ": option "
              << ((optctrls & Options::LONG_ONLY) ? "-" : "--")
              << optspec.LongOpt() << " does NOT take an argument." << endl ;
      }
      optarg = val;     // record the unexpected argument
      nextchar = NULLSTR;  // we've exhausted this argument
      return  optspec.OptChar();
   }

      // Check for argument in this arg
   if (val) {
      optarg = val;      // the argument is right here
      nextchar = NULLSTR;   // we exhausted the rest of this arg
      if (optspec.isList())  listopt = optspec ;  // save the list-spec
      return  optspec.OptChar();
   }

      // Check for argument in next arg
   const char * nextarg = iter.curr();  // find the next argument to parse
   if ((nextarg != NULL)  &&
       (optspec.isValRequired() || (! isOption(optctrls, nextarg)))) {
      optarg = nextarg;        // the argument is right here
      iter.next();             // end of arg - advance
      nextchar = NULLSTR;         // we exhausted the rest of this arg
      if (optspec.isList())  listopt = optspec ;  // save the list-spec
      return  optspec.OptChar();
   }

     // No argument given - if its required, thats an error
   optarg = NULLSTR;
   if (optspec.isValRequired() &&  !(optctrls & Options::QUIET)) {
      const char * longopt = optspec.LongOpt();
      const char * spc = ::strchr(longopt, ' ');
      int  longopt_len;
      if (spc) {
         longopt_len = spc - longopt;
      } else {
         longopt_len = ::strlen(longopt);
      }
      cerr << cmdname << ": argument required for "
           << ((optctrls & Options::LONG_ONLY) ? "-" : "--");
      cerr.write(longopt, longopt_len) << " option." << endl ;
   }
   nextchar = NULLSTR;           // we exhausted the rest of this arg
   return  optspec.OptChar();
}

// ---------------------------------------------------------------------------
// ^FUNCTION: Options::usage - print usage
//
// ^SYNOPSIS:
//    void Options::usage(os, positionals)
//
// ^PARAMETERS:
//    ostream & os -- where to print the usage
//    char * positionals -- command-line syntax for any positional args
//
// ^DESCRIPTION:
//    Print command-usage (using either option or long-option syntax) on os.
//
// ^REQUIREMENTS:
//    os should correspond to an open output file.
//
// ^SIDE-EFFECTS:
//    Prints on os
//
// ^RETURN-VALUE:
//    None.
//
// ^ALGORITHM:
//    Print usage on os, wrapping long lines where necessary.
// ^^-------------------------------------------------------------------------
void
Options::usage(ostream & os, const char * positionals) const {
#ifdef NO_USAGE
   return;
#else
   const char * const * optv = optvec;
   unsigned  cols = 79;
   int  first, nloop;
   char  buf[256] ;

   if ((optv == NULL) || (! *optv))  return;

      // print first portion "usage: progname"
   os << "usage: " << cmdname ;
   unsigned  ll = strlen(cmdname) + 7;

      // save the current length so we know how much space to skip for
      // subsequent lines.
      //
   unsigned  margin = ll + 1;

      // print the options and the positional arguments
   for (nloop = 0, first = 1 ; !nloop ; optv++, first = 0) {
      unsigned  len;
      OptionSpec   optspec = *optv;

         // figure out how wide this parameter is (for printing)
      if (! *optv) {
         len = strlen(positionals);
         ++nloop;  // terminate this loop
      } else {
         if (optspec.isHiddenOpt())  continue;
         len = optspec.Format(buf, optctrls);
      }

      //  Will this fit?
      if ((ll + len + 1) > (cols - first)) {
         os << '\n' ;  // No - start a new line;
#ifdef USE_STDIO
         for (int _i_ = 0; _i_ < margin; ++_i_)  os << " ";
#else
         os.width(margin); os << "" ;
#endif
         ll = margin;
      } else {
         os << ' ' ;  // Yes - just throw in a space
         ++ll;
      }
      ll += len;
      os << ((nloop) ? positionals : buf) ;
   }// for each ad

   os << endl ;
#endif  /* NO_USAGE */
}


// ---------------------------------------------------------------------------
// ^FUNCTION: Options::operator() - get options from the command-line
//
// ^SYNOPSIS:
//   int Options::operator()(iter, optarg)
//
// ^PARAMETERS:
//    OptIter & iter -- option iterator
//    const char * & optarg -- where to store any option-argument
//
// ^DESCRIPTION:
//    Parse the next option in iter (advancing as necessary).
//    Make sure we update the nextchar pointer along the way. Any option
//    we find should be returned and optarg should point to its argument.
//
// ^REQUIREMENTS:
//    None.
//
// ^SIDE-EFFECTS:
//    - iter is advanced when an argument is completely parsed
//    - optarg is modified to point to any option argument
//    - if Options::QUIET is not set, error messages are printed on cerr
//
// ^RETURN-VALUE:
//     0 if all options have been parsed.
//    'c' if the the option or long-option corresponding to the -c option was
//         matched (optarg points to its argument).
//    BADCHAR if the option is invalid (optarg points to the bad option char).
//    BADKWD if the option is invalid (optarg points to the bad long-opt name).
//    AMBIGUOUS if an ambiguous keyword name was given (optarg points to the
//         ambiguous keyword name).
//    POSITIONAL if PARSE_POS was set and the current argument is a positional
//         parameter (in which case optarg points to the positional argument).
//
// ^ALGORITHM:
//    It gets complicated -- follow the comments in the source.
// ^^-------------------------------------------------------------------------
int
Options::operator()(OptIter & iter, const char * & optarg) {
   int parse_opts_only = isOptsOnly(optctrls);
   if (parse_opts_only)  explicit_end = 0;

      // See if we have an option left over from before ...
   if ((nextchar) && *nextchar) {
      return  parse_opt(iter, optarg);
   }

      // Check for end-of-options ...
   const char * arg = NULLSTR;
   int get_next_arg = 0;
   do {
      arg = iter.curr();
      get_next_arg = 0;
      if (arg == NULL) {
         listopt = NULLSTR;
         return  Options::ENDOPTS;
      } else if ((! explicit_end) && isEndOpts(arg)) {
         iter.next();   // advance past end-of-options arg
         listopt = NULLSTR;
         explicit_end = 1;
         if (parse_opts_only)  return  Options::ENDOPTS;
         get_next_arg = 1;  // make sure we look at the next argument.
      }
   } while (get_next_arg);

      // Do we have a positional arg?
   if ( explicit_end || (! isOption(optctrls, arg)) ) {
      if (parse_opts_only) {
         return  Options::ENDOPTS;
      } else {
         optarg = arg;  // set optarg to the positional argument
         iter.next();   // advance iterator to the next argument
         return  Options::POSITIONAL;
      }
   }

   iter.next();  // pass the argument that arg already points to

      // See if we have a long option ...
   if (! (optctrls & Options::SHORT_ONLY)) {
      if ((*arg == '-') && (arg[1] == '-')) {
         nextchar = arg + 2;
         return  parse_longopt(iter, optarg);
      } else if ((optctrls & Options::PLUS) && (*arg == '+')) {
         nextchar = arg + 1;
         return  parse_longopt(iter, optarg);
      }
   }
   if (*arg == '-') {
      nextchar = arg + 1;
      if (optctrls & Options::LONG_ONLY) {
         return  parse_longopt(iter, optarg);
      } else {
         return  parse_opt(iter, optarg);
      }
   }

      // If we get here - it is because we have a list value
   OptionSpec  optspec = listopt;
   optarg = arg ;        // record the list value
   return  optspec.OptChar() ;
}

