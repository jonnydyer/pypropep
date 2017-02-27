/* file getopt.c                                                   */
/*                Author: Peter Wilson                             */
/*                        Catholic University and NIST             */
/*                        pwilson@cme.nist.gov                     */
/*                                                                 */
/* getopt()  from Don Libes "Obfuscated C" */

#include <stdio.h>
#include <string.h>

/* getopt()  -- parse command line arguments */
/* Original Author: AT&T */
/* This version from Don Libes "Obfuscated C and Other Mysteries" */
/* John Wiley & Sons, 1993. Chapter 6 */

#define ERR(s, c)    if (opterr) {\
     char errbuf[3];\
     errbuf[0] = c; errbuf[1] = '\n'; errbuf[2] = '\0'; \
     fprintf(stderr, "%s", argv[0]);\
     fprintf(stderr, "%s", s);\
     fprintf(stderr, "%s", errbuf); }


int   opterr = 1;  /* getopt prints errors if this is one */
int   optind = 1;  /* token pointer */
int   optopt;      /* option character passed back to user */
char *optarg;      /* flag argument (or value) */

/* return option option character, EOF if no more or ? if problem */
int getopt(int argc, char **argv, char *opts) /* opts: option string */
{
  static   int   sp = 1; /* character index in current token */
  register char *cp;     /* pointer into current token */

  if (sp == 1) {
    /* check for more flag-like tokens */
    if(optind >= argc ||
       argv[optind][0] != '-' || argv[optind][1] == '\0')
      return(EOF);
    else if(strcmp(argv[optind], "--") == 0) {
      optind++;
      return(EOF);
    }
  }
  optopt = argv[optind][sp];
  if(optopt == ':' || ( cp = strchr(opts, optopt)) == 0) {
    ERR(": illegal option -- ", optopt);
    /* if no chars left in this token, move to next token */
    if(argv[optind][++sp] == '\0') {
      optind++;
      sp = 1;
    }
    return('?');
  }
  
  if(*++cp == ':') {/* if a value is expected, get it */
    if(argv[optind][sp+1] != '\0')
      /* flag value is rest of current token */
      optarg = &argv[optind++][sp+1];
    else if (++optind >= argc) {
      ERR(": option requires an argument -- ", optopt);
      sp = 1;
      return('?');
    } else
      /* flag value is next token */
      optarg = argv[optind++];
    sp = 1;
  } else {
    /* set up to look at next char in token, next time */
    if(argv[optind][++sp] == '\0') {
      /* no m ore in current token, set up next token */
      sp = 1;
      optind++;
    }
    optarg = 0;
  }
  return(optopt); /* return current flag character found */
}
   




