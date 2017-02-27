
#ifdef BORLAND

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* Compare two strings, ignoring case */
static int strcasecmp(const char *s1, const char *s2)
{
  unsigned int i, len = strlen(s1)+1;
  char c1, c2;

  for (i = 0; i < len; ++i)
  {
    c1 = tolower(s1[i]);
    c2 = tolower(s2[i]);
    if (c1 != c2)
    {
      return ((c1 > c2) ? (1) : (-1));
    }
  }
  return 0;
}

/* Compare two strings up to n characters, ignoring case */
static int strncasecmp(const char *s1, const char *s2, size_t sz)
{
  unsigned int i, len = strlen(s1)+1;
  char c1, c2;

  if (sz < len) len = sz;
  
  for (i = 0; i < len; ++i)
  {
    c1 = tolower(s1[i]);
    c2 = tolower(s2[i]);
    if (c1 != c2)
    {
      return ((c1 > c2) ? (1) : (-1));
    }
  }
  return 0;
}

#endif
