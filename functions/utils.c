// This file contains some general code for use by the rest of carat
#include "config.h"
#include "utils.h"

#include <stdio.h>
#include <string.h>

void *xmalloc_(size_t size, const char *funcname)
{
    void *ptr = malloc(size);
    if (ptr == NULL)
    {
        perror(funcname);
        exit(2);
    }
    return ptr;
}

void *xrealloc_(void *ptr, size_t size, const char *funcname)
{
    ptr = realloc(ptr, size);
    if (ptr == NULL)
    {
        perror(funcname);
        exit(2);
    }
    return ptr;
}

#ifndef HAVE_STRLCPY

size_t strlcpy(char * dst, const char * src, size_t len)
{
    /* Keep a copy of the original src. */
    const char * const orig_src = src;

    /* If a non-empty len was specified, we can actually copy some data. */
    if (len > 0) {
        /* Copy up to len-1 bytes (reserve one for the terminating zero). */
        while (--len > 0) {
            /* Copy from src to dst; if we reach the string end, we are
               done and can simply return the total source string length */
            if ((*dst++ = *src++) == 0) {
                /* return length of source string without the zero byte */
                return src - orig_src - 1;
            }
        }

        /* If we got here, then we used up the whole buffer and len is zero.
           We must make sure to terminate the destination string. */
        *dst = 0;
    }

    /* in the end, we must return the length of the source string, no
       matter whether we completely copied or not; so advance src
       till its terminator is reached */
    while (*src++)
        ;

    /* return length of source string without the zero byte */
    return src - orig_src - 1;
}

#endif /* !HAVE_STRLCPY */


#ifndef HAVE_STRLCAT

size_t strlcat(char * dst, const char * src, size_t len)
{
    /* Keep a copy of the original dst. */
    const char * const orig_dst = dst;

    /* Find the end of the dst string, so that we can append after it. */
    while (*dst != 0 && len > 0) {
        dst++;
        len--;
    }

    /* We can only append anything if there is free space left in the
       destination buffer. */
    if (len > 0) {
        /* One byte goes away for the terminating zero. */
        len--;

        /* Do the actual work and append from src to dst, until we either
           appended everything, or reached the dst buffer's end. */
        while (*src != 0 && len > 0) {
            *dst++ = *src++;
            len--;
        }

        /* Terminate, terminate, terminate! */
        *dst = 0;
    }

    /* Compute the final result. */
    return (dst - orig_dst) + strlen(src);
}

#endif /* !HAVE_STRLCAT */
