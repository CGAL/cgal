//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_NEW_H
#define BOOST_UBLAS_NEW_H

void *operator new (std::size_t size) {
    size += 16;
    void *p = malloc (size);
    char offset = (unsigned long) p % 16;
    if (offset == 0)
        offset = 16;
    (char *) p += offset;
    *((char *) p - 1) = offset;
    return p;
}
void operator delete (void *p) {
    char offset = *((char *) p - 1);
    (char *) p -= offset;
    free (p);
}
void *operator new [] (std::size_t size) {
    size += 16;
    void *p = malloc (size);
    char offset = (unsigned long) p % 16;
    if (offset == 0)
        offset = 16;
    (char *) p += offset;
    *((char *) p - 1) = offset;
    return p;
}
void operator delete [] (void *p) {
    char offset = *((char *) p - 1);
    (char *) p -= offset;
    free (p);
}

#endif