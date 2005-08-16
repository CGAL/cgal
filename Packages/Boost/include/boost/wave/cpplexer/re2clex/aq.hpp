/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/
    
    Copyright (c) 2001 Daniel C. Nuffer.
    Copyright (c) 2001-2005 Hartmut Kaiser. 
    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(AQ_HPP_A21D9145_B643_44C0_81E7_DB346DD67EE1_INCLUDED)
#define AQ_HPP_A21D9145_B643_44C0_81E7_DB346DD67EE1_INCLUDED

#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace cpplexer {
namespace re2clex {

typedef std::size_t aq_stdelement;

typedef struct tag_aq_queuetype
{
    std::size_t head;
    std::size_t tail;
    std::size_t size;
    std::size_t max_size;
    aq_stdelement* queue;
} aq_queuetype;

typedef aq_queuetype* aq_queue;

int aq_enqueue(aq_queue q, aq_stdelement e);
int aq_enqueue_front(aq_queue q, aq_stdelement e);
int aq_serve(aq_queue q, aq_stdelement *e);
int aq_pop(aq_queue q);
#define AQ_EMPTY(q) (q->size == 0)
#define AQ_FULL(q) (q->size == q->max_size)
aq_queue aq_create(void);
void aq_terminate(aq_queue q);
int aq_grow(aq_queue q);

///////////////////////////////////////////////////////////////////////////////
}   // namespace re2clex
}   // namespace cpplexer
}   // namespace wave
}   // namespace boost 

#endif // !defined(AQ_HPP_A21D9145_B643_44C0_81E7_DB346DD67EE1_INCLUDED)
