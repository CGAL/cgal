/**
 * @file   smarter_ptr.h
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#ifndef SMARTER_PTR_H
#define SMARTER_PTR_H

#include <memory>

#include "config.h"
#ifdef DEBUG
    #include "util/SharedPtr.h"
    #define SHARED_PTR util::SharedPtr
#else
    #define SHARED_PTR std::shared_ptr
#endif

#include "util/WeakPtr.h"
#define WEAK_PTR util::WeakPtr

#endif /* SMARTER_PTR_H */
