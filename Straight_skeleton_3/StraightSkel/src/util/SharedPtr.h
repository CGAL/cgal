/**
 * @file   util/SharedPtr.h
 * @author Gernot Walzl
 * @date   2012-02-28
 */

#ifndef UTIL_SHAREDPTR_H
#define UTIL_SHAREDPTR_H

#include "util/StackTrace.h"
#include <iostream>
#include <memory>

namespace util {

template<class T> class SharedPtr : public std::shared_ptr<T> {
public:
    SharedPtr() : std::shared_ptr<T>() {}
    explicit SharedPtr(T* p) : std::shared_ptr<T>(p) {}
    template<class Y> SharedPtr(const SharedPtr<Y>& r) : std::shared_ptr<T>(r) {}
    SharedPtr(const std::shared_ptr<T>& r) : std::shared_ptr<T>(r) {}
    SharedPtr(const std::weak_ptr<T>& r) : std::shared_ptr<T>(r) {}

    T& operator*() const {  
        if (this->use_count() == 0) {
            std::cerr << "Error: Shared Pointer is invalid." << std::endl;
            StackTrace::print(std::cerr);
        }
        return std::shared_ptr<T>::operator*();
    }
    T* operator->() const {
        if (this->use_count() == 0) {
            std::cerr << "Error: Shared Pointer is invalid." << std::endl;
            StackTrace::print(std::cerr);
        }
        return std::shared_ptr<T>::operator->();
    }
};

}

#endif /* UTIL_SHAREDPTR_H */
