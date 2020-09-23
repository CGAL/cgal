// Copyright (c) 2013  The University of Western Sydney, Australia.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Quincy Tse, Weisheng Si

/* Plane_scan_tree_impl.h
 *
 * This header implements the details of the class Plane_scan_tree.
 */

#ifndef CGAL_CONE_SPANNER_PLANE_SCAN_TREE_IMPL_H
#define CGAL_CONE_SPANNER_PLANE_SCAN_TREE_IMPL_H

#include <CGAL/license/Cone_spanners_2.h>


#include <stdexcept>

namespace CGAL {

/* This namespace contains the internal implementation of the tree for building Theta graph.
 * This is not meant to be imported by other codes.
 */
namespace ThetaDetail {

/* Contains the internal data structure for the Plane_scan_tree.
 *
 * Plane_scan_tree is a ternary B+ tree using linked structures.
 */
template < typename Key, typename T, typename Comp, typename VComp >
class Plane_scan_tree;

template < typename Key, typename T, typename Comp=std::less<Key>, typename VComp=std::less<const T> >
class _Node;
template < typename Key, typename T, typename Comp=std::less<Key>, typename VComp=std::less<const T> >
class _Leaf;
template < typename Key, typename T, typename Comp=std::less<Key>, typename VComp=std::less<const T> >
class _Internal;

template < typename Key, typename T, typename Comp=std::less<Key>, typename VComp=std::less<const T> >
class _Iterator;

template < typename Key, typename T, typename Comp=std::less<Key>, typename VComp=std::less<const T> >
class _RIterator;

/* Abstract superclass */
template < typename Key, typename T, typename Comp, typename VComp >
class _Node {
public:
    typedef Key                                 key_type;
    typedef T                                   mapped_type;
    typedef std::pair<const Key, T>             value_type;
    typedef Comp                                key_compare;
    typedef VComp                               value_compare;

    typedef _Node<Key, T, Comp, VComp>          _node_type;
    typedef _Leaf<Key, T, Comp, VComp>          _leaf_type;
    typedef _Internal<Key, T, Comp, VComp>      _internal_type;
    typedef Plane_scan_tree<Key, T, Comp, VComp>  tree_type;

    /* Constructor */
    _Node(const key_compare& less, const value_compare& vless, tree_type *const t)
        : parent(nullptr), less(less), vless(vless), tree(t) {}

    /* Destructor */
    virtual ~_Node() {}

    /* @return true if the node is a leaf node, false otherwise. */
    virtual bool isLeaf() const = 0;

    /* Retrieves the leaf node that handles the specified key.
     *
     * @param k The key being queried
     * @return The leaf node that handles the key.
     */
    virtual _leaf_type* leafNode(const key_type& k) = 0;

    /* Returns the minimum value whose key is greater than x.
     *
     * @param x The threshold key
     * @return the minimum value whose key is greater than x.
     */
    virtual const mapped_type* minAbove(const key_type& x) const = 0;

    friend std::ostream& operator<< (std::ostream& os, const _node_type& n) {
        n.print(os, 0);
        return os;
    }

protected:
    /* Prints a string representation of the node in the ostream. */
    virtual void print(std::ostream&, const size_t) const = 0;

    /* Sets the parent internal node for the node */
    void setParent(_internal_type *const p) {
        parent = p;
    }

    /* Returns the minimum value in the subtree rooted at the current node. */
    virtual const mapped_type* minV() const = 0;

    _internal_type* parent;
    const key_compare& less;
    const value_compare& vless;
    tree_type* tree;

    friend class _Leaf<key_type, mapped_type, key_compare, value_compare>;
    friend class _Internal<key_type, mapped_type, key_compare, value_compare>;
};

/* Leaf node of the B+ tree. */
template < typename Key, typename T, typename Comp, typename VComp >
class _Leaf : public _Node<Key, T, Comp, VComp> {
public:
    typedef _Node<Key, T, Comp, VComp>  _node_type;

    typedef typename _node_type::key_type       key_type;
    typedef typename _node_type::mapped_type    mapped_type;
    typedef typename _node_type::value_type     value_type;
    typedef typename _node_type::key_compare    key_compare;
    typedef typename _node_type::value_compare  value_compare;

    typedef typename _node_type::_leaf_type     _leaf_type;
    typedef typename _node_type::_internal_type _internal_type;
    typedef typename _node_type::tree_type      tree_type;

    _Leaf (const key_compare& less, const value_compare& vless, tree_type *const t,
           _leaf_type *const prev = nullptr,
           _leaf_type *const next = nullptr)
        : _node_type (less, vless, t), prev (prev), next (next) {
        std::memset (values, 0, 2*sizeof(value_type*));
    }

    /* Destructor.
     * Frees memory used for storing key-value pair, thus invalidating any
     * exisitng pointers to any keys and/or values in the tree. During and
     * after destruction, neighbour nodes are not guarenteed to be consistent.
     * Specifically, the linked list along the leaves of the B+ tree is
     * invalidated. */
    virtual ~_Leaf() {
        delete values[0];
        delete values[1];
        values[0] = nullptr;
        values[1] = nullptr;
        prev = nullptr;
        next = nullptr;
    }

    virtual bool isLeaf() const {
        return true;
    }

    virtual _leaf_type* leafNode(const key_type&) {
        return this;
    }

    void add(const key_type& k, const mapped_type& v) {
        if (nullptr == values[0]) {
            // empty
            values[0] = new value_type(k, v);
        } else if (nullptr == values[1]) {
            // Not full;
            if (this->less(k, values[0]->first)) {
                values[1] = values[0];
                values[0] = new value_type(k, v);
                if (this->parent && this->vless(v, values[1]->second))
                    this->parent->updateMin(this);
            } else {
                values[1] = new value_type(k, v);
                if (this->parent && this->vless(v, values[0]->second))
                    this->parent->updateMin(this);
            }
        } else {
            _leaf_type* split
                = new _leaf_type(this->less, this->vless, this->tree, this, next);
            if (nullptr != next) next->prev = split;
            next = split;
            if (this->less(k, values[0]->first)) {
                // k, [0], [1]
                split->values[0] = values[0];
                split->values[1] = values[1];
                values[0] = new value_type(k, v);
                values[1] = nullptr;
            } else if (this->less(k, values[1]->first)) {
                // [0], k, [1]
                split->values[0] = new value_type(k, v);
                split->values[1] = values[1];
                values[1] = nullptr;
            } else {
                split->values[0] = values[1];
                split->values[1] = new value_type(k, v);
                values[1] = nullptr;
            }

            // Update pointer to max leaf (for reverse iterator)
            if (this->tree->m_max == this) this->tree->m_max = split;

            // Create new parent node current node is not root
            if (nullptr == this->parent) {
                this->setParent(new _internal_type(this->less, this->vless, this->tree));
                this->tree->root = this->parent;
            }
            // Promote middle key and get parent to readjust pointers
            this->parent->splitMe (&(split->values[0]->first), this, split);
        }
    }

    virtual const mapped_type* minAbove(const key_type& x) const {
        if ( !this->less(x, values[0]->first) && !this->less(values[0]->first, x) // equals
                && nullptr != values[1]) {
            return &values[1]->second;
        }
        return nullptr;
    }

    virtual const mapped_type* minV() const {
        return (nullptr == values[1]) ? &values[0]->second : &(std::min)(values[0]->second, values[1]->second, this->vless);
    }

protected:
    virtual void print(std::ostream& os, const size_t /* level */) const {
        os << "\t\"" << this << "\"--\"" << &(values[0]->first) << "\" [style=bold];" << std::endl;
        os << "\t" << "{rank=same; \"" << &(values[0]->first) << "\"--\"" << &(values[0]->second) << "\" [style=dotted];}" << std::endl;
        os << "\t\"" << &(values[0]->first) << "\"--\"" << values[0]->first << "\";" << std::endl;
        os << "\t\"" << &(values[0]->second) << "\"--\"" << values[0]->second << "\";" << std::endl;
        if (nullptr != values[1]) {
            os << "\t\"" << this << "\"--\"" << &(values[1]->first) << "\" [style=bold];" << std::endl;
            os << "\t" << "{rank=same; \"" << &(values[1]->first) << "\"--\"" << &(values[1]->second) << "\" [style=dotted];}" << std::endl;
            os << "\t" << "{rank=same; \"" << &(values[0]->second) << "\"--\"" << &(values[1]->first) << "\" [color=white]; rankdir=LR;}" << std::endl;
            os << "\t\"" << &(values[1]->first) << "\"--\"" << values[1]->first << "\";" << std::endl;
            os << "\t\"" << &(values[1]->second) << "\"--\"" << values[1]->second << "\";" << std::endl;
        }
        os << "\t\"" << this << "\" [style=diagonals];" << std::endl;
    }

private:
    /* Key-value pairs */
    value_type* values[2];

    /* Linked list structure of the B+ tree */
    _leaf_type* prev;
    _leaf_type* next;

    friend class _Iterator<key_type, mapped_type, key_compare, value_compare>;
    friend class _RIterator<key_type, mapped_type, key_compare, value_compare>;
};

template < typename Key, typename T, typename Comp, typename VComp >
class _Internal : public _Node<Key, T, Comp, VComp> {
public:
    typedef _Node<Key, T, Comp, VComp>  _node_type;

    typedef typename _node_type::key_type       key_type;
    typedef typename _node_type::mapped_type    mapped_type;
    typedef typename _node_type::value_type     value_type;
    typedef typename _node_type::key_compare    key_compare;
    typedef typename _node_type::value_compare  value_compare;

    typedef typename _node_type::_leaf_type     _leaf_type;
    typedef typename _node_type::_internal_type _internal_type;
    typedef typename _node_type::tree_type      tree_type;

    _Internal (const Comp& less, const VComp& vless, tree_type *const t)
        : _node_type(less, vless, t) {
        std::memset (keys, 0, 2*sizeof(key_type*));
        std::memset (children, 0, 3*sizeof(_node_type*));
        std::memset (vMin, 0, 3*sizeof(mapped_type*));
    }

    virtual ~_Internal() {
        keys[0] = nullptr;
        keys[1] = nullptr;

        delete children[0];
        children[0] = nullptr;
        delete children[1];
        children[1] = nullptr;
        delete children[2];
        children[2] = nullptr;

        vMin[0] = nullptr;
        vMin[1] = nullptr;
        vMin[2] = nullptr;
    }

    virtual bool isLeaf() const {
        return false;
    }

    virtual _leaf_type* leafNode(const key_type& k) {
        int i = 0;
        for (i = 0; (i < 2) && (nullptr != keys[i]); i++) {
            if (this->less(k, *(keys[i]))) return children[i]->leafNode(k);
        }
        return children[i]->leafNode(k);
    }

    /* Update references to min values.
     *
     *  @param ch   The child node involved.
     */
    void updateMin(const _node_type *const ch) {
        int i = 0;
        while (i < 3) {
            if (children[i] == ch) break;
            ++i;
        }
        if (i >= 3)
            throw std::runtime_error("Cannot find child");

        vMin[i] = ch->minV();
        if (this->parent && minV() == vMin[i])
            this->parent->updateMin(this);
    }

    /* Process the splitting of children. This is called by the affected
     * child AFTER splitting had occurred. This function updates relevant
     * pointers, and splits current internal node as necessary.
     *
     * @param k     The key the child was split on.
     * @param left  The (existing) child that lies to the left of k.
     * @param right The (new) child that lies to the right of k.
     */
    void splitMe(const key_type *const k,
                 _node_type *const left,
                 _node_type *const right) {
        if (keys[0] == nullptr) {
            // New root
            keys[0] = k;
            children[0] = left;
            left->setParent(this);
            vMin[0] = left->minV();
            children[1] = right;
            right->setParent(this);
            vMin[1] = right->minV();
        } else if (keys[1] == nullptr) {
            // not full
            if (left == children[0]) {
                // split 0th
                keys[1] = keys[0];
                keys[0] = k;
                children[2] = children[1];
                children[1] = right;
                right->setParent(this);

                vMin[0] = children[0]->minV();
                vMin[1] = children[1]->minV();
                vMin[2] = children[2]->minV();
            } else {
                // split 1st
                keys[1] = k;
                children[2] = right;
                right->setParent(this);

                // did not touch children[0]
                vMin[1] = children[1]->minV();
                vMin[2] = children[2]->minV();
            }
            if (this->parent)
                this->parent->updateMin(this);
        } else {
            // full
            _internal_type* split
                = new _internal_type(this->less, this->vless, this->tree);
            const key_type * toPromote;
            if (left == children[0]) {
                // Split 0th
                split->keys[0] = keys[1];
                toPromote = keys[0];
                keys[1] = nullptr;
                keys[0] = k;

                split->children[1] = children[2];
                split->vMin[1] = split->children[1]->minV();
                split->children[1]->setParent(split);

                split->children[0] = children[1];
                split->vMin[0] = split->children[0]->minV();
                split->children[0]->setParent(split);

                children[1] = right;
                vMin[1] = children[1]->minV();
                children[1]->setParent(this);

                vMin[0] = children[0]->minV();

                children[2] = nullptr;
                vMin[2] = nullptr;
            } else if (left == children[1]) {
                // Split 1st
                split->keys[0] = keys[1];
                toPromote = k;
                keys[1] = nullptr;

                split->children[1] = children[2];
                split->vMin[1] = split->children[1]->minV();
                split->children[1]->setParent(split);

                split->children[0] = right;
                split->vMin[0] = split->children[0]->minV();
                split->children[0]->setParent(split);

                vMin[1] = children[1]->minV();

                children[2] = nullptr;
                vMin[2] = nullptr;
            } else {
                // Split last
                split->keys[0] = k;
                toPromote = keys[1];
                keys[1] = nullptr;

                split->children[1] = right;
                split->vMin[1] = split->children[1]->minV();
                split->children[1]->setParent(split);

                split->children[0] = left;
                split->vMin[0] = split->children[0]->minV();
                split->children[0]->setParent(split);

                children[2] = nullptr;
                vMin[2] = nullptr;
            }

            if (nullptr == this->parent) {
                this->setParent(new _internal_type(this->less, this->vless, this->tree));
                this->tree->root = this->parent;
            }

            this->parent->splitMe (toPromote, this, split);
        }
    }

    virtual const mapped_type* minAbove(const key_type& x) const {
        if (this->less(x, *keys[0])) {
            // x in left tree
            const mapped_type* minFromCh = children[0]->minAbove(x);
            const mapped_type* res =
                (nullptr == minFromCh) ?
                vMin[1] :
                &(std::min) (*children[0]->minAbove(x), *vMin[1], this->vless);
            if (vMin[2] != nullptr)
                res = &(std::min) (*res, *vMin[2], this->vless);
            return res;
        } else if (nullptr == keys[1] || this->less(x, *keys[1])) {
            // x in middle
            const mapped_type* res = children[1]->minAbove(x);
            if (nullptr == res) return vMin[2];

            if (vMin[2] != nullptr)
                res = &(std::min) (*res, *vMin[2], this->vless);
            return res;
        } else {
            return children[2]->minAbove(x);
        }
    }

    virtual const mapped_type* minV() const {
        const mapped_type* res = &(std::min)(*vMin[0], *vMin[1], this->vless);
        if (nullptr != children[2])
            res = &(std::min)(*res, *vMin[2], this->vless);
        return res;
    }

protected:
    virtual void print(std::ostream& os, const size_t level) const {
        os << "\t\"" << this << "\"--\"" << children[0] << "\";" << std::endl;
        os << "\t\"" << this << "\"--\"" << children[1] << "\";" << std::endl;
        if (nullptr != children[2])
            os << "\t\"" << this << "\"--\"" << children[2] << "\";" << std::endl;

        children[0]->print(os, level+1);
        children[1]->print(os, level+1);
        if (nullptr != children[2])
            children[2]->print(os, level+1);

        os << "\t\"" << this << "\"--\"" << vMin[0] << "\" [style=dashed,label=vMin0];" << std::endl;
        os << "\t\"" << this << "\"--\"" << vMin[1] << "\" [style=dashed,label=vMin1];" << std::endl;
        if (nullptr != vMin[2])
            os << "\t\"" << this << "\"--\"" << vMin[2] << "\" [style=dashed,label=vMin2];" << std::endl;

        os << "\t\"" << this << "\"--\"" << keys[0] << "\" [style=dotted,label=keys0];" << std::endl;
        if (nullptr != keys[1])
            os << "\t\"" << this << "\"--\"" << keys[1] << "\" [style=dotted,label=keys1];" << std::endl;
    }

private:
    const key_type* keys[2];
    _node_type* children[3];
    const mapped_type* vMin[3];
};

template <typename Key, typename T, typename Comp, typename VComp >
class _Iterator {
public:
    typedef           _Leaf<Key, T, Comp, VComp>      leaf_type;

    typedef typename  leaf_type::value_type           value_type;
    typedef           _Iterator<Key, T, Comp, VComp>  iterator_type;

    _Iterator(leaf_type* start = nullptr) : cell (0), leaf(start) {}
    _Iterator(leaf_type* start, const Key& key) : cell(0), leaf(start) {
        if (start->values[0]->first != key) {
            if (!start->values[1] || start->values[1]->first != key)
                leaf = nullptr;
            else
                cell = 1;
        }
    }
    ~_Iterator() {
        cell = 0;
        leaf = nullptr;
    }

    value_type& operator*() const {
        return *(leaf->values[cell]);
    }

    value_type* operator->() const {
        return leaf->values[cell];
    }

    iterator_type& operator++() {
        cell++;
        if (nullptr == leaf) return *this;

        if (cell > 1 || nullptr == leaf->values[cell]) {
            cell = 0;
            leaf = leaf->next;
        }
        return *this;
    }

    iterator_type operator++(int) {
        iterator_type res (*this);
        ++(*this);
        return res;
    }

    bool operator==(const iterator_type& x) const {
        return (leaf == x.leaf) && (cell == x.cell);
    }

    bool operator!=(const iterator_type& x) const {
        return (leaf != x.leaf) || (cell != x.cell);
    }

private:
    size_t cell;
    leaf_type* leaf;

    friend class _Node<Key, T, Comp, VComp>;
};

template <typename Key, typename T, typename Comp, typename VComp >
class _RIterator {
public:
    typedef           _Leaf<Key, T, Comp, VComp>      leaf_type;

    typedef typename  leaf_type::value_type           value_type;
    typedef           _RIterator<Key, T, Comp, VComp> iterator_type;

    _RIterator(_Leaf<Key, T, Comp, VComp>* start = nullptr) : cell (1), leaf(start) {
        if (nullptr == start) return;

        if (nullptr == leaf->values[cell]) {
            cell--;
            if (nullptr == leaf->values[cell]) leaf = nullptr;
        }
    }

    ~_RIterator() {
        cell = 0;
        leaf = nullptr;
    }

    value_type& operator*() const {
        return *(leaf->values[cell]);
    }

    value_type* operator->() const {
        return leaf->values[cell];
    }

    iterator_type& operator++() {
        if (nullptr == leaf) return *this;
        if (cell == 1) {
            cell--;
        } else if (cell == 0) {
            cell = 1;
            leaf = leaf->prev;
            if (nullptr == leaf) return *this;

            if (nullptr == leaf->values[cell]) cell--;
        }
        return *this;
    }

    iterator_type operator++(int) {
        iterator_type res (*this);
        ++(*this);
        return res;
    }

    bool operator==(const iterator_type& x) const {
        return (leaf == x.leaf) && (cell == x.cell);
    }

    bool operator!=(const iterator_type& x) const {
        return (leaf != x.leaf) || (cell != x.cell);
    }

private:
    size_t cell;
    leaf_type* leaf;

    friend class _Node<Key, T, Comp, VComp>;
};

}    // namespace ThetaDetail

}   // namespace CGAL

#endif
