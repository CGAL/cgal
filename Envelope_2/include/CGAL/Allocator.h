// ======================================================================
//
// file          : Allocator.h
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// ======================================================================
#ifndef CGAL_ALLOCATOR_H
#define CGAL_ALLOCATOR_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

/*!
 * A memory allocator class which enables allocating objects in large
 * blocks and handles garbage collection.
 * The contained class has only to support the deafult constructor, the
 * assignment operator (operator=) and the destructor.
 */
template <class TYPE>
class Allocator
{
public:

    // Constants.
    enum
    {
        DefaultBlockSize = 2048,
        MinimalObjectsPerBlock = 100
    };

private:

    /*!
     * A block element.
     */
    struct Element
    {
        TYPE        object;         // The contained object.
        Element     *prevP;         // Points to the previous element 
                                    // in the element's linked list.
        Element     *nextP;         // Points to the next element 
                                    // in the element's linked list.
    };

    // Data members:

    int         iBlockSize;         // Size (in bytes) of each block.
    int         iObjectsPerBlock;   // Number of objects per block.
    int         iAllocatedBlocks;   // Number of allocated block pointers.
    int         iUsedBlocks;        // Number of used blocks.
    Element     **blocksPP;         // Pointers to the blocks 
                                    // (iAllocatedBlocks block pointers).
    Element     *firstFreeP;        // Points to the first free element.
    Element     *firstOccupiedP;    // Points to the first occupied element.

    // Copy constructor and assignment operator - not supported.
    Allocator (const Allocator& );
    const Allocator& operator= (const Allocator& );

public:

    /*!
     * Initialize an allocator with a given block size.
     * \param _iBlockSize The desired block size.
     */
    Allocator (const int& _iBlockSize = DefaultBlockSize);

    /*!
     * Destructor.
     */
    ~Allocator ();

    /*!
     * Get a pointer of a newly allocated object.
     * \return A pointer to a new object.
     */
    TYPE* Allocate ();

    /*!
     * Free the given object.
     * Notice that the destructor may not be called immediately.
     * \param objectP The object to dispose.
     * \pre the pointer objectP must be allocated by the Allocate() function.
     */
    void Free (TYPE* objectP);

    /*!
     * Free all objects allocated by the allocator.
     */
    void Reset ();

protected:

    /*!
     * Allocate and initialize a single block.
     */
    void _AllocateBlock ();
};

//---------------------------------------------------------
template <class TYPE>
Allocator<TYPE>::Allocator (const int& _iBlockSize) :
    iBlockSize(_iBlockSize),
    iObjectsPerBlock(iBlockSize / sizeof(Element)),
    iAllocatedBlocks(1),
    iUsedBlocks(0),
    blocksPP(NULL),
    firstFreeP(NULL),
    firstOccupiedP(NULL)
{
    // Allocate the block pointers.
    blocksPP = new Element* [1];

    // Set the number of objects per block.
    if (iObjectsPerBlock < MinimalObjectsPerBlock)
        iObjectsPerBlock = MinimalObjectsPerBlock;

    // Allocate the first block.
    blocksPP[0] = NULL;
    _AllocateBlock();
}

//---------------------------------------------------------
template <class TYPE>
Allocator<TYPE>::~Allocator ()
{
    if (blocksPP != NULL)
    {
        // Go over all blocks and free them.
        for (int i = 0; i < iAllocatedBlocks; i++)
        {
            if (blocksPP[i] != NULL)
                delete[] (blocksPP[i]);
        }

        // Free the block pointers.
        delete blocksPP;
    }
    blocksPP = NULL;
}

//---------------------------------------------------------
template <class TYPE>
TYPE* Allocator<TYPE>::Allocate ()
{
    // In case there are no free elements left, allocate a new block of
    // elements.
    if (firstFreeP == NULL)
        _AllocateBlock();
    
    // Get a pointer to the first free element.
    Element     *elemP = firstFreeP;

    // Remove the element from the free elements' list.
    if (firstFreeP->nextP != NULL)
    {
        // Make the next free element the new free list's head.
        firstFreeP->nextP->prevP = NULL;
        firstFreeP = firstFreeP->nextP;
    }
    else
    {
        // In this case, there are no more free elements.
        firstFreeP = NULL;
    }

    // Add the element to the head of occupied elements' list.
    if (firstOccupiedP != NULL)
    {
        // Push the element before the current list head.
        firstOccupiedP->prevP = elemP;
        elemP->prevP = NULL;
        elemP->nextP = firstOccupiedP;
        firstOccupiedP = elemP;
    }
    else
    {
        // In case this is the only occupied element:
        elemP->prevP = NULL;
        elemP->nextP = NULL;
        firstOccupiedP = elemP;
    }

    // Return a pointer to the contained object.
    return (&(elemP->object));
}

//---------------------------------------------------------
// Free the given object.
//
template <class TYPE>
void Allocator<TYPE>::Free (TYPE* objectP)
{
    // Translate the object pointer to a block element pointer.
    // Note that this is quite risky in general, but we rely here on the fact
    // that the pointer was allocated by the Allocate() function.
    Element     *elemP = reinterpret_cast<Element*>(objectP);

    // Delete the element from the occupied elements' list.
    if (elemP->prevP != NULL)
    {
        // Connect the previous element directly to the next one.
        elemP->prevP->nextP = elemP->nextP;
    }
    else
    {
        // elemP is currently the list head: Assign a new list head.
        firstOccupiedP = elemP->nextP;
    }

    if (elemP->nextP != NULL)
        // Connect the next element directly to the previous one.
        elemP->nextP->prevP = elemP->prevP;

    // Add the element to the head of free elements' list.
    if (firstFreeP != NULL)
    {
        // Push the element before the current list head.
        firstFreeP->prevP = elemP;
        elemP->prevP = NULL;
        elemP->nextP = firstFreeP;
        firstFreeP = elemP;
    }
    else
    {
        // In case this is the only free element:
        elemP->prevP = NULL;
        elemP->nextP = NULL;
        firstFreeP = elemP;
    }

    return;
}

//---------------------------------------------------------
// Free all objects allocated by the allocator.
//
template <class TYPE>
void Allocator<TYPE>::Reset ()
{
    // Free all blocks.
    if (blocksPP != NULL)
    {
        // Go over all blocks and free them.
        for (int i = 0; i < iAllocatedBlocks; i++)
        {
            if (blocksPP[i] != NULL)
                delete[] (blocksPP[i]);
        }

        // Free the block pointers.
        delete blocksPP;
    }

    // Allocate the block pointers.
    blocksPP = new Element* [1];
    iAllocatedBlocks = 1;
    iUsedBlocks = 0;
    firstFreeP = NULL;
    firstOccupiedP = NULL;

    // Allocate the first block.
    blocksPP[0] = NULL;
    _AllocateBlock();

    return;
}

//---------------------------------------------------------
// Allocate and initialize a single block.
//
template <class TYPE>
void Allocator<TYPE>::_AllocateBlock ()
{
    int             i;

    // Check if enough block pointers are allocated.
    // If necessary, double the number of block pointers.
    if (iUsedBlocks == iAllocatedBlocks)
    {
        Element     **newBlocksPP = new Element* [2*iAllocatedBlocks];

        // Copy the used block pointers.
        for (i = 0; i < iAllocatedBlocks; i++)
            newBlocksPP[i] = blocksPP[i];

        // Nullify the rest of the pointers.
        for (i = iAllocatedBlocks; i < 2*iAllocatedBlocks; i++)
            newBlocksPP[i] = NULL;

        // Update the block pointers.
        delete[] blocksPP;
        blocksPP = newBlocksPP;

        iAllocatedBlocks *= 2;
    }

    // Allocate a new block and link all its elements one after the other.
    Element     *blockP = new Element [iObjectsPerBlock];
    Element     *elemP = blockP;

    for (i = 0; i < iObjectsPerBlock; i++)
    {
        elemP->prevP = elemP - 1;
        elemP->nextP = elemP + 1;

        elemP++;
    }

    // This function should be called only when the free elements list is
    // empty, so we assign the list we've just created as the free elements'
    // list.
    firstFreeP = blockP;
    blockP->prevP = NULL;
    (blockP + iObjectsPerBlock - 1)->nextP = NULL;

    // Assign the new block pointer.
    blocksPP[iUsedBlocks] = blockP;
    iUsedBlocks++;

    return;
}

CGAL_END_NAMESPACE

#endif
