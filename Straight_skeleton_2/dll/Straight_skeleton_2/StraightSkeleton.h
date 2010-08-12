#ifndef STRAIGHT_SKELETON_H
#define STRAIGHT_SKELETON_H

/*! \file StraightSkeleton.h
    \brief Functions for computing the straight skeleton of a simple polygon with holes.
    


This library provides functions to compute the straight skeleton of a polygon with holes.

<p>
The polygons must be simple and there must be no pairwise intersections. Note that 
this is not tested by the function. 


The references to pointers which are passed to the function will point
to memory allocated  by the dll. The user is in charge to free this memory
again. If allocation fails the pointer is set to NULL.
*/


typedef void (__stdcall *ProgressCallback) ( int aCurr, int aTotal ) ;


/**
       * Computes the straight skeleton of a simple polygon with holes.
       * @param np is the number of polygons.
       * @param np_i is an array of length np which holds the number of vertices of the i'th polygon.
       * @param xp is an array of length sum(np_i) which holds the x-coordinates of the polygon vertices.
       * @param yp is an array of length sum(np_i) which holds the y-coordinates of the polygon vertices.
       * @param numFaces is a reference parameter. After the call it holds the number of faces of the straight skeleton.
       * @param numVertices is a reference parameter. After the call it holds the number of vertices of all faces.
       * @param numFace_i is a reference parameter. After the call numFace_i[i] holds the number of vertices of the i'th face.
       * @param xf is a reference parameter. After the call it holds the x-ccordinates of the face vertices in an array of length sum(numFaces_i).
       * @param yf is a reference parameter. After the call it holds the y-ccordinates of the face vertices in an array of length sum(numFaces_i).
       * @param dumpEPS is a flag. When != 0 the skeleton is written to the file "dump.eps".
       * @return 0 if the straight skeleton COULD NOT be computed successfully.
       */	

extern "C" 
{

int __declspec (dllexport) StraightSkeleton( int np 
                                           , int* np_i
                                           , double* xp
                                           , double* yp
                                           , int& numFaces
                                           , int& numVertices
                                           , int*& numFace_i
                                           , double*& xf
                                           , double*& yf
                                           , int dumpEPS
                                           , ProgressCallback progress
                                           );


/**
       * Frees the memory allocated by the function StraightSkeleton. The pointers may be NULL, which can happen when the allocation failed.
       * @param numFace_i is a reference parameter. After the call it holds NULL.
       * @param xf is a reference parameter. After the call it holds NULL.
       * @param yf is a reference parameter. After the call it holds NULL.
       */	

void __declspec (dllexport) StraightSkeletonFree(int* numFace_i, double* xf, double* yf );

} // extern "C"

#endif


