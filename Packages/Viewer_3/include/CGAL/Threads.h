// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Threads.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
#include <pthread.h>
#include <iostream.h>


class Synchronizer
{
public: 
 static void initAll();
 
 //to initialize the viewer
 static pthread_mutex_t InitMutex;
  
 //to protect the scene graph
 static pthread_mutex_t sgMutex;
 
 //to return to the debugger  
 static pthread_mutex_t sMutex;
 static pthread_cond_t  sCond;
};


void sendSignal();

void stop();


pthread_mutex_t Synchronizer::InitMutex;

pthread_mutex_t Synchronizer::sgMutex;
pthread_mutex_t Synchronizer::sMutex;
pthread_cond_t  Synchronizer::sCond;

void Synchronizer::initAll()
{
 pthread_mutex_init(&Synchronizer::InitMutex, NULL);
 pthread_mutex_init(&Synchronizer::sgMutex, NULL);
 pthread_mutex_init(&Synchronizer::sMutex, NULL);
 pthread_cond_init(&Synchronizer::sCond, NULL);
}


void sendSignal()
{
 cerr << "condition...";
 pthread_mutex_lock(&Synchronizer::sMutex);
 pthread_cond_signal(&Synchronizer::sCond);
 pthread_mutex_unlock(&Synchronizer::sMutex); 
 cerr << "...signaled\n";
}


void stop()
{
 cerr << "stop...";
 pthread_mutex_lock(&Synchronizer::sMutex);
 pthread_cond_wait(&Synchronizer::sCond, &Synchronizer::sMutex);
 pthread_mutex_unlock(&Synchronizer::sMutex);
 cerr << "...executed\n";
} 

