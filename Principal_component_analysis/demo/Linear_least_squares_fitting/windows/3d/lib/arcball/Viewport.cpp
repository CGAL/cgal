//********************************************
// Viewport.cpp
//********************************************
// class CViewport
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
//********************************************

#include "stdafx.h"
#include "Viewport.h"

#include <math.h>
#include <stdio.h>


//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////


//********************************************
// Set
//********************************************
void
CViewport::Set(const CViewport& v)
{
  sizeX = v.sizeX;
  sizeY = v.sizeY;
  originX = v.originX;
  originY = v.originY;
}
//********************************************
// Set
//********************************************
void
CViewport::Set(const CViewport *pV)
{
  sizeX = pV->sizeX;
  sizeY = pV->sizeY;
  originX = pV->originX;
  originY = pV->originY;
}

//********************************************
// Trace
//********************************************
void
CViewport::Trace() const
{
  TRACE("\n");
  TRACE("** Viewport **\n");
  TRACE("Address      : %x\n",this);
  TRACE("Image Origin : (%i %i)\n",originX,originY);
  TRACE("Image Size   : (%i %i)\n",sizeX,  sizeY);
}

//********************************************
// glDraw
//********************************************
void
CViewport::glDraw() const
{
  glViewport(originX,originY,sizeX,sizeY);
}

// ** EOF **
