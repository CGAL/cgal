//********************************************
// Viewport.cpp
//********************************************
// class CViewport
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
//********************************************

#include "Viewport.h"
#include <math.h>
#include <stdio.h>
#include <qgl.h>

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
// glDraw
//********************************************
void
CViewport::glDraw() const
{
  glViewport(originX,originY,sizeX,sizeY);
}

// ** EOF **
