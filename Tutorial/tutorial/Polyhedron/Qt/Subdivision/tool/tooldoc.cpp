/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// include files for Qt
#include <qdir.h>
#include <qfileinfo.h>
#include <qwidget.h>
#include <qmsgbox.h>
#include <qfiledialog.h>
#include <qlineedit.h>
#include <qtextstream.h>
#include <fstream>
#include <iostream>

using namespace std;

// application specific includes
#include "tooldoc.h"
#include "tool.h"
#include "toolview.h"

// cgal
#include <CGAL/basic.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "quad-triangle.h"
#include "sqrt3.h"
#include "SurfLab/Polyhedron_subdivision.h"

//***************************************
// life cycle
//***************************************
ToolDoc::ToolDoc()
{
  pViewList = new QList<ToolView>;
  pViewList->setAutoDelete(false);
}
ToolDoc::~ToolDoc()
{
  delete pViewList;
  delete m_pMesh;
  pViewList = NULL; // fix
}

void ToolDoc::addView(ToolView *view)
{
  pViewList->append(view);
  changedViewList();
}

void ToolDoc::removeView(ToolView *view)
{
  pViewList->remove(view);
  if(!pViewList->isEmpty())
    changedViewList();
  else
    deleteContents();
}

void ToolDoc::changedViewList(){

  ToolView *w;
  if((int)pViewList->count() == 1){
    w=pViewList->first();
    w->setCaption(m_title);
  }
  else{
    int i;
    for( i=1,w=pViewList->first(); w!=0; i++, w=pViewList->next())
      w->setCaption(QString(m_title+":%1").arg(i));
  }
}

bool ToolDoc::isLastView() {
  return ((int) pViewList->count() == 1);
}


void ToolDoc::updateAllViews(ToolView *sender)
{
  ToolView *w = NULL;
  for(w  = pViewList->first();
      w != 0;
      w  = pViewList->next())
    w->update(sender);
}

void ToolDoc::setPathName(const QString &name)
{
  m_filename=name;
  m_title=QFileInfo(name).fileName();
}

const QString& ToolDoc::pathName() const
{
  return m_filename;
}

void ToolDoc::setTitle(const QString &title)
{
  m_title=title;
}

const QString &ToolDoc::title() const
{
  return m_title;
}


void ToolDoc::closeDocument()
{
  ToolView *w;
  if(!isLastView())
  {
    for(w=pViewList->first(); w!=0; w=pViewList->next())
    {
        if(!w->close())
         break;
    }
  }
  if(isLastView())
  {
    w=pViewList->first();
    w->close();
  }
}

bool ToolDoc::newDocument()
{
  modified=false;
  return true;
}

//************************************************
// open document
//************************************************
bool ToolDoc::openDocument(const QString &filename,
                           const char *format)
{
  // OFF extension
  QFile f(filename);
  if(filename.findRev(".off") != -1)
  {
    // check for valid stream
    fprintf(stderr,"\n");
    fprintf(stderr,"attempt to open file %s...",filename.latin1());
    std::ifstream stream(filename.latin1());
    if(!stream)
    {
      fprintf(stderr,"failed\n");
      return false;
    }


    // eat the mesh
  	// alloc a new mesh
    m_pMesh = new Enriched_polyhedron<My_kernel,Enriched_items>;
    fprintf(stderr,"ok\nfill mesh...");
    stream >> *m_pMesh;

    // print mesh info
    fprintf(stderr,"(%d faces, ",m_pMesh->size_of_facets());
    fprintf(stderr,"%d vertices)\n",m_pMesh->size_of_vertices());

    // compute normals
    m_pMesh->compute_normals();
  }
  else
    fprintf(stderr,"unknown extension\n");

  // update document
  modified = false;
  m_filename = filename;
  m_title = QFileInfo(f).fileName();

  return true;
}

bool ToolDoc::saveDocument(const QString &filename,
                           const char *format /*=0*/)
{
  QFile f(filename);
  std::ofstream stream(filename.latin1());
  if(!stream)
  {
    fprintf(stderr,"save failed\n");
    return false;
  }

  // OFF
  if(filename.findRev(QString(".off"),-1,FALSE) != -1)
  {
    // push the mesh in a .off file
    fprintf(stderr,"saving to %s...",filename.latin1());
    stream << *m_pMesh;
    fprintf(stderr,"ok");
    modified = false;
    m_filename = filename;
    m_title = QFileInfo(f).fileName();
    return true;
  }

  return false;
}

void ToolDoc::deleteContents()
{
}

bool ToolDoc::canCloseFrame(ToolView* pFrame)
{
  if(!isLastView())
    return true;

  bool ret=false;
  if(isModified())
  {
    QString saveName;
    switch(QMessageBox::information(pFrame, title(), tr("The current file has been modified.\n"
                          "Do you want to save it?"),QMessageBox::Yes, QMessageBox::No, QMessageBox::Cancel ))
    {
      case QMessageBox::Yes:
        if(title().contains(tr("Untitled")))
        {
          saveName=QFileDialog::getSaveFileName(0, 0, pFrame);
                                          if(saveName.isEmpty())
                                            return false;
        }
        else
          saveName=pathName();

        if(!saveDocument(saveName))
        {
           switch(QMessageBox::critical(pFrame, tr("I/O Error !"), tr("Could not save the current document !\n"
                                                        "Close anyway ?"),QMessageBox::Yes ,QMessageBox::No))

           {
             case QMessageBox::Yes:
               ret=true;
             case QMessageBox::No:
               ret=false;
           }
        }
        else
          ret=true;
        break;
      case QMessageBox::No:
        ret=true;
        break;
      case QMessageBox::Cancel:
      default:
        ret=false;
        break;
    }
  }
  else
    ret=true;

  return ret;
}


//***********************************************
// drawScene
//***********************************************
void ToolDoc::drawScene(bool move,
                        bool superimposededges,
                        bool superimposedvertices,
                        bool first,
                        bool smooth,
                        bool use_normals,
                        float r, float g, float b)
{
  glColor3f(r, g, b); //set the foreground color
  if(superimposededges || superimposedvertices)
  {
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(2.0,3.0);
  }
  m_pMesh->gl_draw(smooth, use_normals);
  
  if(superimposededges || superimposedvertices)
    ::glDisable(GL_POLYGON_OFFSET_FILL);

  glColor3f(0.0f,0.0f,0.0f);
  if(superimposededges)
    m_pMesh->superimpose_edges();

  glColor3f(1.0f,1.0f,0.0f);
  if(superimposedvertices)
    m_pMesh->superimpose_spheres(0.1);
  glColor3f(r, g, b); //restore the color
}

//***********************************************
// quad-triangle subdivision
//***********************************************
void ToolDoc::subdivisionStamLoop()
{
	// subdivision engine
	CSubdivider_quad_triangle<Enriched_polyhedron<My_kernel,Enriched_items>,My_kernel> subdivider;

	// alloc a new mesh
  Enriched_polyhedron<My_kernel,Enriched_items> *pNewMesh =
		new Enriched_polyhedron<My_kernel,Enriched_items>;

	// subdivide once
	subdivider.subdivide(*m_pMesh,*pNewMesh,true);

	// copy bounding box (approximate, but fast)
	pNewMesh->copy_bounding_box(m_pMesh);

	// delete previous mesh
	delete m_pMesh;

	// set new mesh
	m_pMesh = pNewMesh;

	// compute normals
	m_pMesh->compute_normals();
}

//***********************************************
// sqrt3 subdivision
//***********************************************
void ToolDoc::subdivisionSqrt3()
{
	// subdivision engine
	CSubdivider_sqrt3<Enriched_polyhedron<My_kernel,Enriched_items>,My_kernel> subdivider;
	subdivider.subdivide(*m_pMesh,1); // one iteration
	m_pMesh->compute_normals();

	// compute normals
	m_pMesh->compute_normals();
}

//***********************************************
// sqrt3 subdivision
//***********************************************
void ToolDoc::subdivisionSqrt3Twice()
{
	// subdivision engine
	CSubdivider_sqrt3<Enriched_polyhedron<My_kernel,Enriched_items>,My_kernel> subdivider;
	subdivider.subdivide(*m_pMesh,2); // two iterations
	m_pMesh->compute_normals();

	// compute normals
	m_pMesh->compute_normals();
}

//***********************************************
// Doo-Sabin subdivision
//***********************************************
void ToolDoc::subdivisionDooSabin()
{
	// subdivision engine
  Polyhedron_subdivision<Enriched_polyhedron<My_kernel,Enriched_items> >::DooSabin_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
}

//***********************************************
// Loop subdivision
//***********************************************
void ToolDoc::subdivisionLoop()
{
	// subdivision engine
  Polyhedron_subdivision<Enriched_polyhedron<My_kernel,Enriched_items> >::Loop_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
}

//***********************************************
// Catmull-Clark subdivision
//***********************************************
void ToolDoc::subdivisionCatmullClark()
{
	// subdivision engine
  Polyhedron_subdivision<Enriched_polyhedron<My_kernel,Enriched_items> >::CatmullClark_subdivision(*m_pMesh,1);

	// compute normals
	m_pMesh->compute_normals();
}





#include "tooldoc.moc"

