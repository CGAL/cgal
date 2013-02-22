#include <boost/config.hpp>

#if defined(BOOST_MSVC)
 // Avoid warning concerning spatial_sort(QList::begin(), QList.end() QT "bug" 
#  pragma warning(disable: 4267 )
#  pragma warning(disable: 4244 )
#endif

#include "Viewer.h"
#include <CGAL/glu.h>



using namespace std;

#include "Viewer.moc" // .moc will be the output from moc preprocessor

void Viewer::init()
{
  /* Initial timer for playing incremental construction */
  m_pTimer = new QTimer(this);
  connect(m_pTimer, SIGNAL(timeout()), this, SLOT(incremental_insert()));

  /* Scene inits */
  setBackgroundColor(::Qt::white);
  // scene are defined by a sphere of 2.0, camera at the center, i.e. (0, 0, 0)
  setSceneCenter( qglviewer::Vec(-0.,-0.,-0.) );
  setSceneRadius( 2. );
  // show text message
  setTextIsEnabled(true);
  setForegroundColor(::Qt::red);
  setFont(QFont("Arial Black", 16, QFont::Bold));

  /* OpenGL inits */
  // Increase the material shininess, so that the difference between
  // the two versions of the spiral is more visible.
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);
  GLfloat specular_color[4] = { 0.8f, 0.8f, 0.8f, 1.0 };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  specular_color);
  // Set Smooth Shading
  ::glShadeModel(GL_SMOOTH);

  // depth buffer setup 
  ::glClearDepth(1.0f);
  ::glEnable(GL_DEPTH_TEST);
  ::glDepthFunc(GL_LEQUAL);
  ::glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  // enable semi-transparent culling planes
  ::glEnable(GL_BLEND);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // anti-aliasing, i.e. reduce jaggedness (if the OpenGL driver permits that)
  ::glEnable(GL_POINT_SMOOTH);
  ::glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  ::glEnable(GL_LINE_SMOOTH);
  ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  /* Add mouse and key description */
  setKeyDescription( Qt::CTRL + Qt::Key_G, tr("Generate points") );
  setKeyDescription( Qt::CTRL + Qt::Key_O, tr("Load points") );
  setKeyDescription( Qt::CTRL + Qt::Key_S, tr("Save points") );
  setKeyDescription( Qt::CTRL + Qt::Key_Comma, tr("Preference") );
  setKeyDescription( Qt::CTRL + Qt::Key_H, tr("Hide Kernel Demo") );
  setKeyDescription( Qt::CTRL + Qt::Key_Q, tr("Quit Kernel Demo") );
  setKeyDescription( Qt::Key_Return,
         tr("Insert new point to triangulation in <u>Input-Point</u> mode") );
  setKeyDescription( Qt::Key_Escape,
         tr("Cancel insertion in <u>Input-Point</u> mode;<br>")
         + tr("Cancel current selection in <u>Select</u> mode") );
  setKeyDescription( Qt::Key_Delete, tr("Delete selected vertices in <u>Select</u> mode") );

  setMouseBindingDescription( Qt::LeftButton,
         tr("Hold to move new point in <u>Input-Point</u> mode;<br>")
         + tr("Hold to move a vertex in <u>Move</u> mode") );
  setMouseBindingDescription( Qt::SHIFT + Qt::LeftButton,
         tr("Click to insert a vertex in <u>Input-Vertex</u> mode;<br>")
         + tr("Click to insert a point in <u>Input-Point</u> mode;<br>")
         + tr("Click or Drag to select multiple points in <u>Select</u> mode;<br>")
         + tr("Click to place a query point in <u>Find-Nearest-Neighbor</u> mode;<br>")
         + tr("Click to place a query point in <u>Show-Empty-Sphere</u> mode") );
  setMouseBindingDescription( Qt::CTRL + Qt::LeftButton,
         tr("Drag to add vertices to current selection in <u>Select</u> mode") );
}

QString Viewer::helpString() const
{
  QString text("<h1>3D Triangulation Demo</h1>");

  text += "This example illustrates a generic interactive demo for 3D Triangulation in CGAL. ";
  text += "This demo could be used as a simple skeleton ";
  text += "for potential demos of other 3D packages or for teaching CGAL.<br><br>";

  text += "The key feature is to edit vertices/points with mouse.";
  text += "There are several modes:<br><br>";

  text += " - <u>Normal Mode</u>: ";
  text += "Rotate, zoom, or translate camera using mouse.<br>";
  text += " - <u>Insert Vertex</u>: ";
  text += "Insert a vertex on the surface of the trackball ";
  text += "and the triangulation will be updated correspondingly.<br>";
  text += " - <u>Insert Point</u>: ";
  text += "Insert a point on the surface of the trackball. ";
  text += "Its conflict region will be highlighted. ";
  text += "When the new point is moving, ";
  text += "its conflict region will be updated correspondingly.<br>";
  text += " - <u>Select</u>: ";
  text += "Click or drag mouse left button to select multiple points.<br>";
  text += " - <u>Move</u>: Hold mouse left button to move a vertex ";
  text += "and the triangulation will be updated correspondingly.<br>";
  text += " - <u>Find Nearest Neighbor</u>: ";
  text += "Place a query point and its nearest neighbor will be highlighted.<br>";
  text += " - <u>Show Empty Sphere</u>: ";
  text += "Place a query point, locate the point in a cell ";
  text += "and then show the empty sphere of that cell. ";
  text += "An empty sphere of a cell is a sphere ";
  text += "with all four vertices of the cell lying on it ";
  text += "and no other vertices inside it.<br><br>";
  text += "<b>Shift+Wheel</b> to resize the trackball when it exists. ";
  text += "See <b>Mouse</b> page for more details.<br><br>";

  text += "Other basic features include:<br>";
  text += " - Randomly generate points,<br>";
  text += " - Read/Write files,<br>";
  text += " - Show vertices, Voronoi edges, Delaunay edges, and/or facets,<br>";
  text += " - Incremental Construct: ";
  text += "Re-construct the current triangulation incrementally. ";
  text += "If no triangulation exists yet, randomly generate 100 points ";
  text += "and construct a Delaunay triangulation of those points.<br>";

  return text;
}

/*************************************************************/
/*  Draw functions */

void Viewer::draw()
{
  if( m_pScene == NULL ) return;

  QFont fontPrompt("Arial", 14);

  if( m_showAxis ) {
    qglColor(::Qt::black);
    drawAxis( sceneRadius() );
  }

  /* Draw vertices */
  if ( m_showVertex && m_pScene->m_dt.number_of_vertices()>0 ) {
    for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
        vit < m_pScene->m_vhArray.end(); ++vit) {
      if( m_curMode == SELECT && (*vit)->isSeled() )  continue;
      if( (*vit) == m_nearestNb ) continue;
      drawVertex( (*vit)->point(), m_colorVertex, m_fSizeVertex );
    }//end-for-points
  }//end-if-points

  /* Draw all points during incremental mode */
  if( !m_incrementalPts.isEmpty() ) {
    /* draw the rest to-be-inserted vertices */
    for(QList<Point_3>::iterator pit=m_incrementalPts.begin();
        pit < m_incrementalPts.end(); ++pit) {
      drawVertex( (*pit), ::Qt::gray, m_fSizeVertex );
    }

  	switch( m_curStep ) {
  	case NEWPT:
      /* Show prompt messages */
      qglColor( ::Qt::black );
      drawText( 10, 20, tr("Highlight the next-to-insert point"), fontPrompt );
      /* Highlight the next-to-insert point */
      drawVertex( m_curIncPt, ::Qt::red, m_fSizeVertex );
      break;
    case CELL:  // show the tetrahedron that contains the point
      /* Show prompt messages */
      qglColor( ::Qt::black );
      drawText( 10, 20, tr("Show the tetrahedron containing the point"), fontPrompt );
      drawText( 10, 40, tr("(Only finite facets are drawn)"), fontPrompt );
      /* Highlight the next-to-insert vertex */
      drawVertex( m_curIncPt, ::Qt::red, m_fSizeVertex );
      /* Draw the cell containing that point */
      for(int i=0; i<4; ++i) {
        if( m_pScene->m_dt.is_infinite(m_cellContain, i) )  continue;
        drawFacet( m_pScene->m_dt.triangle( m_cellContain, i ), m_colorFacet );
      }//end-for-facets
      break;
    case CONFLICT:  // show the conflict region
      /* Show prompt messages */
      qglColor( ::Qt::black );
      drawText( 10, 20, tr("Show the conflict region"), fontPrompt );
      /* Highlight the next-to-insert vertex */
      drawVertex( m_curIncPt, ::Qt::red, m_fSizeVertex );
      /* Draw conflict region */
      for(QList<Facet>::iterator fit = m_boundaryFacets.begin();
          fit < m_boundaryFacets.end(); ++fit) {
        if( m_pScene->m_dt.is_infinite(*fit) )  continue;
        drawFacet( m_pScene->m_dt.triangle(*fit), QColor(215, 80, 0, 96) ); //semi-transparent purple
      }//end-for-facets
      break;
    default:
      break;
  	}//end-of=switch
  }//end-if-incpts

  /* Draw Delaunay edges */
  if( m_showDEdge ) {
    for(edges_iterator eit = m_pScene->m_dt.finite_edges_begin();
        eit != m_pScene->m_dt.finite_edges_end(); ++eit) {
      Segment_3 seg = m_pScene->m_dt.segment(*eit);
      drawEdge( seg.vertex(0), seg.vertex(1), m_colorDEdge, m_fSizeDEdge );
    }//end-for-edges
  }//end-if-dt

  /* Draw Voronoi edges */
  if( m_showVEdge ) {
    for(facets_iterator fit = m_pScene->m_dt.finite_facets_begin();
        fit != m_pScene->m_dt.finite_facets_end(); ++fit) {
      Object_3 o = m_pScene->m_dt.dual(*fit);
      if (const Segment_3 *s = CGAL::object_cast<Segment_3>(&o)) {
        drawEdge( s->vertex(0), s->vertex(1), m_colorVEdge, m_fSizeVEdge );
      } else if (const Ray_3 *r = CGAL::object_cast<Ray_3>(&o)) {
        drawEdge( r->point(0),  // the source of the ray
                  r->point(1),  // another point on the ray, different from the source
                  m_colorVEdge, m_fSizeVEdge );
      }
    }//end-for-edges
  }//end-if-vd

  /* Draw facets */
  if( m_showFacet ) {
    for(facets_iterator fit = m_pScene->m_dt.finite_facets_begin();
        fit != m_pScene->m_dt.finite_facets_end(); ++fit) {
      drawFacet( m_pScene->m_dt.triangle(*fit), m_colorFacet );
    }//end-for-facets
  }//end-if-facets

  /* Insert vertex mode */
  if( m_curMode == INSERT_V ) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Shift+Left: Insert a vertex"), fontPrompt );
    drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
    /* Draw the trackball */
    drawSphere( m_fRadius, m_colorTrackball );
  }//end-if-insv

  /* Insert point mode */
  else if( m_curMode == INSERT_PT ) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Shift+Left: Insert a point"), fontPrompt );
    drawText( width()-200, 40, tr("Hold Left: Move the point"), fontPrompt );
    drawText( width()-200, 60, tr("Return: Insert to DT"), fontPrompt );
    drawText( width()-200, 80, tr("Escape: Cancel insertion"), fontPrompt );
    drawText( width()-200, 100, tr("Shift+Wheel: Resize trackball"), fontPrompt );

    /* Draw the trackball */
    drawSphere( m_fRadius, m_colorTrackball );

    if( m_hasNewPt ) {
      /* Draw the newly inserted point */
      drawVertex( m_newPt, ::Qt::red, m_fSizeVertex );
      /* Draw conflict region */
      for(QList<Facet>::iterator fit = m_boundaryFacets.begin();
          fit < m_boundaryFacets.end(); ++fit) {
        if( m_pScene->m_dt.is_infinite(*fit) )  continue;
        drawFacet( m_pScene->m_dt.triangle(*fit), QColor(215, 80, 0, 96) ); //semi-transparent purple
      }//end-for-facets
    }//end-if-shown
  }//end-if-inspt

  /* Select mode */
  else if( m_curMode == SELECT) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Shift+Left: Select"), fontPrompt );
    drawText( width()-200, 40, tr("Ctrl+Left: Add selection"),
              QFont("Arial", 14) );
    drawText( width()-200, 60, tr("Escape: Cancel selection"), fontPrompt );
    drawText( width()-200, 80, tr("DEL: Delete selected"), fontPrompt );
    /* Highlight the selected vertices */
    for(QList<int>::iterator vit=m_vidSeled.begin(); vit<m_vidSeled.end(); ++vit) {
      drawVertex( m_pScene->m_vhArray.at(*vit)->point(), ::Qt::red, m_fSizeVertex );
    }//end-for-seledpts
    /* Draw the multiple selection window */
    if( m_isPress ) {
      ::glDisable( GL_LIGHTING );
      startScreenCoordinatesSystem();
      qglColor( QColor(80, 180, 180, 64) );
      ::glBegin(GL_QUADS);
      ::glVertex2i(m_rectSel.left(), m_rectSel.top());
      ::glVertex2i(m_rectSel.right(), m_rectSel.top());
      ::glVertex2i(m_rectSel.right(), m_rectSel.bottom());
      ::glVertex2i(m_rectSel.left(), m_rectSel.bottom());
      ::glEnd();
      stopScreenCoordinatesSystem();
      ::glEnable( GL_LIGHTING );
    }//end-if-press
  }//end-if-sel

  /* Move mode */
  else if( m_curMode == MOVE ) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Left Click: Select"), fontPrompt );

    if( m_isMoving ) {
      drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
      /* Draw the trackball */
      drawSphere( m_fRadius, m_colorTrackball );
      /* Highlight the moving point */
      drawVertex( m_pScene->m_vhArray.at( m_vidMoving )->point(), ::Qt::red, m_fSizeVertex );
    }//end-if-v
  }//end-if-move

  /* FindNb mode */
  else if( m_curMode == FINDNB ) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
    drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
    /* Draw the trackball */
    drawSphere( m_fRadius, m_colorTrackball );
    /* Draw the nearest neighbor */
    if( m_nearestNb != NULL ) {
      drawVertex( m_queryPt, ::Qt::red, m_fSizeVertex );
      drawVertex( m_nearestNb->point(), ::Qt::red, m_fSizeVertex );
    }
  }//end-if-findnb

  /* EmptySphere mode */
  else if( m_curMode == EMPTYSPH ) {
    /* Show prompt messages */
    qglColor( ::Qt::black );
    drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
    drawText( width()-200, 40, tr("Press S: Show/Hide trackball"), fontPrompt );
    drawText( width()-200, 60, tr("Shift+Wheel: Resize trackball"), fontPrompt );
    /* Draw the trackball */
    if( m_showTrackball )
      drawSphere( m_fRadius, m_colorTrackball );

    if( m_hasEmptyS ) {
      /* Draw the query point */
      drawVertex( m_queryPt, ::Qt::red, m_fSizeVertex );
      /* Draw the cell containing that point */
      for(int i=0; i<4; ++i) {
        if( m_pScene->m_dt.is_infinite(m_cellContain, i) )  continue;
        drawFacet( m_pScene->m_dt.triangle( m_cellContain, i ), m_colorFacet );
      }//end-for-facets
      /* Draw the sphere */
      drawSphere( m_fREmptyS, m_colorEmptySphere, m_centerPt );
    }
  }//end-if-emptyS
}

void Viewer::drawVertex(const Point_3& p, const QColor& clr, float r)
{
  /* Draw regular points */
  if( m_isFlat ) {
    // disable lighting
    ::glDisable( GL_LIGHTING );

    ::glPointSize(8.0);
    qglColor( clr );

    ::glBegin(GL_POINTS);
    ::glVertex3f( p.x(), p.y(), p.z() );
    ::glEnd();

    // resume lighting
    ::glEnable( GL_LIGHTING );

    return;
  }

  /* Draw vertices as 3D balls */
  GLboolean lighting, colorMaterial;
  ::glGetBooleanv( GL_LIGHTING, &lighting );
  ::glGetBooleanv( GL_COLOR_MATERIAL, &colorMaterial );
  ::glEnable( GL_LIGHTING );
  ::glDisable(GL_COLOR_MATERIAL);

  float color[4];
  color[0] = clr.redF();
  color[1] = clr.greenF();
  color[2] = clr.blueF();
  color[3] = clr.alphaF();

  // move to the point
  ::glPushMatrix();
  ::glTranslatef( p.x(), p.y(), p.z() );

  // draw
  GLUquadricObj* quadratic = ::gluNewQuadric();	// Create A Pointer To The Quadric Object
  ::gluQuadricNormals( quadratic, GLU_SMOOTH );	// Create Smooth Normals
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color );
  ::gluSphere( quadratic, r, 16, 16 );

  // move back to origin
  ::glPopMatrix();

  if ( colorMaterial )
    ::glEnable( GL_COLOR_MATERIAL );
  if ( !lighting )
    ::glDisable( GL_LIGHTING );
}

void Viewer::drawEdge(const Point_3& from, const Point_3& to, const QColor& clr, float r)
{
  /* Draw regular lines */
  if( m_isFlat ) {
    // disable lighting
    ::glDisable( GL_LIGHTING );

    ::glLineWidth(1.0);
  	qglColor( clr );

    ::glBegin(GL_LINES);
    ::glVertex3f( from.x(), from.y(), from.z() );
    ::glVertex3f( to.x(), to.y(), to.z() );
    ::glEnd();

    // resume lighting
    ::glEnable( GL_LIGHTING );

    return;
  }

  /* Draw edges as 3D cylinders */
  GLboolean lighting, colorMaterial;
  ::glGetBooleanv( GL_LIGHTING, &lighting );
  ::glGetBooleanv( GL_COLOR_MATERIAL, &colorMaterial );
  ::glEnable( GL_LIGHTING );
  ::glDisable(GL_COLOR_MATERIAL);

  float color[4];
  color[0] = clr.redF();
  color[1] = clr.greenF();
  color[2] = clr.blueF();
  color[3] = clr.alphaF();

  Vector_3 v = to - from;

  // compute the length of the edge
  // method 1:
//  float length = sqrt( CGAL::squared_distance( from, to ) );
  // method 2:
  float length = sqrt( v.squared_length() );

  // normalize
  v = v / length;
  // compute the angle: cos theta = v.z/1.0
  GLfloat angle = acos( v.z() ) / 3.1415927 * 180;

  ::glPushMatrix();

  // move to "from" point
  ::glTranslatef( from.x(), from.y(), from.z() );
  // rotate from z-axis to from-->to
  //  axis: cross product of z-axis and from-->to
  ::glRotatef( angle, -v.y(), v.x(), 0.0f );
  // draw
  GLUquadricObj* quadratic = ::gluNewQuadric();	// Create A Pointer To The Quadric Object
  ::gluQuadricNormals( quadratic, GLU_SMOOTH );	// Create Smooth Normals
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color );
  // gluCylinder draws a cylinder oriented along the z-axis
  ::gluCylinder( quadratic, r, r, length, 16, 4 );

  // move back to origin
  ::glPopMatrix();

  if ( colorMaterial )
    ::glEnable( GL_COLOR_MATERIAL );
  if ( !lighting )
    ::glDisable( GL_LIGHTING );
}

void Viewer::drawFacet(const Triangle_3& t, const QColor& /*clr*/)
{
  // disable lighting
  ::glDisable( GL_LIGHTING );

  // disable depth buffer writing
  ::glDepthMask( GL_FALSE );

  qglColor( m_colorFacet );

  ::glBegin(GL_TRIANGLES);
  Point_3 p0 = t.vertex(0);
  Point_3 p1 = t.vertex(1);
  Point_3 p2 = t.vertex(2);
  ::glVertex3f( p0.x(), p0.y(), p0.z() );
  ::glVertex3f( p1.x(), p1.y(), p1.z() );
  ::glVertex3f( p2.x(), p2.y(), p2.z() );
  ::glEnd();

  // resume depth buffer writing
  ::glDepthMask( GL_TRUE );

  // resume lighting
  ::glEnable( GL_LIGHTING );
}

void Viewer::drawSphere(float r, const QColor& clr, const Point_3& center)
{
  GLboolean lighting, colorMaterial;
  ::glGetBooleanv( GL_LIGHTING, &lighting );
  ::glGetBooleanv( GL_COLOR_MATERIAL, &colorMaterial );
  ::glEnable( GL_LIGHTING );
  ::glDisable(GL_COLOR_MATERIAL);

  float color[4];
  color[0] = clr.redF();
  color[1] = clr.greenF();
  color[2] = clr.blueF();
  color[3] = clr.alphaF();

  ::glPushMatrix();

  // move to the point
  if( center != CGAL::ORIGIN )  ::glTranslatef( center.x(), center.y(), center.z() );

  // disable depth buffer writing
  ::glDepthMask( GL_FALSE );
  // draw
  GLUquadricObj* quadratic = ::gluNewQuadric();	// Create A Pointer To The Quadric Object
  ::gluQuadricNormals( quadratic, GLU_SMOOTH );	// Create Smooth Normals
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color );
  ::gluSphere( quadratic, r, 32, 32 );
  // resume depth buffer writing
  ::glDepthMask( GL_TRUE );

  // move back to origin
  ::glPopMatrix();

  if ( colorMaterial )
    ::glEnable( GL_COLOR_MATERIAL );
  if ( !lighting )
    ::glDisable( GL_LIGHTING );
}

/*************************************************************/
/*  Select functions */

void Viewer::drawWithNames()
{
  for(int i=0; i<m_pScene->m_vhArray.size(); ++i) {
    // push a name for each point onto the name stack
    // note: it can NOT be used between glBegin and glEnd
    ::glPushName( i );

    // draw the point
    ::glBegin(GL_POINTS);
    Point_3& p = m_pScene->m_vhArray.at(i)->point();
    ::glVertex3f(p.x(), p.y(), p.z());
    ::glEnd();

    // pop one name off the top of the name stack
    ::glPopName();
  }//end-for-points

  // push a name for the newly inserted point
  if( m_curMode == INSERT_PT && m_hasNewPt ) {
    ::glPushName( ::GLuint(-1) );
    ::glBegin(GL_POINTS);
    ::glVertex3f(m_newPt.x(), m_newPt.y(), m_newPt.z());
    ::glEnd();
    ::glPopName();
  }//end-if-newPt
}

void Viewer::endSelection(const QPoint& /*point*/)
{
  // flush GL buffers
  ::glFlush();

  // reset GL_RENDER mode (was GL_SELECT) and get the number of selected points
  size_t nSel = ::glRenderMode(GL_RENDER);

  /* No selection */
  if( nSel <= 0 ) {
    if( m_curMode == SELECT )
      m_isPress = false;
  }//end-if-notselected

  // each hit record has 4 data: # of names in name stack, min and max depth of old hits,
  //  name stack contents [see glSelectBuffer man page for more details]
  //  i.e. (selectBuffer())[4*i+3] is the id pushed on the stack

  /* Check whether the new point is clicked on */
  else if( m_curMode == INSERT_PT ) {
  	if( m_hasNewPt && (selectBuffer())[3] == ::GLuint(-1) )
      m_isMoving = true;
  }//end-if-inspt

  /* Check whether vertex is clicked on */
  else if( m_curMode == MOVE ) {
    m_isMoving = true;
    m_vidMoving = (selectBuffer())[3];
    // compute the corresponding size of trackball, i.e. selectedV is on the ball
    Point_3 p = m_pScene->m_vhArray.at( m_vidMoving )->point();
    m_fRadius = sqrt( p.x()*p.x() + p.y()*p.y() + p.z()*p.z() );
  }//end-if-move

  /* Store current selections */
  else { // m_curMode == SELECT
    if( m_selMode == NORMAL ) {
      // remove the old selections
      for(QList<int>::iterator vit=m_vidSeled.begin();
          vit < m_vidSeled.end(); ++vit) {
        m_pScene->m_vhArray.at(*vit)->setSeled( false );
      }
      m_vidSeled.clear();

      // record the new selections
      for(std::size_t i=0; i<nSel; ++i) {
        m_vidSeled.push_back( (selectBuffer())[4*i+3] );
        m_pScene->m_vhArray.at( m_vidSeled.back() )->setSeled();
      }
    } else {
      for(std::size_t i=0; i<nSel; ++i) {
        if( !m_vidSeled.contains( (selectBuffer())[4*i+3] ) ) {
          m_vidSeled.push_back( (selectBuffer())[4*i+3] );
          m_pScene->m_vhArray.at( (selectBuffer())[4*i+3] )->setSeled();
        }//end-if-contain
      }//end-for
    }//end-if-add
  }//end-if-sel
}

/*************************************************************/
/*  Mouse and Keyboard functions */

void Viewer::mousePressEvent(QMouseEvent *event)
{
  // button() holds the button that caused the event
  //  note: for mouse move event, button() always return Qt::NoButton
  // modifiers() holds the keyboard modifier flags at the time of the event
  // buttons() holds the button state when the event was generated,
  //  i.e. all buttons that are pressed down
  // pos() holds the mouse cursor's position relative to the receiving widget

  // Get event modifiers key
#if QT_VERSION < 0x040000
  // Bug in Qt : use 0x0f00 instead of Qt::KeyButtonMask with Qt versions < 3.1
  const Qt::ButtonState modifiers = (Qt::ButtonState)(event->state() & Qt::KeyButtonMask);
#else
  const Qt::KeyboardModifiers modifiers = event->modifiers();
#endif

  if( m_curMode == INSERT_V
      && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
  	m_isPress = true;
  }//end-if-insv

  else if(m_curMode == INSERT_PT && event->button() == Qt::LeftButton ) {
    /* shift+left to insert */
    if( modifiers == Qt::SHIFT ) {
      if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
        m_isPress = true;
      else
        displayMessage( tr("There exists no triangulation yet.") );
      m_hasNewPt = false;
    } else { /* left button to move */
      m_isMoving = false;
      // define selection window (default was 3)
      setSelectRegionWidth( 10 );
      setSelectRegionHeight( 10 );
      // perform the selection
      select( event->pos() );
      if( m_isMoving )
        // redraw window
        updateGL();
      else
        // if no point is selected, then regular action (rotation) will be performed
        QGLViewer::mousePressEvent(event);
    }//end-if-shift
  }//end-if-inspt

  else if( m_curMode == SELECT && event->button() == Qt::LeftButton ) {
    // set the selection mode
    switch( modifiers ) {
    case Qt::SHIFT : // select
      m_isPress = true;
      m_selMode = NORMAL;
      // initialize multiple selection window
      m_rectSel = QRect( event->pos(), event->pos() );
      // redraw window
      updateGL();
      break;
    case Qt::CTRL : // add selection
      m_isPress = true;
      m_selMode = ADD;
      // initialize multiple selection window
      m_rectSel = QRect( event->pos(), event->pos() );
      // redraw window
      updateGL();
      break;
    default: // rotate
      QGLViewer::mousePressEvent(event);
      break;
    }
  }//end-if-select

  else if(m_curMode == MOVE && event->button() == Qt::LeftButton ) {
  	m_isMoving = false;
    // define selection window (default was 3)
    setSelectRegionWidth( 10 );
    setSelectRegionHeight( 10 );
    // perform the selection
    select( event->pos() );
    if( m_isMoving ) // redraw window
      updateGL();
    else // if no point is selected, then regular action (rotation) will be performed
      QGLViewer::mousePressEvent(event);
  }//end-if-move

  else if( m_curMode == FINDNB
      && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
    if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
      m_isPress = true;
    else
      displayMessage( tr("There exists no triangulation yet.") );
  }//end-if-findnb

  else if( m_curMode == EMPTYSPH
      && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
    if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
      m_isPress = true;
    else
      displayMessage( tr("There exists no triangulation yet.") );
    m_hasEmptyS = false;
  }//end-if-emptyS

  else
    QGLViewer::mousePressEvent(event);
}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
  if( m_curMode == INSERT_PT && m_isMoving ) {
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      m_newPt = Point_3(pt.x, pt.y, pt.z);
      // compute the conflict hole induced by point p
      computeConflict( m_newPt );
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-inspt

  else if( m_curMode == SELECT && m_isPress ) {
    // update multiple selection window
    m_rectSel.setBottomRight( event->pos() );
    // redraw
    updateGL();
  }//end-if-sel

  else if( m_curMode == MOVE && m_isMoving ) {
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      // note: QList::operator[] return a modifiable reference;
      //   while QList::at return a const reference
#if CGAL_VERSION_NR < 1030700000
      // move_point moves the point stored in v to p while preserving the Delaunay property
      // it calls remove(v) followed by insert(p) and return the new handle
      // it supposely faster when the point has not moved much
      m_pScene->m_vhArray[m_vidMoving] = m_pScene->m_dt.move_if_no_collision(
                                 m_pScene->m_vhArray.at( m_vidMoving ),
                                 Point_3( pt.x, pt.y, pt.z ) );
#else
      // move_if_no_collision moves the point stored in v to pt
      //  if there is not already another vertex placed on pt,
      //  the triangulation is modified s.t. the new position of v is pt;
      //  otherwise, the vertex at point pt is returned.
      Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                                 m_pScene->m_vhArray.at( m_vidMoving ),
                                 Point_3( pt.x, pt.y, pt.z ) );
      int id1 = m_pScene->m_vhArray.indexOf( vh );
      int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
      // remove the duplicate in vhArray
      if( id1 != m_vidMoving )
        m_pScene->m_vhArray.removeAt( id1 );
      else if( id2 != -1 )
        m_pScene->m_vhArray.removeAt( id2 );
      m_pScene->m_vhArray[m_vidMoving] = vh;
#endif
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-move

  else
    QGLViewer::mouseMoveEvent(event);
}

void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
  /* INS_V mode - Shift+Left: compute and insert a vertex */
  if( m_curMode == INSERT_V && m_isPress ) {
    m_isPress = false;
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      m_pScene->m_vhArray.push_back( m_pScene->m_dt.insert( Point_3( pt.x, pt.y, pt.z ) ) );
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-ins

  /* INS_PT mode - Shift+Left: compute and insert a point */
  else if( m_curMode == INSERT_PT && m_isPress ) {
    m_isPress = false;
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      m_hasNewPt = true;
      m_newPt = Point_3(pt.x, pt.y, pt.z);
      // compute the conflict hole induced by point p
      computeConflict( m_newPt );
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-inspt

  /* INS_PT mode - Left: compute and insert a point */
  else if( m_curMode == INSERT_PT && m_isMoving ) {
    m_isMoving = false;
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      m_newPt = Point_3(pt.x, pt.y, pt.z);
      // compute the conflict hole induced by point p
      computeConflict( m_newPt );
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-inspt

  /* SEL mode - Left: terminate multiple point selection */
  else if( m_curMode == SELECT && m_isPress ) {
    // might swap left/right and top/bottom to make rectanle valid
#if QT_VERSION < 0x040000
    m_rectSel = m_rectSel.normalize();
#else
    m_rectSel = m_rectSel.normalized();
#endif

    if( m_rectSel.width() == 1 && m_rectSel.height() == 1 ) { /* select a point */
      // set a default selection window
      setSelectRegionWidth( 10 );
      setSelectRegionHeight( 10 );
      // compute rectangle center and perform selection
      select( m_rectSel.center() );
      if( m_isPress ) {
        m_isPress = false;
      } else {
        displayMessage( tr("No point is selected.") );
      }
    } else {  /* select multiple points, ie. selection window > 1 */
      // define selection window
      if( m_rectSel.width() < 10 )
        setSelectRegionWidth( 10 );
      else
        setSelectRegionWidth( m_rectSel.width() );
      if( m_rectSel.height() < 10 )
        setSelectRegionHeight( 10 );
      else
        setSelectRegionHeight( m_rectSel.height() );
      // compute rectangle center and perform selection
      select( m_rectSel.center() );
      if( m_isPress ) {
        m_isPress = false;
        displayMessage( QString::number(m_vidSeled.size()) + tr(" points are selected") );
      } else { // empty window will cancel the current selection
        for(QList<int>::iterator iit = m_vidSeled.begin(); iit < m_vidSeled.end(); ++iit)
          m_pScene->m_vhArray.at(*iit)->setSeled( false );
  	    m_vidSeled.clear();
      }
    }//end-if-selwindow

    // update display to show
    updateGL();
  }//end-if-select

  /* MOVE mode - Left: terminate point moving */
  else if( m_curMode == MOVE && m_isMoving ) {
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      // note: QList::operator[] return a modifiable reference;
      //   while QList::at return a const reference
#if CGAL_VERSION_NR < 1030700000
      // move_point moves the point stored in v to p while preserving the Delaunay property
      // it calls remove(v) followed by insert(p) and return the new handle
      // it supposely faster when the point has not moved much
      m_pScene->m_vhArray[m_vidMoving] = m_pScene->m_dt.move_if_no_collision(
                                 m_pScene->m_vhArray.at( m_vidMoving ),
                                 Point_3( pt.x, pt.y, pt.z ) );
#else
      // move_if_no_collision moves the point stored in v to pt
      //  if there is not already another vertex placed on pt,
      //  the triangulation is modified s.t. the new position of v is pt;
      //  otherwise, the vertex at point pt is returned.
      Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                                 m_pScene->m_vhArray.at( m_vidMoving ),
                                 Point_3( pt.x, pt.y, pt.z ) );
      int id1 = m_pScene->m_vhArray.indexOf( vh );
      int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
      // remove the duplicate in vhArray
      if( id1 != m_vidMoving )
        m_pScene->m_vhArray.removeAt( id1 );
      else if( id2 != -1 )
        m_pScene->m_vhArray.removeAt( id2 );
      m_pScene->m_vhArray[m_vidMoving] = vh;
#endif
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-move

  /* FindNb mode - Shift+Left: find the nearest neighbor of the point */
  else if( m_curMode == FINDNB && m_isPress ) {
    m_isPress = false;
    Vec pt;
    if( computeIntersect( event->pos(), pt ) ) {
      m_queryPt = Point_3( pt.x, pt.y, pt.z );
      m_nearestNb = m_pScene->m_dt.nearest_vertex( m_queryPt );
    }//end-if-compute

    // redraw
    updateGL();
  }//end-if-findnb

  /* EmptySphere mode - Shift+Left: show the empty sphere of the cell */
  else if( m_curMode == EMPTYSPH && m_isPress ) {
    m_isPress = false;
    Vec pt;
    m_hasEmptyS = computeIntersect( event->pos(), pt );
    if( m_hasEmptyS ) {
      m_queryPt = Point_3( pt.x, pt.y, pt.z );
      // find the cell that contains point p in its interior
      m_cellContain = m_pScene->m_dt.locate( m_queryPt );
      // show error if point is outside the convex hull
      if( m_pScene->m_dt.is_infinite( m_cellContain ) ) {
        m_hasEmptyS = false;
        displayMessage( tr("Query point is outside the convex hull!") );
      } else { /* compute the empty sphere */
        // find the circumcenter of the four vertices of c
        m_centerPt = m_pScene->m_dt.dual( m_cellContain );
        // compute the radius of the empty sphere
        m_fREmptyS = sqrt( CGAL::squared_distance( m_centerPt,
                             m_cellContain->vertex(0)->point() ) );
      }
    }//end-if-compute
    // redraw
    updateGL();
  }//end-if-emptysphere

  else
    QGLViewer::mouseReleaseEvent(event);
}

void Viewer::wheelEvent(QWheelEvent *event)
{
  // Get event modifiers key
#if QT_VERSION < 0x040000
  // Bug in Qt : use 0x0f00 instead of Qt::KeyButtonMask with Qt versions < 3.1
  const Qt::ButtonState modifiers = (Qt::ButtonState)(event->state() & Qt::KeyButtonMask);
#else
  const Qt::KeyboardModifiers modifiers = event->modifiers();
#endif

  if( (m_curMode == INSERT_V || m_curMode == FINDNB || m_curMode == EMPTYSPH )
     && modifiers == Qt::SHIFT ) {
    // delta() returns the distance that the wheel is rotated, in eighths of a degree.
    //  note: most mouse types work in steps of 15 degrees
    //  positive value: rotate forwards away from the user;
    //  negative value: rotate backwards toward the user.
    m_fRadius += (event->delta()*1. / m_iStep ); // inc-/decrease by 0.1 per step
    if( m_fRadius < 0.1 )
      m_fRadius = 0.1f;

    // redraw
    updateGL();
  }//end-if-insv

  else if( m_curMode == INSERT_PT && modifiers == Qt::SHIFT ) {
    // delta() returns the distance that the wheel is rotated, in eighths of a degree.
    //  note: most mouse types work in steps of 15 degrees
    //  positive value: rotate forwards away from the user;
    //  negative value: rotate backwards toward the user.
  	float origR = m_fRadius;
    m_fRadius += (event->delta()*1. / m_iStep ); // inc-/decrease by 0.1 per step
    if( m_fRadius < 0.1 )
      m_fRadius = 0.1f;
    // update the new point and its conflict region
    if( m_hasNewPt ) {
      origR = m_fRadius / origR;
      m_newPt = Point_3( m_newPt.x()*origR, m_newPt.y()*origR, m_newPt.z()*origR );
      // compute the conflict hole induced by point p
      computeConflict( m_newPt );
    }//end-if-conflict

    // redraw
    updateGL();
  }//end-if-inspt

  // resize the trackball when moving a point
  else if( m_curMode == MOVE && modifiers == Qt::SHIFT && m_isMoving ) {
  	float origR = m_fRadius;
    m_fRadius += (event->delta()*1. / m_iStep ); // inc-/decrease by 0.1 per step
    if( m_fRadius < 0.1 )
      m_fRadius = 0.1f;
    origR = m_fRadius / origR;
    Point_3 pt = m_pScene->m_vhArray.at( m_vidMoving )->point();
    // note: QList::operator[] return a modifiable reference;
    //   while QList::at return a const reference
#if CGAL_VERSION_NR < 1030700000
    // move_point moves the point stored in v to p while preserving the Delaunay property
    // it calls remove(v) followed by insert(p) and return the new handle
    // it supposely faster when the point has not moved much
    m_pScene->m_vhArray[m_vidMoving] = m_pScene->m_dt.move_if_no_collision(
                               m_pScene->m_vhArray.at( m_vidMoving ),
                               Point_3( pt.x()*origR, pt.y()*origR, pt.z()*origR ) );
#else
    // move_if_no_collision moves the point stored in v to pt
    //  if there is not already another vertex placed on pt,
    //  the triangulation is modified s.t. the new position of v is pt;
    //  otherwise, the vertex at point pt is returned.
    Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                               m_pScene->m_vhArray.at( m_vidMoving ),
                               Point_3( pt.x()*origR, pt.y()*origR, pt.z()*origR ) );
    int id1 = m_pScene->m_vhArray.indexOf( vh );
    int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
    // remove the duplicate in vhArray
    if( id1 != m_vidMoving )
      m_pScene->m_vhArray.removeAt( id1 );
    else if( id2 != -1 )
      m_pScene->m_vhArray.removeAt( id2 );
    m_pScene->m_vhArray[m_vidMoving] = vh;
#endif

    // redraw
    updateGL();
  }//end-if-move

  else
    QGLViewer::wheelEvent(event);
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
  // Get event modifiers key
#if QT_VERSION < 0x040000
  // Bug in Qt : use 0x0f00 instead of Qt::KeyButtonMask with Qt versions < 3.1
  const Qt::ButtonState modifiers = (Qt::ButtonState)(event->state() & Qt::KeyButtonMask);
#else
  const Qt::KeyboardModifiers modifiers = event->modifiers();
#endif

  /* Insert the newly inserted point as a vertex */
  if( m_curMode == INSERT_PT && m_hasNewPt
     && ( event->key()==Qt::Key_Return || event->key()==Qt::Key_Enter )
     && modifiers==Qt::NoButton ) {
    Facet& f = m_boundaryFacets.first(); // a boundary facet, i.e. a pair (cell_handle, i)
    // insert_in_hole will create a new vertex by starring a hole
    //   i.e. delete all conflict cells, create a new vertex,
    //        and for each boundary facet, create a new cell with the new vertex
    // it takes in an iterator range of conflict cells which specifies a hole
    //   and (begin, i) is a boundary facet that begin is one of the conflict cell
    //   but begin->neighbor(i) is not
    // it returns the handle of the new vertex
    m_pScene->m_vhArray.push_back( m_pScene->m_dt.insert_in_hole( m_newPt, // the point
                                     m_conflictCells.begin(), // cell_begin
                                     m_conflictCells.end(), // cell_end
                                     f.first, // cell_handle begin
                                     f.second ) ); // integer i

    m_hasNewPt = false;
    // erase old conflict hole info
    m_boundaryFacets.clear();
    m_conflictCells.clear();

  	// redraw
    updateGL();
  }//end-if-insVertex

  /* Cancel the newly inserted point and its conflict region */
  else if( m_curMode == INSERT_PT && m_hasNewPt
     && event->key()==Qt::Key_Escape && modifiers==Qt::NoButton ) {
    m_hasNewPt = false;
    // erase old conflict hole info
    m_boundaryFacets.clear();
    m_conflictCells.clear();

  	// redraw
    updateGL();
  }//end-if-escapeIns

  /* Delete selected points */
  else if( m_curMode == SELECT
     && event->key()==Qt::Key_Delete && modifiers==Qt::NoButton ) {
    // sort selected id's in descending order
    qSort(m_vidSeled.begin(), m_vidSeled.end(), qGreater<int>());
    for(QList<int>::iterator vit=m_vidSeled.begin(); vit<m_vidSeled.end(); ++vit) {
      // remove the selected point from DT and vertex_handle_array
      // note: QList::takeAt will removes the item at index position i and returns it.
      m_pScene->m_dt.remove( m_pScene->m_vhArray.takeAt( *vit ) );
  	}
  	// clear the selection buffer
  	m_vidSeled.clear();

  	// redraw
    updateGL();
  }//end-if-del

  /* Cancel the selection */
  else if( m_curMode == SELECT
     && event->key()==Qt::Key_Escape && modifiers==Qt::NoButton ) {
  	// clear the selection buffer
    for(QList<int>::iterator iit=m_vidSeled.begin(); iit<m_vidSeled.end(); ++iit) {
      m_pScene->m_vhArray.at(*iit)->setSeled( false );
    }
  	m_vidSeled.clear();

    // redraw
    updateGL();
  }//end-if-escapeSel


  /* Show/hide the trackball when drawing the empty sphere */
  else if( m_curMode == EMPTYSPH
     && event->key()==Qt::Key_S && modifiers==Qt::NoButton ) {
    m_showTrackball = !m_showTrackball;
    // redraw
    updateGL();
  }//end-if-showBall

  else
    QGLViewer::keyPressEvent(event);
}

/*************************************************************/
/*  Computation functions */

bool Viewer::computeIntersect( const QPoint & pos, Vec & pt )
{
  Vec eye, dir;
  // Compute eye position and direction to the clicked point,
  //  used to draw a representation of the intersecting line
  camera()->convertClickToLine( pos, eye, dir );
  // Compute the intersection point with the sphere
  //  note that the center of the sphere is at the origin (0, 0, 0)
  //  thus, (1) pt = eye + t*dir and (2) dist( pt, origin ) = radius
  //  i.e. (x_eye + t*x_dir)^2 + (y_eye + t*y_dir)^2 + (z_eye + t*z_dir)^2 = r^2
  //  --> t^2( dir*dir ) + 2t( eye*dir ) + eye*eye - r^2 = 0
  //      where "dir*dir" is the dot product of vector dir
  //  we need to solve t and the smaller t (nearer to eye position) is what we want
  float a = dir*dir;
  float b = eye*dir;
  float c = eye*eye - m_fRadius*m_fRadius;
  float delta = b*b - a*c;
  if( delta < 0 ) {
    displayMessage( tr("Point is not on the sphere!") );
    return false;
  } else {
    float t = ( (-1.)*b - sqrt(delta) ) / a;
    pt = eye + t*dir;
    return true;
  }
}

void Viewer::computeConflict( Point_3 pt )
{
  // find the cell that contains point p in its interior
  m_cellContain = m_pScene->m_dt.locate( pt );
  // erase old conflict hole info
  m_boundaryFacets.clear();
  m_conflictCells.clear();
  // show msg if point is outside the convex hull
  if( m_pScene->m_dt.is_infinite( m_cellContain ) )
    displayMessage( tr("Note: point is outside the convex hull.") );
  // compute the conflict hole induced by point p
  m_pScene->m_dt.find_conflicts( pt, // the point
                  m_cellContain, // starting cell that must be in conflict
                  std::back_inserter(m_boundaryFacets), // the facets on the boundary
                  std::back_inserter(m_conflictCells) ); // the cells in conflict
}

/*************************************************************/
/*  Animation functions */

void Viewer::toggleIncremental(bool on) {
  if( on ) {  // play
    if( m_incrementalPts.isEmpty() ) {
      /* start play */
      if( m_pScene->m_dt.number_of_vertices() == 0 ) {
        CGAL::Random_points_in_cube_3<Point_3> pts_generator(1.0);
        CGAL::cpp11::copy_n( pts_generator, 100, std::back_inserter(m_incrementalPts) );
      } else {
        for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
            vit < m_pScene->m_vhArray.end(); ++vit) {
          m_incrementalPts.push_back( (*vit)->point() );
        }//end-for
        // erase existing vertices
        initClean();
      }//end-if-pts
      // sorts points in a way that improves space locality
      CGAL::spatial_sort( m_incrementalPts.begin(), m_incrementalPts.end() );
      // set the current to "hightlight the new point"
      m_curStep = INIT;
    }/* else resume play */

    // set up the timer
    m_pTimer->start(1000);
  } else { // pause
    m_pTimer->stop();
  }

  // redraw
  updateGL();
}

void Viewer::stopIncremental() {
  if( !m_incrementalPts.isEmpty() ) {
    // will call toggleIncremental to stop the timer
    emit( stopIncAnimation() );

    // insert the rest points
    for(QList<Point_3>::iterator pit=m_incrementalPts.begin();
        pit < m_incrementalPts.end(); ++pit) {
      Vertex_handle hint;
      if( m_pScene->m_vhArray.isEmpty() ) {
        hint = m_pScene->m_dt.insert( *pit );
      } else {
        hint = m_pScene->m_vhArray.last();
        hint = m_pScene->m_dt.insert( *pit, hint );
      }
      m_pScene->m_vhArray.push_back( hint );
    }
    m_incrementalPts.clear();
  }

  // redraw
  updateGL();
}

void Viewer::incremental_insert() {
  Vertex_handle hint;
  if( !m_incrementalPts.isEmpty() ) {
  	switch( m_curStep ) {
  	case INIT:  // end of INIT: get the next-to-insert point
      m_curIncPt = m_incrementalPts.at(0);
      m_curStep = NEWPT;
  	  break;
    case NEWPT:  // end of NEWPT: locate the cell containing the point
      if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 ) {
        computeConflict( m_curIncPt );
        m_curStep = CELL;
      }
      else {
        // popup the first point and insert it
        m_curIncPt = m_incrementalPts.takeFirst();
        if( m_pScene->m_vhArray.isEmpty() ) {
          hint = m_pScene->m_dt.insert( m_curIncPt );
        }
        else {
          hint = m_pScene->m_vhArray.last();
          hint = m_pScene->m_dt.insert( m_curIncPt, hint );
        }
        m_pScene->m_vhArray.push_back( hint );
        m_curStep = INIT;
      }
      break;
    case CELL:  // end of CELL: compute the conflict region
      m_curStep = CONFLICT;
      break;
    case CONFLICT:  // end of CONFLICT: do the insertion and go back to INIT
      // popup the first point and insert it
      m_curIncPt = m_incrementalPts.takeFirst();
      if( m_pScene->m_vhArray.isEmpty() ) {
        hint = m_pScene->m_dt.insert( m_curIncPt );
      } else {
        hint = m_pScene->m_vhArray.last();
        hint = m_pScene->m_dt.insert( m_curIncPt, hint );
      }
      m_pScene->m_vhArray.push_back( hint );
      m_curStep = INIT;
      break;
  	}//end-of-switch
  } else {
    /* if finished, then start over */
    for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
        vit < m_pScene->m_vhArray.end(); ++vit) {
      m_incrementalPts.push_back( (*vit)->point() );
    }//end-for
    // erase existing vertices
    initClean();
    // set the current to "hightlight the new point"
    m_curStep = INIT;
  }

  // redraw
  updateGL();
}
