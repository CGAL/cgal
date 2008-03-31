// GyrovizView_2.cpp: implementation of the CGyrovizView_2 class
//

#include "stdafx.h"
#include "GyrovizDoc_2.h"
#include "Gyroviz.h"
#include "GyrovizView_2.h"
#include <limits.h>

#include <GL/gl.h>


#ifdef _DEBUG
#define new DEBUG_NEW
#endif




// CGyrovizView_2

IMPLEMENT_DYNCREATE(CGyrovizView_2, CView)

BEGIN_MESSAGE_MAP(CGyrovizView_2, CView)
  ON_WM_CREATE()
  ON_WM_SIZE()
  ON_WM_DESTROY()
  ON_WM_PAINT()
  ON_WM_LBUTTONUP()
  ON_WM_LBUTTONDOWN()
  ON_WM_MOUSEMOVE()
  ON_WM_RBUTTONUP()
  ON_WM_RBUTTONDOWN()
  ON_COMMAND(ID_RENDER_DELAUNAYEDGES, OnViewDelaunayedges)
  ON_UPDATE_COMMAND_UI(ID_RENDER_DELAUNAYEDGES, OnUpdateViewDelaunayedges)
  ON_WM_ERASEBKGND()
  ON_COMMAND(ID_RENDER_ZOOM, OnViewZoom)
  ON_UPDATE_COMMAND_UI(ID_RENDER_ZOOM, OnUpdateViewZoom)
  ON_COMMAND(ID_RENDER_ZOOM1, OnViewZoom1)
  ON_COMMAND(ID_RENDER_ZOOMIN, OnViewZoomin)
  ON_COMMAND(ID_RENDER_ZOOMOUT, OnViewZoomout)
  ON_WM_TIMER()
  ON_COMMAND(ID_PAN_LEFT, OnPanLeft)
  ON_COMMAND(ID_PAN_RIGHT, OnPanRight)
  ON_COMMAND(ID_PAN_TOP, OnPanTop)
  ON_COMMAND(ID_PAN_BOTTOM, OnPanBottom)
  ON_COMMAND(ID_RENDER_CONSTRAINEDEDGES, OnViewConstrainededges)
  ON_UPDATE_COMMAND_UI(ID_RENDER_CONSTRAINEDEDGES, OnUpdateViewConstrainededges)
  ON_WM_MOUSEWHEEL()
  ON_COMMAND(ID_RENDER_ALLCDTEDGES, OnViewAllcdtedges)
  ON_UPDATE_COMMAND_UI(ID_RENDER_ALLCDTEDGES, OnUpdateViewAllcdtedges)
  ON_COMMAND(ID_RENDER_GRAYSCALED, &CGyrovizView_2::OnRenderGrayscaled)
  ON_UPDATE_COMMAND_UI(ID_RENDER_GRAYSCALED, &CGyrovizView_2::OnUpdateRenderGrayscaled)
  ON_COMMAND(ID_RENDER_FILTERED, &CGyrovizView_2::OnRenderFiltered)
  ON_UPDATE_COMMAND_UI(ID_RENDER_FILTERED, &CGyrovizView_2::OnUpdateRenderFiltered)
  ON_COMMAND(ID_RENDER_SEGMENTED, &CGyrovizView_2::OnRenderSegmented)
  ON_UPDATE_COMMAND_UI(ID_RENDER_SEGMENTED, &CGyrovizView_2::OnUpdateRenderSegmented)
  ON_COMMAND(ID_RENDER_ORIGINALIMAGE, &CGyrovizView_2::OnRenderOriginalimage)
  ON_UPDATE_COMMAND_UI(ID_RENDER_ORIGINALIMAGE, &CGyrovizView_2::OnUpdateRenderOriginalimage)
  ON_COMMAND(ID_RENDER_POINTS, &CGyrovizView_2::OnRenderPoints)
  ON_UPDATE_COMMAND_UI(ID_RENDER_POINTS, &CGyrovizView_2::OnUpdateRenderPoints)
  ON_COMMAND(ID_RENDER_TRIANGLES, &CGyrovizView_2::OnRenderTriangles)
  ON_UPDATE_COMMAND_UI(ID_RENDER_TRIANGLES, &CGyrovizView_2::OnUpdateRenderTriangles)
  ON_COMMAND(ID_RENDER_BORDERVERTICES, &CGyrovizView_2::OnRenderBordervertices)
  ON_UPDATE_COMMAND_UI(ID_RENDER_BORDERVERTICES, &CGyrovizView_2::OnUpdateRenderBordervertices)
END_MESSAGE_MAP()





// CGyrovizView_2 construction/destruction

CGyrovizView_2::CGyrovizView_2()
{
  // mouse
  m_LeftButtonDown = false;
  m_RightButtonDown = false;
  m_Moving = false;
  m_has_moved = false;
  m_view_zoom_in = false;

  // Render menu

  // Image processing sub menu
  m_view_original = true;
  m_view_grayscale = false;
  m_view_filtering = false;
  m_view_segmentation = false;

  // CG sub menu
  m_view_all_2D_vertices = false;
  m_view_border_2D_vertices = false;
  m_view_2D_delaunay_triangulation = false;
  m_view_delaunay_edges = false;
  m_view_all_cdt_edges = false;
  m_view_constrained_edges = false;
}

CGyrovizView_2::~CGyrovizView_2()
{
}

BOOL CGyrovizView_2::PreCreateWindow(CREATESTRUCT& cs)
{
  // TODO: Modify the Window class or styles here by modifying
  //  the CREATESTRUCT cs

  return CView::PreCreateWindow(cs);
}

// CGyrovizView_2 drawing

void CGyrovizView_2::OnDraw(CDC* /*pDC*/)
{
  CGyrovizDoc_2* pDoc = GetDocument();
  ASSERT_VALID(pDoc);
  if (!pDoc)
    return;
  // TODO: add draw code for native data here
}


// CGyrovizView_2 diagnostics
#ifdef _DEBUG
void CGyrovizView_2::AssertValid() const
{
  CView::AssertValid();
}

void CGyrovizView_2::Dump(CDumpContext& dc) const
{
  CView::Dump(dc);
}

CGyrovizDoc_2* CGyrovizView_2::GetDocument() const // non-debug version is inline
{
  ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGyrovizDoc_2)));
  return (CGyrovizDoc_2*)m_pDocument;
}
#endif //_DEBUG

BOOL CGyrovizView_2::SetWindowPixelFormat(HDC hDC)
{
  PIXELFORMATDESCRIPTOR pixelDesc;

  pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
  pixelDesc.nVersion = 1;
  pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | 
    PFD_SUPPORT_OPENGL |
    PFD_DOUBLEBUFFER |
    PFD_STEREO_DONTCARE;
  pixelDesc.iPixelType = PFD_TYPE_RGBA;
  pixelDesc.cColorBits = 32;
  pixelDesc.cRedBits = 8;
  pixelDesc.cRedShift = 16;
  pixelDesc.cGreenBits = 8;
  pixelDesc.cGreenShift = 8;
  pixelDesc.cBlueBits = 8;
  pixelDesc.cBlueShift = 0;
  pixelDesc.cAlphaBits = 0;
  pixelDesc.cAlphaShift = 0;
  pixelDesc.cAccumBits = 64;
  pixelDesc.cAccumRedBits = 16;
  pixelDesc.cAccumGreenBits = 16;
  pixelDesc.cAccumBlueBits = 16;
  pixelDesc.cAccumAlphaBits = 0;
  pixelDesc.cDepthBits = 32;
  pixelDesc.cStencilBits = 8;
  pixelDesc.cAuxBuffers = 0;
  pixelDesc.iLayerType = PFD_MAIN_PLANE;
  pixelDesc.bReserved = 0;
  pixelDesc.dwLayerMask = 0;
  pixelDesc.dwVisibleMask = 0;
  pixelDesc.dwDamageMask = 0;

  m_GLPixelIndex = ChoosePixelFormat(hDC,&pixelDesc);
  if(m_GLPixelIndex == 0) // Choose default
  {
    m_GLPixelIndex = 1;
    if(DescribePixelFormat(hDC,m_GLPixelIndex,
      sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
      return FALSE;
  }

  if(SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc)==FALSE)
    return FALSE;

  return TRUE;
}


BOOL CGyrovizView_2::CreateViewGLContext(HDC hDC)
{
  m_hGLContext = wglCreateContext(hDC);
  if(m_hGLContext==NULL)
    return FALSE;
  if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
    return FALSE;
  return TRUE;
}
int CGyrovizView_2::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
  if (CView::OnCreate(lpCreateStruct) == -1)
    return -1;

  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);

  if(SetWindowPixelFormat(hDC)==FALSE)
    return 0;
  if(CreateViewGLContext(hDC)==FALSE)
    return 0;

  COLORREF color = ::GetSysColor(COLOR_3DFACE);
  glClearColor((float)GetRValue(color)/255.0f,
    (float)GetGValue(color)/255.0f,
    (float)GetBValue(color)/255.0f,
    1.0);
  glClearColor(1.0f,1.0f,1.0f,1.0f);

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  return 0;
}

void CGyrovizView_2::OnSize(UINT nType, int cx, int cy)
{
  resize(nType,cx,cy);
}

// On window resize
void CGyrovizView_2::resize(UINT nType, int cx, int cy)
{
  CView::OnSize(nType, cx, cy);

  // Set viewport size
  glViewport(0,0,(GLsizei)cx,(GLsizei)cy);

  // Get displayed image's size
  CGyrovizDoc_2* pDoc = GetDocument();
  double image_dim_x = pDoc->m_cimg_interm_image.dimx(), 
    image_dim_y = pDoc->m_cimg_interm_image.dimy();

  // Set camera matrix
  double dx = (double)cx;
  double dy = (double)cy;
  if(dx/dy  > image_dim_x/image_dim_y)
  {
    double ratio = (dx/dy) / (image_dim_x/image_dim_y); // dx/dy
    m_xmin = image_dim_x*(1.0-ratio)/2.0; //0.5-ratio/2
    m_xmax = image_dim_x*(1.0+ratio)/2.0; //0.5+ratio/2;
    m_ymin = 0;
    m_ymax = image_dim_y;//1
  }

  else
  {
    double ratio = (dy/dx) / (image_dim_y/image_dim_x); //dy/dx;
    m_xmin = 0;
    m_xmax = image_dim_x;//1;
    m_ymin = image_dim_y*(1.0-ratio)/2.0;//0.5-ratio/2;
    m_ymax = image_dim_y*(1.0+ratio)/2.0;//0.5+ratio/2;
  }
  // glMatrixMode(GL_PROJECTION);
  // glLoadIdentity();
  // gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  // // invert the y axis, down is positive
  //glScalef(1, -1, 1);
  //// move the origin from the bottom left corner
  //// to the upper left corner
  //glTranslatef(0, -image_dim_y, 0);
  // 
  // glMatrixMode(GL_MODELVIEW);
  // glLoadIdentity();  
  // glDrawBuffer(GL_BACK);

  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnDestroy()
{
  if(wglGetCurrentContext() != NULL)
    wglMakeCurrent(NULL,NULL);

  if(m_hGLContext != NULL)
  {
    wglDeleteContext(m_hGLContext);
    m_hGLContext = NULL;
  }
  CView::OnDestroy();
}



void CGyrovizView_2::OnPaint()
{
  // device context for painting
  CPaintDC dc(this);

  CGyrovizDoc_2* pDoc = GetDocument();

  // useful in multidoc templates
  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);
  wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);

  ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set camera matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  // invert the y axis, down is positive
  glScalef(1, -1, 1);
  // move the origin from the bottom left corner
  // to the upper left corner
  glTranslatef(0, (float)-pDoc->m_cimg_interm_image.dimy(), 0);

  // Init world (model) matrix 
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();  
  glDrawBuffer(GL_BACK);



  // if(m_view_delaunay_edges)
  //   pDoc->m_cdt.gl_draw_unconstrained_edges(1.0f,0,0,0);


  //if(m_view_all_cdt_edges)
  //	pDoc->m_cdt.gl_draw_all_edges(1.0f,5,5,5);


  // if(m_view_constrained_edges)
  //	pDoc->m_cdt.gl_draw_constrained_edges(2.5f,180,0,0);


  if(m_view_original)
  {
    // draw original image
    glPolygonMode(GL_BACK,GL_FILL);
    ::glColor3ub(255,255,255);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pDoc->m_cimg_interm_image.dimx(), pDoc->m_cimg_interm_image.dimy(), 0, GL_RGB, GL_UNSIGNED_BYTE, pDoc->m_original_image);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);
    glTexCoord2d(0,0); glVertex2d(0,0);
    glTexCoord2d(1,0); glVertex2d(pDoc->m_cimg_interm_image.dimx(),0);
    glTexCoord2d(1,1); glVertex2d(pDoc->m_cimg_interm_image.dimx(),pDoc->m_cimg_interm_image.dimy());
    glTexCoord2d(0,1); glVertex2d(0,pDoc->m_cimg_interm_image.dimy());
    glEnd();
  }


  if(m_view_grayscale)
  {
    // draw grayscaled Image
    glPolygonMode(GL_BACK,GL_FILL);
    ::glColor3ub(255,255,255);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pDoc->m_cimg_gray_image.dimx(), pDoc->m_cimg_gray_image.dimy(), 0, GL_RGB, GL_UNSIGNED_BYTE, pDoc->m_grayscaled_image);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);
    glTexCoord2d(0,0); glVertex2d(0,0);
    glTexCoord2d(1,0); glVertex2d(pDoc->m_cimg_gray_image.dimx(),0);
    glTexCoord2d(1,1); glVertex2d(pDoc->m_cimg_gray_image.dimx(),pDoc->m_cimg_gray_image.dimy());
    glTexCoord2d(0,1); glVertex2d(0,pDoc->m_cimg_gray_image.dimy());
    glEnd();
  }


  if(m_view_filtering)
  {
    // draw filtered image Gauss 3x3 Mask
    glPolygonMode(GL_BACK,GL_FILL);
    ::glColor3ub(255,255,255);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pDoc->m_cimg_filt_image.dimx(), pDoc->m_cimg_filt_image.dimy(), 0, GL_RGB, GL_UNSIGNED_BYTE, pDoc->m_filtered_image);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);
    glTexCoord2d(0,0); glVertex2d(0,0);
    glTexCoord2d(1,0); glVertex2d(pDoc->m_cimg_filt_image.dimx(),0);
    glTexCoord2d(1,1); glVertex2d(pDoc->m_cimg_filt_image.dimx(),pDoc->m_cimg_filt_image.dimy());
    glTexCoord2d(0,1); glVertex2d(0,pDoc->m_cimg_filt_image.dimy());
    glEnd();
  }


  if(m_view_segmentation)
  {
    // draw segmented image Frei-Chen Mask
    glPolygonMode(GL_BACK,GL_FILL);
    ::glColor3ub(255,255,255);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pDoc->m_cimg_seg_image.dimx(), pDoc->m_cimg_seg_image.dimy(), 0, GL_RGB, GL_UNSIGNED_BYTE, pDoc->m_segmented_image);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);
    glTexCoord2d(0,0); glVertex2d(0,0);
    glTexCoord2d(1,0); glVertex2d(pDoc->m_cimg_seg_image.dimx(),0);
    glTexCoord2d(1,1); glVertex2d(pDoc->m_cimg_seg_image.dimx(),pDoc->m_cimg_seg_image.dimy());
    glTexCoord2d(0,1); glVertex2d(0,pDoc->m_cimg_seg_image.dimy());
    glEnd();
  }


  if(m_view_all_2D_vertices)
  {
    ::glColor3ub(0,0,255);
    pDoc->get_dt2().gl_draw_2D_vertices(0,0,255,4);
  }
  
  
 if(m_view_border_2D_vertices)
  {
    ::glColor3ub(255,0,0);
    pDoc->get_dt2().gl_draw_on_border_2D_vertices(255,0,0,4,pDoc->m_cimg_seg_image);
  }

  if(m_view_2D_delaunay_triangulation)
  {
    ::glColor3ub(0,255,0);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    pDoc->get_dt2().gl_draw_2D_delaunay_triangles(0,255,0);
  }

  if(m_view_zoom_in)
  {
    double m_xmin = min(m_corner2.x(),m_corner1.x());
    double m_xmax = max(m_corner2.x(),m_corner1.x());
    double m_ymin = min(m_corner2.y(),m_corner1.y());
    double m_ymax = max(m_corner2.y(),m_corner1.y());
    ::glColor3ub(0,100,0);
    ::glBegin(GL_LINE_LOOP);
    ::glVertex2d(m_xmin,m_ymin);
    ::glVertex2d(m_xmax,m_ymin);
    ::glVertex2d(m_xmax,m_ymax);
    ::glVertex2d(m_xmin,m_ymax);
    ::glEnd();
  }

  // Double buffer
  ::glFlush();
  SwapBuffers(dc.m_ps.hdc);
}

void CGyrovizView_2::OnLButtonDown(UINT nFlags, CPoint point)
{
  m_mouse_point = convert(point);
  CGyrovizDoc_2* pDoc = GetDocument();
  if(m_view_zoom_in)
  {
    m_corner1 = convert(point);
    m_corner2 = convert(point);
  }
  SetCapture();
  m_LeftButtonDown = TRUE;
  CView::OnLButtonDown(nFlags, point);
}
void CGyrovizView_2::OnRButtonDown(UINT nFlags, CPoint point)
{
  //SetCapture();
  //m_RightButtonDown = TRUE;

  // // snap vertex
  // CGyrovizDoc_2* pDoc = GetDocument();
  //pDoc->m_selected_vertex = pDoc->m_cdt.select_vertex(convert(point));
  //InvalidateRect(NULL,FALSE);

  // CView::OnRButtonDown(nFlags, point);
}

void CGyrovizView_2::OnLButtonUp(UINT nFlags, CPoint point)
{
  ReleaseCapture();
  m_LeftButtonDown = FALSE;

  CGyrovizDoc_2* pDoc = GetDocument();

  if(m_view_zoom_in)
  {
    m_corner2 = convert(point);
    double xmin = min(m_corner2.x(),m_corner1.x());
    double xmax = max(m_corner2.x(),m_corner1.x());
    double ymin = min(m_corner2.y(),m_corner1.y());
    double ymax = max(m_corner2.y(),m_corner1.y());
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double mx = 0.5*(xmax + xmin);
    double my = 0.5*(ymax + ymin);

    CRect rect;
    GetClientRect(rect);
    double w = (double)rect.Width();
    double h = (double)rect.Height();

    double nw,nh;
    if(dx/dy < w/h)
    {
      nw = w/h*dy;
      nh = dy;
    }
    else
    {
      nw = dx;
      nh = h/w*dx;
    }

    m_xmin = mx-nw/2;
    m_xmax = mx+nw/2;
    m_ymin = my-nh/2;
    m_ymax = my+nh/2;

    //::glMatrixMode(GL_PROJECTION);
    //::glLoadIdentity();
    //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
    //::glMatrixMode(GL_MODELVIEW);

    m_view_zoom_in = false;
    InvalidateRect(NULL,FALSE);
  }

  InvalidateRect(NULL,FALSE);

  CView::OnLButtonUp(nFlags, point);
}
void CGyrovizView_2::OnRButtonUp(UINT nFlags, CPoint point)
{
  ReleaseCapture();
  m_RightButtonDown = FALSE;

  // add seed point
  // CGyrovizDoc_2* pDoc = GetDocument();
  // pDoc->m_seeds.push_back(convert(point));

  InvalidateRect(NULL,FALSE);

  CView::OnRButtonUp(nFlags, point);
}

void CGyrovizView_2::OnMouseMove(UINT nFlags, 
                                 CPoint point)
{
  CGyrovizDoc_2* pDoc = GetDocument();

  if(m_view_zoom_in && m_LeftButtonDown)
  {
    m_corner2 = convert(point);
    InvalidateRect(NULL,FALSE);
  }
  else if(m_LeftButtonDown)
  {
    //if(m_drawing_mode)
    //{
    //  pDoc->m_component.push_back(convert(point));
    //  m_has_moved = true;
    // InvalidateRect(NULL,FALSE);
    //}
    //else
    //{
    // panning
    /*
    Point_2 p = convert(point);
    Vector_2 v = p - m_mouse_point;
    m_mouse_point = p;
    m_xmin += 1.0 * v.x();
    m_xmax += 1.0 * v.x();
    m_ymin += 1.0 * v.y();
    m_ymax += 1.0 * v.y();

    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();
    ::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
    ::glMatrixMode(GL_MODELVIEW);
    */

    //if(m_spray)
    //{
    //	pDoc->m_cdt.insert_spay(convert(point),pDoc->m_spray,0.01);
    //	pDoc->m_cdt.tag_faces_inside_outside(pDoc->m_seeds);
    //}
    //else
    //{
    //	// std insertion
    //	pDoc->m_cdt.insert(convert(point));
    //	pDoc->m_cdt.tag_faces_inside_outside(pDoc->m_seeds);
    //}
    InvalidateRect(NULL,FALSE);
  }
  else if(m_RightButtonDown && m_LeftButtonDown)
  {
  }
  else if(m_RightButtonDown)
  {
    //  if(pDoc->m_selected_vertex != NULL)
    //  {
    //	pDoc->m_cdt.remove(pDoc->m_selected_vertex);
    //		pDoc->m_selected_vertex = pDoc->m_cdt.insert(convert(point));
    //pDoc->m_cdt.tag_faces_inside_outside(pDoc->m_seeds);
    // 	InvalidateRect(NULL,FALSE);
    //  }
  }
  CView::OnMouseMove(nFlags, point);
}

// Convert 2D MFC point to OpenGL
Point_2 CGyrovizView_2::convert(const CPoint& point)
{
  CRect rect;
  GetClientRect(rect);
  FT cx = (FT)point.x;
  FT cy = (FT)point.y;
  FT w = (FT)rect.Width();
  FT h = (FT)rect.Height();
  double xrange = m_xmax-m_xmin;
  double yrange = m_ymax-m_ymin;
  FT x = m_xmin+cx/w*xrange;
  FT y = m_ymin+(1-cy/h)*yrange;
  return Point_2(x,y);
}


BOOL CGyrovizView_2::OnEraseBkgnd(CDC* pDC)
{
  return TRUE; // avoids flickering
}

void CGyrovizView_2::OnActivateView(BOOL bActivate, CView* pActivateView, CView* pDeactiveView)
{
  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);
  wglMakeCurrent(hDC,m_hGLContext);
  ::ReleaseDC(hWnd,hDC);
  CView::OnActivateView(bActivate, 
    pActivateView, 
    pDeactiveView);
}

void CGyrovizView_2::OnViewZoom()
{
  m_view_zoom_in = true;
  m_corner2 = Point_2(0,0);
  m_corner1 = Point_2(0,0);
}
void CGyrovizView_2::OnUpdateViewZoom(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_zoom_in);
}

void CGyrovizView_2::OnViewZoom1()
{
  CRect rect;
  GetClientRect(rect);
  resize(SIZE_RESTORED,rect.Width(),rect.Height());
}


void CGyrovizView_2::OnViewZoomin()
{
  double w = m_xmax - m_xmin;
  double h = m_ymax - m_ymin;
  double r = 0.01;

  m_xmin += w*r;
  m_xmax -= w*r;
  m_ymin += h*r;
  m_ymax -= h*r;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  m_view_zoom_in = false;
  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnViewZoomout()
{
  double w = m_xmax - m_xmin;
  double h = m_ymax - m_ymin;
  double r = 0.01;

  m_xmin -= w*r;
  m_xmax += w*r;
  m_ymin -= h*r;
  m_ymax += h*r;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  m_view_zoom_in = false;
  InvalidateRect(NULL,FALSE);
}

// zoom with mouse wheel
BOOL CGyrovizView_2::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
  double nb = (double)zDelta / 120.0;
  if(nb > 0)
    for(int i=0;i<10;i++)
      OnViewZoomin();
  else
    for(int i=0;i<10;i++)
      OnViewZoomout();
  return CView::OnMouseWheel(nFlags, zDelta, pt);
}


void CGyrovizView_2::OnInitialUpdate()
{
  CView::OnInitialUpdate();
  OnViewZoom1();
}


void CGyrovizView_2::OnPanLeft()
{
  double dx = m_xmax - m_xmin;
  double ratio = 0.01;
  m_xmin -= ratio * dx;
  m_xmax -= ratio * dx;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnPanRight()
{
  double dx = m_xmax - m_xmin;
  double ratio = 0.01;
  m_xmin += ratio * dx;
  m_xmax += ratio * dx;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnPanTop()
{
  double dy = m_ymax - m_ymin;
  double ratio = 0.01;
  m_ymin += ratio * dy;
  m_ymax += ratio * dy;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnPanBottom()
{
  double dy = m_ymax - m_ymin;
  double ratio = 0.01;
  m_ymin -= ratio * dy;
  m_ymax -= ratio * dy;

  //::glMatrixMode(GL_PROJECTION);
  //::glLoadIdentity();
  //::gluOrtho2D(m_xmin,m_xmax,m_ymin,m_ymax);
  //::glMatrixMode(GL_MODELVIEW);

  InvalidateRect(NULL,FALSE);
}


void CGyrovizView_2::OnViewDelaunayedges()
{
  m_view_delaunay_edges = !m_view_delaunay_edges;
  InvalidateRect(NULL,FALSE);
}
void CGyrovizView_2::OnUpdateViewDelaunayedges(CCmdUI *pCmdUI)
{
}
void CGyrovizView_2::OnViewAllcdtedges()
{
  m_view_all_cdt_edges = !m_view_all_cdt_edges;
  InvalidateRect(NULL,FALSE);
}


void CGyrovizView_2::OnUpdateViewAllcdtedges(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_all_cdt_edges);
}


void CGyrovizView_2::OnViewConstrainededges()
{
  m_view_constrained_edges = !m_view_constrained_edges;
  InvalidateRect(NULL,FALSE);
}
void CGyrovizView_2::OnUpdateViewConstrainededges(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_constrained_edges);
}


void CGyrovizView_2::OnRenderOriginalimage()
{
  m_view_original = !m_view_original;
  InvalidateRect(NULL,FALSE);
}
void CGyrovizView_2::OnUpdateRenderOriginalimage(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck( m_view_original);
}

void CGyrovizView_2::OnRenderGrayscaled()
{
  m_view_grayscale = !m_view_grayscale;
  InvalidateRect(NULL,FALSE);
}
void CGyrovizView_2::OnUpdateRenderGrayscaled(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck( m_view_grayscale);
}


void CGyrovizView_2::OnRenderFiltered()
{
  m_view_filtering = !m_view_filtering;
  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnUpdateRenderFiltered(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_filtering);
}

void CGyrovizView_2::OnRenderSegmented()
{
  m_view_segmentation = !m_view_segmentation;
  InvalidateRect(NULL,FALSE);
}


void CGyrovizView_2::OnUpdateRenderSegmented(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_segmentation);
}



void CGyrovizView_2::OnRenderPoints()
{
  m_view_all_2D_vertices = !m_view_all_2D_vertices;
  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnUpdateRenderPoints(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_all_2D_vertices);
}


void CGyrovizView_2::OnRenderBordervertices()
{

  m_view_border_2D_vertices = ! m_view_border_2D_vertices;
  InvalidateRect(NULL,FALSE);

}

void CGyrovizView_2::OnUpdateRenderBordervertices(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_border_2D_vertices);
}


void CGyrovizView_2::OnRenderTriangles()
{
  m_view_2D_delaunay_triangulation = !m_view_2D_delaunay_triangulation;
  InvalidateRect(NULL,FALSE);
}

void CGyrovizView_2::OnUpdateRenderTriangles(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_2D_delaunay_triangulation);
}

