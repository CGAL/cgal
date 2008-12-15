// PoissonView.cpp : implementation of the CPoissonView class
//

#include "stdafx.h"
#include "Poisson.h"
#include "PoissonDoc.h"
#include "PoissonView.h"

#include <cassert>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CPoissonView

IMPLEMENT_DYNCREATE(CPoissonView, CView)

BEGIN_MESSAGE_MAP(CPoissonView, CView)
    ON_WM_CREATE()
    ON_WM_PAINT()
    ON_WM_LBUTTONDOWN()
    ON_WM_LBUTTONUP()
    ON_WM_MBUTTONDOWN()
    ON_WM_MBUTTONUP()
    ON_WM_RBUTTONDOWN()
    ON_WM_RBUTTONUP()
    ON_WM_MOUSEMOVE()
    ON_WM_MOUSEWHEEL()
    ON_WM_SIZE()
    ON_COMMAND(ID_RENDER_DELAUNAYEDGES, OnRenderDelaunayedges)
    ON_COMMAND(ID_RENDER_CONTOUR, OnRenderContour)
    ON_UPDATE_COMMAND_UI(ID_RENDER_POINTS, OnUpdateRenderPoints)
    ON_COMMAND(ID_RENDER_POINTS, OnRenderPoints)
    ON_UPDATE_COMMAND_UI(ID_RENDER_DELAUNAYEDGES, OnUpdateRenderDelaunayedges)
    ON_UPDATE_COMMAND_UI(ID_RENDER_CONTOUR, OnUpdateRenderContour)
    ON_COMMAND(ID_RENDER_NORMALS, OnRenderNormals)
    ON_UPDATE_COMMAND_UI(ID_RENDER_NORMALS, OnUpdateRenderNormals)
    ON_COMMAND(ID_RENDER_SURFACE, OnRenderSurface)
    ON_UPDATE_COMMAND_UI(ID_RENDER_SURFACE, OnUpdateRenderSurface)
    ON_COMMAND(ID_ARCBALL_RESET, OnArcballReset)
    ON_COMMAND(ID_RENDER_ARCBALL, OnRenderArcball)
    ON_UPDATE_COMMAND_UI(ID_RENDER_ARCBALL, OnUpdateRenderArcball)
    ON_COMMAND(ID_RENDER_ORIGINAL_NORMALS, OnRenderOriginalnormals)
    ON_UPDATE_COMMAND_UI(ID_RENDER_ORIGINAL_NORMALS, OnUpdateRenderOriginalnormals)
END_MESSAGE_MAP()

// CPoissonView construction/destruction

CPoissonView::CPoissonView()
{
  // OpenGL
  m_hGLContext = NULL;
  m_GLPixelIndex = 0;

  // mouse
  m_LeftButtonDown = false;
  m_MiddleButtonDown = false;
  m_RightButtonDown = false;
  m_moving = false;

  // arcballl
  m_pArcball = NULL;
  first_paint = true;

  // options
  m_view_vertices = true;
  m_view_delaunay_edges = false;
  m_view_contour = true;
  m_view_normals = true; 
  m_view_original_normals = false; 
  m_view_surface = true;
  m_view_arcball = false;
}

CPoissonView::~CPoissonView()
{
  delete m_pArcball;
}

BOOL CPoissonView::PreCreateWindow(CREATESTRUCT& cs)
{
  // TODO: Modify the Window class or styles here by modifying
  //  the CREATESTRUCT cs

  return CView::PreCreateWindow(cs);
}

// CPoissonView drawing

void CPoissonView::OnDraw(CDC* /*pDC*/)
{
  CPoissonDoc* pDoc = GetDocument();
  ASSERT_VALID(pDoc);
  if (!pDoc)
    return;

  // TODO: add draw code for native data here
}


// CPoissonView diagnostics

#ifdef _DEBUG
void CPoissonView::AssertValid() const
{
  CView::AssertValid();
}

void CPoissonView::Dump(CDumpContext& dc) const
{
  CView::Dump(dc);
}

CPoissonDoc* CPoissonView::GetDocument() const // non-debug version is inline
{
  ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CPoissonDoc)));
  return (CPoissonDoc*)m_pDocument;
}
#endif 


//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CPoissonView::SetWindowPixelFormat(HDC hDC)
{
  PIXELFORMATDESCRIPTOR pixelDesc;
  
  pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
  pixelDesc.nVersion = 1;
  
  pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
    PFD_DOUBLEBUFFER | PFD_STEREO_DONTCARE;
  
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
  
  if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
    return FALSE;

  return TRUE;
}

//********************************************
// CreateViewGLContext
// Create an OpenGL rendering context
//********************************************
BOOL CPoissonView::CreateViewGLContext(HDC hDC)
{
  m_hGLContext = wglCreateContext(hDC);

  if(m_hGLContext==NULL)
    return FALSE;

  if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
    return FALSE;

  return TRUE;  
}


// CPoissonView message handlers

int CPoissonView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
  if(CView::OnCreate(lpCreateStruct) == -1)
    return -1;

  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);

  if(SetWindowPixelFormat(hDC)==FALSE)
    return -1;
  
  if(CreateViewGLContext(hDC)==FALSE)
    return -1;

  // Activate Z-buffer
  glEnable(GL_DEPTH_TEST);
  
  // Background color
  glClearColor(1.0f,1.0f,1.0f,1.0f);

  // activate lighting
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
  glLightfv(GL_LIGHT0,GL_POSITION,lpos);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

  ::glEnable(GL_NORMALIZE);

  ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  m_pArcball = new Arcball(100,100,1.0f,0.0f,0.0f,0.0f);

  return 0;
}


void CPoissonView::OnPaint()
{
  // device context for painting
  CPaintDC dc(this); 

  // model is stored in Document
  CPoissonDoc *pDoc = GetDocument();

  // useful in multidoc templates
  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);
  wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  assert(pDoc->points()->begin() != pDoc->points()->end());

  // Scene's region of interest (= bounding sphere minus outliers)
  Sphere region_of_interest;
  if (pDoc->edit_mode() == CPoissonDoc::POINT_SET)
    region_of_interest = pDoc->points()->region_of_interest();
  else if (pDoc->edit_mode() == CPoissonDoc::POISSON)
    region_of_interest = pDoc->poisson_function()->region_of_interest();
  else if (pDoc->edit_mode() == CPoissonDoc::APSS)
    region_of_interest = pDoc->apss_function()->region_of_interest();

  // Do points have normals?
  bool points_have_normals = (pDoc->points()->begin()->normal() != CGAL::NULL_VECTOR);
  bool normals_are_oriented = pDoc->points()->begin()->normal().is_oriented();
  bool points_have_original_normals = (pDoc->points()->begin()->original_normal() != CGAL::NULL_VECTOR);

  if(first_paint)
  {
    // allocate new arcball
    CRect rect;
    GetClientRect(&rect);
    int width = rect.Width();
    int height = rect.Height();
    float size = (float)sqrt(region_of_interest.squared_radius());
    float cx = (float)region_of_interest.center().x();
    float cy = (float)region_of_interest.center().y();
    float cz = (float)region_of_interest.center().z();
    delete m_pArcball;
    m_pArcball = new Arcball(width,height,size,cx,cy,cz);
    first_paint = false;
  }

  // shading option
  glShadeModel(GL_FLAT);

  // test
  glMatrixMode(GL_MODELVIEW);

  // drawing
  glPushMatrix();
  
    ::glDisable(GL_LIGHTING);

    // Draw arcball's "scene" sphere
    if ( m_view_arcball && (m_KeyboardFlags & MK_SHIFT) )
    {
        m_pArcball->setCamera_no_offset_no_scale();
        m_pArcball->draw_sphere_scene();
    }

    // setup viewpoint from current arcball
    m_pArcball->setCamera();

    // Draw arcball's "object" sphere
    if ( m_view_arcball && !(m_KeyboardFlags & MK_SHIFT) )
        m_pArcball->draw_sphere();

    // draw vertices
    if(m_view_vertices)
    {
      if (pDoc->edit_mode() == CPoissonDoc::POINT_SET || pDoc->edit_mode() == CPoissonDoc::APSS)
        pDoc->points()->gl_draw_vertices(0,0,0 /*color*/, 2.0f /*size*/);
      else if (pDoc->edit_mode() == CPoissonDoc::POISSON)
        pDoc->poisson_function()->triangulation().gl_draw_delaunay_vertices(0,0,0 /*color*/, 2.0f /*size*/);
    }

    // draw Delaunay edges
    if(m_view_delaunay_edges && pDoc->edit_mode() == CPoissonDoc::POISSON)
        pDoc->poisson_function()->triangulation().gl_draw_delaunay_edges(0,0,0 /*color*/, 1.0f /*size*/);

    // draw normals
    if(m_view_normals && points_have_normals)
    {
      float length = (float)sqrt(region_of_interest.squared_radius() / 5000.0f);
      if (pDoc->edit_mode() == CPoissonDoc::POINT_SET || pDoc->edit_mode() == CPoissonDoc::APSS)
        pDoc->points()->gl_draw_normals(0,255,0 /*color*/, length);
      else if (pDoc->edit_mode() == CPoissonDoc::POISSON)
        pDoc->poisson_function()->triangulation().gl_draw_normals(0,255,0 /*color*/, length);
    }
    if(m_view_original_normals && points_have_original_normals)
    {
      float length = (float)sqrt(region_of_interest.squared_radius() / 5000.0f);
      if (pDoc->edit_mode() == CPoissonDoc::POINT_SET || pDoc->edit_mode() == CPoissonDoc::APSS)
        pDoc->points()->gl_draw_original_normals(0,0,255 /*color*/, length, 0.5 /*width*/);
    }

    // draw surface reconstructed by marching tet
    if(m_view_contour && pDoc->edit_mode() == CPoissonDoc::POISSON)
    {
      static GLfloat grey[4] = {0.9f, 0.9f, 0.9f, 0.0f };
      ::glEnable(GL_LIGHTING);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,grey);
      ::glColor3ub(0,0,0);
      pDoc->marching_tet_countour()->gl_draw_surface();
    }

    // draw surface reconstructed by Surface Mesher
    if(m_view_surface && 
       (pDoc->edit_mode() == CPoissonDoc::POISSON || pDoc->edit_mode() == CPoissonDoc::APSS))
    {
      ::glEnable(GL_LIGHTING);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      ::glColor3ub(100,100,255);
      ::glEnable(GL_POLYGON_OFFSET_FILL);
      ::glPolygonOffset(3.0f,-3.0f);
      pDoc->surface_mesher_surface()->gl_draw_surface();

      ::glDisable(GL_LIGHTING);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      pDoc->surface_mesher_surface()->gl_draw_surface();
    }

  glPopMatrix();

  // Double buffer
  SwapBuffers(dc.m_ps.hdc);
  glFlush();

  // release hdc
  ::ReleaseDC(hWnd,hDC);
}

void CPoissonView::OnLButtonDown(UINT nFlags, CPoint point)
{
  SetCapture();
  m_LeftButtonDown = TRUE;
  m_prev_mouse_pos = point;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnLButtonDown(nFlags, point);
}

void CPoissonView::OnLButtonUp(UINT nFlags, CPoint point)
{
  ReleaseCapture();
  m_LeftButtonDown = FALSE;
  m_moving = false;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnLButtonUp(nFlags, point);
}

void CPoissonView::OnMButtonDown(UINT nFlags, CPoint point)
{
  SetCapture();
  m_MiddleButtonDown = TRUE;
  m_prev_mouse_pos = point;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnMButtonDown(nFlags, point);
}

void CPoissonView::OnMButtonUp(UINT nFlags, CPoint point)
{
  ReleaseCapture();
  m_MiddleButtonDown = FALSE;
  m_moving = false;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnMButtonUp(nFlags, point);
}

void CPoissonView::OnRButtonDown(UINT nFlags, CPoint point)
{
  SetCapture();
  m_RightButtonDown = TRUE;
  m_prev_mouse_pos = point;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnRButtonDown(nFlags, point);
}

void CPoissonView::OnRButtonUp(UINT nFlags, CPoint point)
{
  ReleaseCapture();
  m_RightButtonDown = FALSE;
  m_moving = false;
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
    m_KeyboardFlags = nFlags;
    InvalidateRect(NULL,FALSE);
  }
  CView::OnRButtonUp(nFlags, point);
}

// Mouse move callback
void CPoissonView::OnMouseMove(UINT nFlags, CPoint point)
{
  if (m_KeyboardFlags != nFlags) // update keyboard status and redraw arcball sphere
  {
      m_KeyboardFlags = nFlags;
      InvalidateRect(NULL,FALSE);
  }

  // Drag with middle button (or both left and right buttons): zoom
  if( m_MiddleButtonDown || (m_LeftButtonDown && m_RightButtonDown) )
  {
      m_moving = true;
      int dy = point.y - m_prev_mouse_pos.y;
      float factor = (float)(1.0 - dy * 0.01);
      m_pArcball->zoom(factor);
      m_prev_mouse_pos = point;
      InvalidateRect(NULL,FALSE);
  }
  // Drag with left button: rotate
  else if(m_LeftButtonDown)
  {
      m_moving = true;
      if (m_KeyboardFlags & MK_SHIFT)
          m_pArcball->rotate_scene(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
      else
          m_pArcball->rotate_obj(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
      m_prev_mouse_pos = point;
      InvalidateRect(NULL,FALSE);
  }
  // Drag withn right button: translate
  else if(m_RightButtonDown)
  {
      m_moving = true;
      m_pArcball->translate(m_prev_mouse_pos.x,m_prev_mouse_pos.y,point.x,point.y);
      m_prev_mouse_pos = point;
      InvalidateRect(NULL,FALSE);
  }
  else
      m_moving = false;
        
  CView::OnMouseMove(nFlags, point);
}

// Mouse wheel callback: zoom
BOOL CPoissonView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
  m_KeyboardFlags = nFlags;
  float nb = (float)zDelta / 120.0f;
  float factor = (float)(1.0 + 0.3 * nb);
  m_pArcball->zoom(factor);
  InvalidateRect(NULL,FALSE);
  return CView::OnMouseWheel(nFlags, zDelta, pt);
}


void CPoissonView::OnSize(UINT nType, int cx, int cy)
{
  CView::OnSize(nType, cx, cy);

  HWND hWnd = GetSafeHwnd();
  HDC hDC = ::GetDC(hWnd);
  wglMakeCurrent(hDC,m_hGLContext);

  glMatrixMode(GL_MODELVIEW);
  m_pArcball->resize(cx,cy);
  glDrawBuffer(GL_BACK);
  ::ReleaseDC(hWnd,hDC);
}

void CPoissonView::OnRenderDelaunayedges()
{
  m_view_delaunay_edges = !m_view_delaunay_edges;
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderDelaunayedges(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_delaunay_edges);
}

void CPoissonView::OnRenderContour()
{
  m_view_contour = !m_view_contour;
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderContour(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_contour);
}

void CPoissonView::OnRenderPoints()
{
  m_view_vertices = !m_view_vertices;
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderPoints(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_vertices);
}

void CPoissonView::OnRenderNormals()
{
  m_view_normals = !m_view_normals; 
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderNormals(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_normals);
}

void CPoissonView::OnRenderOriginalnormals()
{
  m_view_original_normals = !m_view_original_normals; 
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderOriginalnormals(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_original_normals);
}

void CPoissonView::OnRenderSurface()
{
  m_view_surface = !m_view_surface; 
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderSurface(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_surface);
}

void CPoissonView::OnArcballReset()
{
  m_pArcball->reset();
  InvalidateRect(NULL,FALSE);
}

void CPoissonView::OnRenderArcball()
{
  m_view_arcball = !m_view_arcball; 
  InvalidateRect(NULL,FALSE);
}
void CPoissonView::OnUpdateRenderArcball(CCmdUI *pCmdUI)
{
  pCmdUI->SetCheck(m_view_arcball);
}
