///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CViewport                                                     //
//                                                                       //
//  Viewport class representing the image plane and corresponding window.//
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _VIEWPORT_
#define _VIEWPORT_

class CViewport
{

private :

protected:

 int originX, originY;
 int sizeX, sizeY;

public :

  // Constructors
  CViewport() { }
  CViewport(const int xRes, const int yRes)
    { sizeX = xRes; sizeY = yRes; }
  CViewport(const CViewport &v)  { Set(v); }
  CViewport(const CViewport *pV) { Set(pV); }

  virtual ~CViewport() { }

  // Data setting
  void Clear()
    { SetSize(640,480); SetOrigin(0,0); }
  void Set(const CViewport& v);
  void Set(const CViewport *pV);

  void SetOrigin(const int x, const int y)
    { originX = x; originY = y; }
  void SetSize(const int xRes, const int yRes)
    { sizeX = xRes; sizeY = yRes; }

  // Per coordinate (explicit inline functions)
  void Width(int w)  { sizeX = w; }
  void Height(int h) { sizeY = h; }
  void xRes(int x)   { sizeX = x; }
  void yRes(int y)   { sizeY = y; }
  void oX(int ox)    { originX = ox; }
  void oY(int oy)    { originY = oy; }

  // Data access (explicit inline functions)
  void GetSize(float& x, float& y) 
    { x = (float)sizeX; y = (float)sizeY; }
  float GetAspectRatio() const
    { return (float)sizeX / (float)sizeY; }

  // Per coordinate data access (explicit inline functions)
  int Width()   const { return sizeX; }
  int Height()  const { return sizeY; }
  int xRes()    const { return sizeX; }
  int yRes()    const { return sizeY; }
  int oX()      const { return originX; }
  int oY()      const { return originY; }

  // Debug
  void Trace() const;

  // OpenGL
  void glDraw() const;
};

#endif // _VIEWPORT_

