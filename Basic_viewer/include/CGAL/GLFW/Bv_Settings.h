#pragma once

/*************OPENGL WINDOW PARAMS*************/

#ifndef WINDOW_WIDTH_INIT
#define WINDOW_WIDTH_INIT 500
#endif

#ifndef WINDOW_HEIGHT_INIT
#define WINDOW_HEIGHT_INIT 450
#endif

#ifndef WINDOW_SAMPLES
#define WINDOW_SAMPLES 4
#endif

/*********************************************/
#ifndef CLIPPING_PLANE_RENDERING_TRANSPARENCY
#define CLIPPING_PLANE_RENDERING_TRANSPARENCY 0.5f
#endif


#ifndef SIZE_POINTS
#define SIZE_POINTS 7.0f
#endif

#ifndef SIZE_EDGES
#define SIZE_EDGES  3.1f
#endif

#ifndef SIZE_RAYS
#define SIZE_RAYS   3.1f
#endif

#ifndef SIZE_LINES
#define SIZE_LINES  3.1f
#endif


#ifndef FACES_MONO_COLOR
#define FACES_MONO_COLOR    {60, 60, 200}
#endif

#ifndef VERTICES_MONO_COLOR
#define VERTICES_MONO_COLOR {200, 60, 60}
#endif

#ifndef EDGES_MONO_COLOR
#define EDGES_MONO_COLOR    {0, 0, 0}
#endif

#ifndef RAYS_MONO_COLOR
#define RAYS_MONO_COLOR     {0, 0, 0}
#endif

#ifndef LINES_MONO_COLOR
#define LINES_MONO_COLOR    {0, 0, 0}
#endif


#ifndef LIGHT_POSITION
#define LIGHT_POSITION {0.0f, 0.0f, 0.0f, 0.0f}
#endif


#ifndef AMBIENT_COLOR
#define AMBIENT_COLOR  {0.6f, 0.5f, 0.5f, 1.f}
#endif

#ifndef DIFFUSE_COLOR
#define DIFFUSE_COLOR  {0.9f, 0.9f, 0.9f, 1.0f}
#endif

#ifndef SPECULAR_COLOR
#define SPECULAR_COLOR {0.0f, 0.0f, 0.0f, 1.0f}
#endif

#ifndef SHININESS
#define SHININESS 0.5f
#endif

#ifndef CAM_MOVE_SPEED
#define CAM_MOVE_SPEED 0.4f 
#endif

#ifndef CAM_ROT_SPEED
#define CAM_ROT_SPEED 0.05f 
#endif

#ifndef CLIPPING_PLANE_MOVE_SPEED
#define CLIPPING_PLANE_MOVE_SPEED 0.04f 
#endif

#ifndef CLIPPING_PLANE_ROT_SPEED
#define CLIPPING_PLANE_ROT_SPEED 4.f 
#endif


#ifndef SCENE_ROT_SPEED
#define SCENE_ROT_SPEED 0.5f
#endif