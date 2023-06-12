
#version 330

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 normal;

//out vec4 vCol;
//out vec3 vpos;
flat out vec3 vNormal;

uniform mat4 MVP; 

void main()
{
	//vpos = pos;
  vNormal = normal;
	gl_Position = MVP * vec4(pos.xyz, 1);
  //gl_Position = vec4(pos.xyz, 1);
}