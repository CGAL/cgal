
#version 330

layout (location = 0) in vec3 a_pos;
layout (location = 1) in vec3 a_normal;


uniform mat4 u_mvp; 

void main()
{
	gl_Position = u_mvp * vec4(a_pos.xyz, 1);
}