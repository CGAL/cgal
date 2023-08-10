
#version 330

layout (location = 0) in vec3 a_pos;
layout (location = 1) in vec3 a_normal;

//out vec4 v_col;
out vec3 v_pos;
out vec3 v_normal;

uniform mat4 u_mvp; 

void main()
{
	v_pos = a_pos;
  v_normal = a_normal;
	gl_Position = u_mvp * vec4(a_pos.xyz, 1);
  //gl_Position = vec4(pos.xyz, 1);
}