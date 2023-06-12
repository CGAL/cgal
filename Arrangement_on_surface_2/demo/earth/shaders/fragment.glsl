
#version 330

//in vec4 vCol;
flat in vec3 vNormal;

out vec4 color;

void main()
{
	const vec3 lightDir = normalize(vec3(1,.5,.5));

	//float c = clamp(dot(lightDir,triNormal), 0, 1);
	vec3 n = normalize(vNormal);
  float c = abs( dot(lightDir, n) );
	color = vec4(.2, .2,0,1) + 0.8*vec4(c,c,0,1);

	//color = vec4(1,1,0,1);
  //color = vCol;
}