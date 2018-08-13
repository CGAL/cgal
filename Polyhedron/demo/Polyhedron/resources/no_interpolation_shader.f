#version 150

in GS_OUT
{
  vec4 fP;
  vec3 fN;
  flat vec4 color[4];
  vec2 uv;
  flat vec4 prob[4];
} fs_in;

uniform bool is_two_side;
uniform vec4 light_pos;
uniform vec4 light_diff;
uniform vec4 light_spec;
uniform vec4 light_amb;
uniform float spec_power ;

out vec4 out_color;

void main(void)
{
  vec4 color;
  //find base color of pixel
  vec4 m1 = mix(fs_in.prob[0], fs_in.prob[1], fs_in.uv.y);
  vec4 m2 = mix(fs_in.prob[3], fs_in.prob[2], fs_in.uv.y);
  vec4 m3 = mix(m1,m2, fs_in.uv.x);
  float maximum =
    max(m3[0],
          max(m3[1],
            max(m3[2],m3[3])));
  if(maximum == m3[0])
    color = fs_in.color[0];
  else if(maximum == m3[1])
    color = fs_in.color[1];
  else if(maximum == m3[2])
    color = fs_in.color[2];
  else //(maximum == m3[3])
    color = fs_in.color[3];

  //compute and apply light effects
  vec3 L = light_pos.xyz - fs_in.fP.xyz;
  vec3 V = -fs_in.fP.xyz;

  vec3 N;
  if(fs_in.fN == vec3(0.0,0.0,0.0))
      N = vec3(0.0,0.0,0.0);
  else
      N = normalize(fs_in.fN);
  L = normalize(L);
  V = normalize(V);

  vec3 R = reflect(-L, N);
  vec4 diffuse;
  if(!is_two_side)
      diffuse = max(dot(N,L),0) * light_diff*color;
  else
      diffuse = max(abs(dot(N,L)),0) * light_diff*color;
  vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;
  out_color = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1.0);
}
