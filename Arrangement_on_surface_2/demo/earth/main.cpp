

#include <stdio.h>
#include <iostream>
using namespace std;

#include <string>
//using string::literals;


#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <glm\gtc\type_ptr.hpp>

//const GLint width = 800, height = 600;
const GLint width = 800, height = 800;

GLuint vao, vbo, ibo, shader;
GLuint uniformMVP; // ModelViewProjection

bool direction = true;
float triOffset = 0;
float triMaxOffset = 0.7;
float triIncrement = 0.0005;


// vertex shader
const char* vShader = R"vs(
#version 330

layout (location = 0) in vec3 pos;
//out vec4 vCol;
out vec3 vpos;

uniform mat4 MVP; 

void main()
{
	vpos = pos;
	gl_Position = MVP * vec4(pos.xyz, 1);
}
)vs";


// GEOMETRY SHADER
// * I am using the geometry shader to compute the face-normals in the GPU on the fly
const char* gShader = R"gs(
#version 330

in vec3 vpos[];
out vec4 vCol;

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

void main()
{ 
	const vec3 lightDir = normalize(vec3(1,.5,.5));

	// compute the normal for the current triangle
	vec3 triNormal = normalize(cross(vpos[1]-vpos[0], vpos[2]-vpos[0]));
	float c = clamp(dot(lightDir,triNormal), 0, 1);
	vCol = vec4(.2, .2,0,1) + vec4(c,c,0,1);

	gl_Position = gl_in[0].gl_Position; EmitVertex();
	gl_Position = gl_in[1].gl_Position; EmitVertex();
	gl_Position = gl_in[2].gl_Position; EmitVertex();
	EndPrimitive();
}
)gs";


// FRAGMENT SHADER
static const char* fShader = R"fs(
#version 330

in vec4 vCol;
out vec4 color;

void main()
{
	color = vCol;  //vec4(1,0,0,1);
}

)fs";


void createBox()
{
	const float c = 2.5;
	GLfloat vertices[] = {
		// front
		-c, -c,  c,
		 c, -c,  c,
		 c,  c,  c,
		-c,  c,  c,
		// back
		-c, -c, -c,
		 c, -c, -c,
		 c,  c, -c,
		-c,  c, -c
	};

	GLuint indices[] = {
		// front
		0, 1, 2,
		2, 3, 0,
		// right
		1, 5, 6,
		6, 2, 1,
		// back
		7, 6, 5,
		5, 4, 7,
		// left
		4, 0, 3,
		3, 7, 4,
		// bottom
		4, 5, 1,
		1, 0, 4,
		// top
		3, 2, 6,
		6, 7, 3
	};

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	{
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		{
			glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

			GLint index = 0;
			glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(index);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);

	}
	glBindVertexArray(0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


#include <cmath>
#include <vector>
using namespace std;
using namespace glm;
#define M_PI       3.14159265358979323846


int numIndices = 0;
void createSphere(int numSlices, int numStacks, float r)
{
	numStacks = glm::max<int>(2, numStacks);
	vector<vec3> vertices, normals;

	// NORTH POLE
	vertices.push_back(vec3(0, 0, r));
	normals.push_back(vec3(0, 0, 1));

	// SOUTH POLE
	vertices.push_back(vec3(0, 0, -r));
	normals.push_back(vec3(0, 0, -1));
	int startingIndexOfMiddleVertices = vertices.size();

	for (int j = 1; j < numStacks; j++)
	{
		// Calculate the latitude (vertical angle) for the current stack
		float lat = M_PI * j / numStacks;
		float rxy = r * sin(lat);
		float z = r * cos(lat);

		for (int i = 0; i < numSlices; i++)
		{
			// Calculate the longitude (horizontal angle) for the current slice
			float lon = 2 * M_PI * i / numSlices;

			// Convert spherical coordinates to Cartesian coordinates
			float x = rxy * cos(lon);
			float y = rxy * sin(lon);

			auto p = vec3(x, y, z);
			auto n = p / length(p);
			vertices.push_back(p);
			normals.push_back(n);
		}
	}

	// add the indices for all trinagles
	vector<GLuint> indices;

	// NORTH CAP
	const int northVertexIndex = 0;
	const int northCapVertexIndexStart = startingIndexOfMiddleVertices;
	for (int i = 0; i < numSlices; i++)
	{
		indices.push_back(northVertexIndex);
		indices.push_back(northCapVertexIndexStart + i);
		indices.push_back(northCapVertexIndexStart + (i + 1) % numSlices);
	}

	// 0 = NORTH VERTEX
	// 1 = SOUTH VERTEX
	// [2, 2 + (numSlices-1)] = bottom vertices of the stack #1
	// [2+numSlices, 2 + (2*numSlices - 1)] = bottom vertices of the stack #2
	// ...
	// [2+(k-1)*numSlices, 2 + (k*numSlices -1) ] = bottom vertices of the stack #k
	// ..
	// [2+(numStacks-1)*numSlices, 2+(numStacks*numSlices-1)] = bottom vertices of the last stack (# numStacks)

	// SOUTH CAP
	const int southVertexIndex = 1;
	const int southCapIndexStart = startingIndexOfMiddleVertices + (numStacks - 2) * numSlices;
	for (int i = 0; i < numSlices; i++)
	{
		const auto vi0 = southVertexIndex;
		const auto vi1 = southCapIndexStart + i;
		const auto vi2 = southCapIndexStart + (i + 1) % numSlices;
		indices.push_back(vi2);
		indices.push_back(vi1);
		indices.push_back(vi0);
	}

	// MIDDLE TRIANGLES
	for (int k = 0; k < numStacks - 2; k++)
	{
		const int stackStartIndex = startingIndexOfMiddleVertices + k * numSlices;
		const int nextStackStartIndex = stackStartIndex + numSlices;
		for (int i = 0; i < numSlices; i++)
		{
			//int vi0 = stackStartIndex + i;
			//int vi1 = nextStackStartIndex + i;
			//int vi2 = nextStackStartIndex + (i + 1) % numSlices;
			//int vi3 = stackStartIndex + (i + 1) % numSlices;
			int vi0 = stackStartIndex + i;
			int vi1 = stackStartIndex + (i + 1) % numSlices;
			int vi2 = nextStackStartIndex + i;
			int vi3 = nextStackStartIndex + (i + 1) % numSlices;

			indices.push_back(vi0);
			indices.push_back(vi2);
			indices.push_back(vi1);
			//
			indices.push_back(vi2);
			indices.push_back(vi3);
			indices.push_back(vi1);
		}
	}

	numIndices = indices.size();


	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	{
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size(), reinterpret_cast<const void*>(indices.data()), GL_STATIC_DRAW);

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		{
			glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * vertices.size(), reinterpret_cast<const void*>(vertices.data()), GL_STATIC_DRAW);

			GLint index = 0;
			glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(index);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);

	}
	glBindVertexArray(0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void addShader(GLuint theProgram, const char* shaderCode, GLenum shaderType)
{
	GLuint theShader = glCreateShader(shaderType);

	const GLchar* theCode[] = { shaderCode };
	GLint codeLength[] = { strlen(shaderCode) };

	glShaderSource(theShader, 1, theCode, codeLength);
	glCompileShader(theShader);


	GLint result = 0;
	GLchar elog[1024] = { 0 };
	glGetShaderiv(theShader, GL_COMPILE_STATUS, &result);
	if (!result)
	{
		string shaderTypeName;
		switch (shaderType)
		{
		case GL_VERTEX_SHADER:   shaderTypeName = "VERTEX"; break;
		case GL_GEOMETRY_SHADER: shaderTypeName = "GEOMETRY"; break;
		case GL_FRAGMENT_SHADER: shaderTypeName = "FRAGMENT"; break;
		}
		glGetShaderInfoLog(theShader, sizeof(elog), NULL, elog);
		cout << "! error compiling the " << shaderTypeName << " shader:\n" << elog << endl;
		return;
	}

	glAttachShader(theProgram, theShader);
}

void compileShaders()
{
	shader = glCreateProgram();
	if (!shader) {
		cout << "error creating shader program!\n";
		return;
	}

	addShader(shader, vShader, GL_VERTEX_SHADER);
	addShader(shader, gShader, GL_GEOMETRY_SHADER);
	addShader(shader, fShader, GL_FRAGMENT_SHADER);

	GLint result = 0;
	GLchar elog[1024] = { 0 };

	glLinkProgram(shader);
	glGetProgramiv(shader, GL_LINK_STATUS, &result);
	if (!result)
	{
		glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
		cout << "! error linking program:\n" << elog << endl;
		return;
	}

	glValidateProgram(shader);
	glGetProgramiv(shader, GL_VALIDATE_STATUS, &result);
	if (!result)
	{
		glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
		cout << "! error validating program:\n" << elog << endl;
		return;
	}

	uniformMVP = glGetUniformLocation(shader, "MVP");
	cout << "uniform loc = " << uniformMVP << endl;
}

using namespace glm;

int bufferWidth, bufferHeight;
glm::mat4 View, Projection; // we need access to these two matrices inside the scroll_callback()


ostream& operator << (ostream& stream, const vec3& v)
{
	stream << v.x << ", " << v.y << ", " << v.z;
	return stream;
}

float clamp(const float v, const float vmin, const float vmax)
{
	if (v < vmin) return vmin;
	if (v > vmax) return vmax;
	return v;
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	static bool scrolling = false;
	static double prevx = -1, prevy = -1;
	//cout << xoffset << ", " << yoffset << endl;

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	//cout << "cursor pos = " << xpos << ", " << ypos << endl;

	if ((prevx != xpos) && (prevy != ypos))
	{
		prevx = xpos;
		prevy = ypos;
		cout << "NEW SCROLL STARTED\n";
	}

	// ZOOM PART
	if (1)
	{
		const auto view = View;
		const auto w = bufferWidth;
		const auto h = bufferHeight;

		vec3 pmin3 = glm::unProject(vec3(0, 0, 0), view, Projection, vec4(0, 0, w, h));
		vec3 pmax3 = glm::unProject(vec3(w, h, 0), view, Projection, vec4(0, 0, w, h));
		pmin3 = view * vec4(pmin3, 1);
		pmax3 = view * vec4(pmax3, 1);

		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		ypos = bufferHeight - ypos; // transform from GLFW-win-coords to OpenGL-win-coords
		xpos = clamp(xpos, 0, bufferWidth);
		ypos = clamp(ypos, 0, bufferHeight);

		vec3 zoomCenter3 = glm::unProject(vec3(xpos, ypos, 0), view, Projection, vec4(0, 0, bufferWidth, bufferHeight));
		zoomCenter3 = view * vec4(zoomCenter3, 1);


		vec2 zoomCenter(zoomCenter3.x, zoomCenter3.y);
		vec2 pmin(pmin3.x, pmin3.y);
		vec2 pmax(pmax3.x, pmax3.y);
		const float c = yoffset < 0 ? 1.01 : 0.99;
		auto dmin = c * (pmin - zoomCenter);
		auto dmax = c * (pmax - zoomCenter);
		pmin = zoomCenter + dmin;
		pmax = zoomCenter + dmax;

		//projection = glm::frustum(pmin.x,pmax.x, pmin.y, pmax.y, 1.f, 100.f);
		float near = -pmin3.z;
		float far = near + 100.f;
		Projection = glm::frustum(pmin.x, pmax.x, pmin.y, pmax.y, near, far);
	}

}


ostream& operator << (ostream& stream, const mat4& m4)
{
	auto* m = glm::value_ptr(m4);
	stream << m[0] << ", " << m[4] << ", " << m[8] << ", " << m[12] << endl;
	stream << m[1] << ", " << m[5] << ", " << m[9] << ", " << m[13] << endl;
	stream << m[2] << ", " << m[6] << ", " << m[10] << ", " << m[14] << endl;
	stream << m[3] << ", " << m[7] << ", " << m[11] << ", " << m[15] << endl;
	return stream;
}


int main()
{

	if (!glfwInit())
	{
		cout << "GLFW init failed!" << endl;
		glfwTerminate();
		return 1;
	}

	// set up glfw window properties
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	// core profile = no backwards compatinility
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	GLFWwindow* mainWindow = glfwCreateWindow(width, height, "Test Window", NULL, NULL);
	if (!mainWindow)
	{
		cout << "GLFW window creation failed!" << endl;
		glfwTerminate();
		return 1;
	}


	glfwSetScrollCallback(mainWindow, scroll_callback);


	// get buffer size information
	glfwGetFramebufferSize(mainWindow, &bufferWidth, &bufferHeight);

	// set context for GLFW to use
	glfwMakeContextCurrent(mainWindow);

	// allow modern extension features
	glewExperimental = GL_TRUE;

	if (glewInit() != GLEW_OK)
	{
		cout << "glew initialization failed!" << endl;
		glfwDestroyWindow(mainWindow);
		glfwTerminate();
		return 1;
	}


	glEnable(GL_DEPTH_TEST);

	glViewport(0, 0, bufferWidth, bufferHeight);

	createBox();
	createSphere(20, 10, 3);
	compileShaders();


	cout << "buffer width = " << bufferWidth << endl;
	cout << "buffer height = " << bufferHeight << endl;


	Projection = glm::perspective(45.f, (float)bufferWidth / bufferHeight, 1.f, 100.f);

	// loop until window is closed
	while (!glfwWindowShouldClose(mainWindow))
	{
		// get & handle user input events
		glfwPollEvents();

		if (direction)
			triOffset += triIncrement;
		else
			triOffset -= triIncrement;

		if (abs(triOffset) >= triMaxOffset)
			;// direction = !direction;


		// clear window
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		glUseProgram(shader);
		{
			// MODEL
			auto model = glm::mat4(1);
			model = glm::rotate(model, triOffset, glm::vec3(0, 1, 0));

			// CAMERA
			//View = glm::lookAt(vec3(0, 0, 4), vec3(0, 0, 0), vec3(0, 1, 0));
			View = glm::lookAtRH(vec3(0, 10, 10), vec3(0, 0, 0), vec3(0, 1, 0));



			auto mvp = Projection * View * model;
			glUniformMatrix4fv(uniformMVP, 1, GL_FALSE, glm::value_ptr(mvp));


			// DRAW CUBE
			glBindVertexArray(vao);
			{
				//glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
				//cout << "num indices = " << numIndices << endl;
				glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
			}
			glBindVertexArray(0);
		}
		glUseProgram(0);

		glfwSwapBuffers(mainWindow);
	}

	return 0;
}
