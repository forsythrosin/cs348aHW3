#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GL/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<CurvatureInfo> curvature;
Mesh mesh;

bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = {0,0,4,1};
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
bool showSurface = true, showAxes = true, showCurvature = false, showNormals = false;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float ambient[] = { 0.1, 0.1, 0.1, 1.0 };
float diffuse[] = {1.0,0.8,0.6,1.0};
float shininess[] = { 50.0 };

void renderMesh() {
	if (!showSurface) glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); // render regardless to remove hidden lines
	
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);

	glDepthRange(0.001,1);
	glEnable(GL_NORMALIZE);

	// WRITE CODE HERE TO RENDER THE TRIANGLES OF THE MESH
	// ---------------------------------------------------------
	OpenMesh::Vec3f point, normal;
	Mesh::FaceIter f_it, f_end(mesh.faces_end());
	glBegin(GL_TRIANGLES);
	for (f_it = mesh.faces_begin(); f_it != f_end; ++f_it){
	  Mesh::FaceHandle fh = f_it.handle();
	  Mesh::FaceVertexIter fv_it;
	  for (fv_it = mesh.fv_iter(fh); fv_it; ++fv_it) {
	    Mesh::VertexHandle vh = fv_it.handle();
	    normal = mesh.normal(vh);
	    glNormal3f(normal.values_[0], normal.values_[1], normal.values_[2]);
	    point = mesh.point(vh);
	    glVertex3f(point.values_[0], point.values_[1], point.values_[2]);
	  }
	}
	glEnd();
	// -------------------------------------------------------------------------------------------------------------
	
	if (!showSurface) glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	
	glDisable(GL_LIGHTING);
	glDepthRange(0,0.999);
	
	if (showCurvature) {
		// WRITE CODE HERE TO RENDER THE PRINCIPAL DIRECTIONS YOU COMPUTED ---------------------------------------------
	        glBegin(GL_LINES);
		for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it
		       != mesh.vertices_end(); ++it){
		       CurvatureInfo info = mesh.property(curvature, it);
		       Vec3f p = mesh.point(it.handle()),
			 d1 = info.directions[0].normalized(),
			 d2 = info.directions[1].normalized();
		       float k1 = info.curvatures[0],
			 k2 = info.curvatures[1];
		       float vecLength = 0.02f;
		       Vec3f p10 = p - vecLength/2 * d1,
			 p11 = p + vecLength/2 * d1,
			 p20 = p - vecLength/2 * d2,
			 p21 = p + vecLength/2 * d2;
                       /*Vec3f p10 = p - vecLength/2 * d1,
			 p11 = p + vecLength/2 * d1,
			 p20 = p - (k2/k1)*vecLength/2 * d2,
			 p21 = p + (k2/k1)*vecLength/2 * d2;*/

		       glColor3f(1,0,0); // maximum curvature direction
		       glVertex3f(p10[0],p10[1],p10[2]);
		       glVertex3f(p11[0],p11[1],p11[2]);
		       glColor3f(0,0,1); // minimum curvature direction
		       glVertex3f(p20[0],p20[1],p20[2]);
		       glVertex3f(p21[0],p21[1],p21[2]);
		}
		glEnd();
		// -------------------------------------------------------------------------------------------------------------
	}
	
	if (showNormals) {
		glBegin(GL_LINES);
		glColor3f(0,1,0);
		for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
			Vec3f n = mesh.normal(it.handle());
			Vec3f p = mesh.point(it.handle());
			Vec3f d = p + n*.01;
			glVertex3f(p[0],p[1],p[2]);
			glVertex3f(d[0],d[1],d[2]);
		}
		glEnd();
	}
	
	glDepthRange(0,1);
}

void display() {
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

        glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
	glEnable(GL_LIGHT0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,windowWidth,windowHeight);
	
	float ratio = (float)windowWidth / (float)windowHeight;
	gluPerspective(50, ratio, 0.1, 1000); // 50 degree vertical viewing angle, zNear = 0.1, zFar = 1000
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPos[0]+pan[0], cameraPos[1]+pan[1], cameraPos[2]+pan[2], pan[0], pan[1], pan[2], up[0], up[1], up[2]);
	
	// Draw mesh
	renderMesh();

	// Draw axes
	if (showAxes) {
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glLineWidth(1);
			glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(1,0,0); // x axis
			glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,1,0); // y axis
			glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,1); // z axis
		glEnd(/*GL_LINES*/);
	}
	
	glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
	else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
	else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);
	
	lastPos[0] = x;
	lastPos[1] = y;
}

void mouseMoved(int x, int y) {
	int dx = x - lastPos[0];
	int dy = y - lastPos[1];
	Vec3f curCamera(cameraPos[0],cameraPos[1],cameraPos[2]);
	Vec3f curCameraNormalized = curCamera.normalized();
	Vec3f right = up % curCameraNormalized;

	if (leftDown) {
		// Assume here that up vector is (0,1,0)
		Vec3f newPos = curCamera - 30*(float)((float)dx/(float)windowWidth) * right + 30*(float)((float)dy/(float)windowHeight) * up;
		newPos = newPos.normalized() * curCamera.length();
		
		up = up - (up | newPos) * newPos / newPos.sqrnorm();
		up.normalize();
		
		for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
	}
	else if (rightDown) for (int i = 0; i < 3; i++) cameraPos[i] *= pow(1.1,dy*.1);
	else if (middleDown) {
		pan += -30*(float)((float)dx/(float)windowWidth) * right + 30*(float)((float)dy/(float)windowHeight) * up;
	}

	
	lastPos[0] = x;
	lastPos[1] = y;
	
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);

	if (key == 's' || key == 'S') showSurface = !showSurface;
	else if (key == 'a' || key == 'A') showAxes = !showAxes;
	else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
	else if (key == 'n' || key == 'N') showNormals = !showNormals;
	else if (key == 'q' || key == 'Q') exit(0);
	glutPostRedisplay();
}

void reshape(int width, int height) {
	windowWidth = width;
	windowHeight = height;
	glutPostRedisplay();
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " mesh_filename\n";
		exit(0);
	}
	
	IO::Options opt;
	opt += IO::Options::VertexNormal;
	opt += IO::Options::FaceNormal;
	
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	
	cout << "Reading from file " << argv[1] << "...\n";
	if ( !IO::read_mesh(mesh, argv[1], opt )) {
		cout << "Read failed.\n";
		exit(0);
	}
	
	mesh.update_normals();
	mesh.add_property(curvature);
	
	// Move center of mass to origin
	Vec3f center(0,0,0);
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) center += mesh.point(vIt);
	center /= mesh.n_vertices();
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) -= center;

	// Fit in the unit sphere
	float maxLength = 0;
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) maxLength = max(maxLength, mesh.point(vIt).length());
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) /= maxLength;
	
	computeCurvature(mesh,curvature);

	cout << "Mesh stats:\n";
	cout << '\t' << mesh.n_vertices() << " vertices.\n";
	cout << '\t' << mesh.n_edges() << " edges.\n";
	cout << '\t' << mesh.n_faces() << " faces.\n";
	
	up = Vec3f(0,1,0);
	pan = Vec3f(0,0,0);

	glutInit(&argc, argv); 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 
	glutInitWindowSize(windowWidth, windowHeight); 
	glutCreateWindow(argv[0]);

	glutDisplayFunc(display);
	glutMotionFunc(mouseMoved);
	glutMouseFunc(mouse);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);

	glutMainLoop();
	
	return 0;
}
