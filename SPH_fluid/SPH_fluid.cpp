//// "Smooth particle-based fluid simulation" -SPH BY Kiril Hristov
//https://www.youtube.com/watch?v=HzFatL3WT6g

#include <stdio.h>
#include <GL/glut.h>
#include <iostream>
#include <vector>

#include "Particle_Physics.cuh"


#include <Eigen/Dense>




void InitSPH(void) {

	std::cout << "initialzing dam break" << DAM_PARTICLES << "particles" << std::endl;
	for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H)
		for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H)
			if (particles.size() < DAM_PARTICLES) {

				float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				particles.push_back(SIGNLE_Particle(x + jitter, y));
			}
}


void Update(void)
{
	DensityPressure();
	ComputeForces();
	Integrate();

	glutPostRedisplay();
}

void InitGL(void)
{
	glClearColor(0.9f, 0.9f, 0.9f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(H / 2.f);
	glMatrixMode(GL_PROJECTION);
}

void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

	glColor4f(0.0f, 0.0f, 3.f, 2);
	glBegin(GL_POINTS);
	for (auto& p : particles)
	{
		glVertex2f(p.x(0), p.x(1));
	}
	glEnd();

	glutSwapBuffers();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{
	switch (c)
	{
	case 'r':
	case 'R':
		particles.clear();
		InitSPH();
		break;
	}
}

// https://www.youtube.com/watch?v=idGYpwvxBTo
// https://www.youtube.com/watch?v=GCCplaZSZBE&t=310s
void mouse(int button, int state, int mousex, int mousey)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{

		puts("Left button clicked");
		if (particles.size() >= MAX_PARTICLES)
		{
			std::cout << "maximum number of particles reached" << std::endl;
		}
		else
		{
			//VIEW_WIDTH
			// VIEW_HEIGHT
			float x, y;
			unsigned int placed = 0;
			x = mousex+150;

			y = VIEW_HEIGHT -mousey-150;

			
			placed++;
			while (placed++ < BLOCK_PARTICLES && particles.size() < MAX_PARTICLES)
			{             // ideally to the kernel raidus H . It will burst a certain number of particles when rightclick

				particles.push_back(SIGNLE_Particle(x+ rand() %(16 + 1 ), y + rand() %(16 + 1)));
			}
		}
	}
}



int main(int argc, char** argv)
{
	cuda_hello();
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInit(&argc, argv);
	glutCreateWindow("KIRIL SPH");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(mouse);
	InitGL();
	InitSPH();

	glutMainLoop();
	return 0;
}

// Bibliography 

// https://www.cs.cornell.edu/%7Ebindel/class/cs5220-f11/code/sph-derive.pdf
// https://lucasschuermann.com/writing/particle-based-fluid-simulation
// https://math.stackexchange.com/questions/3898821/laplacian-of-viscosity-kernel-function-m%C3%BCller-sph
// 
// 
// 
// setup cuda
// https://medium.com/@aviatorx/c-and-cuda-project-visual-studio-d07c6ad771e3