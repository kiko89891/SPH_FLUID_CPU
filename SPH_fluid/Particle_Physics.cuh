


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <vector>
#include <iostream>
# define M_PI          3.141592653589793238462643383279502884L /* pi */
using namespace std;

#define __attribute__(x) //  for the keyboard


#include <Eigen/Dense>
using namespace Eigen;
__global__ void cuda_hello() {
	printf("Hello World from GPU!\n");
}


// constexpr runs at compile time and  is generally faster than const.

Vector2d G(0.f, -9.81f);   // external (gravitational) forces
constexpr float REST_DENS = 300.f;   // rest density
constexpr float GAS_CONST = 300.f;   // const for equation of state
constexpr float H = 5.f;			   // kernel radius( radius of a particle). if the H decreases so does the size of the particle
constexpr float HSQ = H * H;		   // radius^2 for optimization
constexpr float MASS = 20.f;		   // assume all particles have the same mass
constexpr float VISCOSITY = 100.f;   // viscosity constant
constexpr float dt = 0.001f;	       // integration timestep. should be adjusted according to parameters

// smoothing kernels defined in 
const float POLY6 = 315.f / (64.f * (M_PI * pow(H, 9.f)));
const float SPIKY_GRAD = -15.f / (M_PI * pow(H, 5.f));
const float VISC_LAP = 40.f / (M_PI * pow(H, 5.f));

// simulation parameters
constexpr float EPS = H; // boundary epsilon
float BOUND_DAMPING = -0.5f;
// interaction
constexpr int MAX_PARTICLES = 10000;
constexpr int DAM_PARTICLES = 1000;
constexpr int BLOCK_PARTICLES = 250;

// rendering projection parameters
constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 600;
constexpr double VIEW_WIDTH = 1.5 * 800.f;
constexpr double VIEW_HEIGHT = 1.5 * 600.f;



// SPH VALUES
// 
// particle data structure
// stores position, velocity, and force for integration
// stores density (rho) 
// pressure values  

struct SIGNLE_Particle
{
	SIGNLE_Particle(float _x, float _y) : x(_x, _y), v(0.f, 0.f), f(0.f, 0.f), rho(0), p(0.f) {}
	Vector2d x, v, f; //position   ,  velocity   , force
	//vector<vector<float>> x,y,z;
	float rho, p;

};
// solver data
vector<SIGNLE_Particle> particles;




void Integrate(void)
{
	for (auto& p : particles)
	{

		// https://matthias-research.github.io/pages/publications/sca03.pdf equation(8) on page 3 
		// forward Euler integration

		p.v += dt * p.f / p.rho;
		p.x += dt * p.v;

		// enforce boundary conditions
		if (p.x(0) - EPS < 0.f)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = EPS;
		}
		if (p.x(0) + EPS > VIEW_WIDTH)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = VIEW_WIDTH - EPS;
		}
		if (p.x(1) - EPS < 0.f)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = EPS;
		}
		if (p.x(1) + EPS > VIEW_HEIGHT)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = VIEW_HEIGHT - EPS;
		}
	}
}

void DensityPressure(void)
{
	for (auto& pi : particles)
	{
		pi.rho = 0.f;
		for (auto& pj : particles)
		{
			Vector2d rij = pj.x - pi.x;
			float r2 = rij.squaredNorm();

			if (r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
			}
		}
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}

void ComputeForces(void)
{
	for (auto& p_i : particles)

	{   // initalize a new vector (x,y) for each particle up to all particles

		Vector2d F_pressure(0.f, 0.f);
		Vector2d F_viscosity(0.f, 0.f);

		// 
		for (auto& p_j : particles)
		{
			if (&p_i == &p_j)
			{
				continue;
			}
			//  r - square root of the sum of the absolute squares of its elements( Frobenius Norm )
			Vector2d rij = p_j.x - p_i.x;
			float r = rij.norm();

			if (r < H)
			{
				// compute pressure force contribution
				F_pressure += -rij.normalized() * MASS * (p_i.p + p_j.p) / (2.f * p_j.rho) * SPIKY_GRAD * pow(H - r, 3.f);
				// compute viscosity force contribution
				F_viscosity += VISCOSITY * MASS * (p_j.v - p_i.v) / p_j.rho * VISC_LAP * (H - r);
			}
			else {
				F_viscosity(0.f, 0.f);
			}
		}
		Vector2d F_gravity = G * MASS / p_i.rho;
		p_i.f = F_pressure + F_viscosity + F_gravity;
	}
}
