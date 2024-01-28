#include <iostream>
#include <stdio.h>
#include <vector>

# define M_PI          3.1415L /* pi */
using namespace std;

#define __attribute__(x) //  for the keyboard



// constexpr runs at compile time and  is generally faster than const.

std::vector<float> G = { 0.f, -10.f };  // external (gravitational) forces



constexpr float REST_DENS = 1000.f;   // rest density
constexpr float GAS_CONST = 200.f;   // const for equation of state
constexpr float H = 16.f;			   // kernel radius( radius of a particle). if the H decreases so does the size of the particle
constexpr float HSQ = H * H;		   // radius^2 for optimization
constexpr float MASS = 2.5f;		   // assume all particles have the same mass
constexpr float VISCOSITY = 200.f;   // viscosity constant
constexpr float dt = 0.0007f;	       // integration timestep. should be adjusted according to parameters

// smoothing kernels defined in 
const float POLY6 = 315.f / (64.f * (M_PI * pow(H, 8.f)));
const float SPIKY_GRAD = -10.f / (M_PI * pow(H, 5.f));
const float VISC_LAP = 40.f / (M_PI * pow(H, 5.f));

// simulation parameters
constexpr float EPS = H; // boundary epsilon
float BOUND_DAMPING = -0.5f;
// interaction
constexpr int MAX_PARTICLES = 25000;
constexpr int DAM_PARTICLES = 500;
constexpr int BLOCK_PARTICLES = 250;

// rendering projection parameters
constexpr double VIEW_WIDTH = 1.5 * 800.f;
constexpr double VIEW_HEIGHT = 1.5 * 600.f;



// SPH VALUES
// 
// particle data structure
// stores position, velocity, and force for integration
// stores density (rho) 
// pressure values  

class SINGLE_Particle {
public:
	std::vector<float> x, v, f;  // position, velocity, force
	float rho, p;

	// Constructor
	SINGLE_Particle(float x_pos, float y_pos) : rho(0), p(0.f) {
		x = { x_pos, y_pos };
		v = { 0.f, 0.f };
		f = { 0.f, 0.f };
	}
};


// solver data
vector<SINGLE_Particle> particles;

//thrust::device_vector<SINGLE_Particle> particles(numParticles);

void Integrate(void) {
	for (auto& p : particles) {
		// Forward Euler integration (equation 8 from the provided link)

		// Update velocity
		for (size_t i = 0; i < p.v.size(); ++i) {
			p.v[i] += dt * p.f[i] / p.rho;
			p.x[i] += dt * p.v[i];
		}

		// Enforce boundary conditions
		for (size_t i = 0; i < p.x.size(); ++i) {
			if (p.x[i] - EPS < 0.f) {
				p.v[i] *= BOUND_DAMPING;
				p.x[i] = EPS;
			}
			if (p.x[i] + EPS > (i == 0 ? VIEW_WIDTH : VIEW_HEIGHT)) {
				p.v[i] *= BOUND_DAMPING;
				p.x[i] = (i == 0 ? VIEW_WIDTH : VIEW_HEIGHT) - EPS;
			}
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
			std::vector<float> rij(2);
			rij[0] = pj.x[0] - pi.x[0];
			rij[1] = pj.x[1] - pi.x[1];

			float r2 = rij[0] * rij[0] + rij[1] * rij[1];

			if (r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS * POLY6 * std::pow(HSQ - r2, 3.f);
			}
		}
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}

void ComputeForces(void)
{
	for (auto& p_i : particles)
	{
		std::vector<float> F_pressure(2, 0.f);
		std::vector<float> F_viscosity(2, 0.f);

		for (auto& p_j : particles)
		{
			if (&p_i == &p_j)
			{
				continue;
			}

			std::vector<float> rij(2);
			rij[0] = p_j.x[0] - p_i.x[0];
			rij[1] = p_j.x[1] - p_i.x[1];

			float r = std::sqrt(rij[0] * rij[0] + rij[1] * rij[1]);

			if (r < H)
			{
				// compute pressure force contribution
				F_pressure[0] += -rij[0] / r * MASS * (p_i.p + p_j.p) / (2.f * p_j.rho) * SPIKY_GRAD * std::pow(H - r, 3.f);
				F_pressure[1] += -rij[1] / r * MASS * (p_i.p + p_j.p) / (2.f * p_j.rho) * SPIKY_GRAD * std::pow(H - r, 3.f);

				// compute viscosity force contribution
				F_viscosity[0] += VISCOSITY * MASS * (p_j.v[0] - p_i.v[0]) / p_j.rho * VISC_LAP * (H - r);
				F_viscosity[1] += VISCOSITY * MASS * (p_j.v[1] - p_i.v[1]) / p_j.rho * VISC_LAP * (H - r);
			}
		}

		std::vector<float> F_gravity = { 0.f, G[1] * (MASS / p_i.rho) };
		p_i.f[0] = F_pressure[0] + F_viscosity[0] + F_gravity[0];
		p_i.f[1] = F_pressure[1] + F_viscosity[1] + F_gravity[1];
	}
}
void InitSPH(void) {

	std::cout << "initialzing dam break" << DAM_PARTICLES << "particles" << std::endl;
	for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H)
		for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H)
			if (particles.size() < DAM_PARTICLES) {

				float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				particles.push_back(SINGLE_Particle(x + jitter, y));
			}
}


