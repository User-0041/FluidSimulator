#include "Simulation.h"
#include <algorithm>
Simulation::Simulation(float denstiy, int numX, int numY, float h)
{
	this->numX = numX + 2;
	this->numY = numY + 2;
	this->denstiy = denstiy;
	this->h = h;
	this->u = new std::vector<float>(this->numX * this->numY, 0.0f);
	this->v = new std::vector<float>(this->numX * this->numY, 0.0f);
	this->newU = new std::vector<float>(this->numX * this->numY, 0.0f);
	this->newV = new std::vector<float>(this->numX * this->numY, 0.0f);
	this->p = new std::vector<float>(this->numX * this->numY, 0.0f);
	this->s = new std::vector<float>(this->numX * this->numY,0.0f);
	this->m = new std::vector<float>(this->numX * this->numY, 1.0f);
	this->newM = new std::vector<float>(this->numX * this->numY);
}
void Simulation::initializeS() {
	for (int i = 0; i < numX; ++i) {
		for (int j = 0; j < numY; ++j) {
			if (i == 0 || i == numX - 1 || j == 0 || j == numY - 1) {
				// Set border values to 0
				s->at(i * numY + j) = 0.0f;
			}
			else {
				// Set inner values to 1
				s->at(i * numY + j) = 1.0f;
			}
		}
	}
}
int Simulation::getSize()
{
	return this->numX * this->numY;
}

void Simulation::integrate(float dt, float gravity)
{
	SPLICE;
	for (int i = 1; i < numX; i++)
	{
		for (int j = 1; j < numY-1; j++)
		{
			if (this->s->at(i * n + j) != 0.0 && this->s->at(i * n + j - 1) != 0.0) {
				this->v->at(i * n + j) +=(float) dt * gravity; 
			}
		}
	}
}


void Simulation::solveIncompressibility(int numIters, float dt)
{
	int n = this->numY;
	float cp = this->denstiy * this->h / dt;
	for (int itr = 0; itr < numIters; itr++) {
		for (int i = 1; i < numX -1; i++)
		{
			for (int j = 1; j < numY -1; j++)
			{
				if (this->s->at(i * n + j) == 0.0f) { continue; }
				float sx0 = this->s->at((i - 1) * n + j);
				float sx1 = this->s->at((i + 1) * n + j);
				float sy0 = this->s->at(i * n + j - 1);
				float sy1 = this->s->at(i * n + j + 1);
				float s = sx0 + sx1 + sy0 + sy1;
				if (s == 0.0f) { continue; };
				float div = this->u->at((i + 1) * n + j) - this->u->at(i * n + j) + this->v->at(i * n + j + 1) - this->v->at(i * n + j);
				
				float p = -div / s;
				p *= 1.9f;
				this->p->at(n * i + j) += cp * p;

				this->u->at(i * n + j) -= sx0 * p;
				this->u->at((i + 1) * n + j) += sx1 * p;
				this->v->at(i * n + j) -= sy0 * p;
				this->v->at(i * n + j + 1) += sy1 * p;
			}
		}
	}
}

void Simulation::extrapolate()
{
	SPLICE;
	for (int i = 0; i < numX; i++)
	{
		this->u->at(i * n + 0) = this->u->at(i * n + 1);
		this->u->at(i * n + this->numY-1) = this->u->at(i * n + this->numY - 2);

	}

	for (int j = 0; j < numY; j++)
	{
		this->v->at(0 * n + j) = this->v->at(1 * n + j);
		this->v->at( (this->numX-1) * n + j) = this->v->at((this->numX - 2) * n + j);
	}
}

float Simulation::sampleField(float x, float y, int FILED)
{

	SPLICE;
	float h = this->h;
	float h1 = 1.0f / h;
	float h2 = 0.5f / h1;
	using namespace std;
	float dx = 0.0f;
	float dy = 0.0f;
	std::vector<float>* f;
	switch (FILED)
	{
	case 0:
		f = this->u;
		dy = h2;
		break;
	case 1:
		f = this->v;
		dx = h2;
		break;
	case 2:
		dx = h2;
		dy = h2;
		f = this->m;
		break;
	default:
		f = nullptr;
		break;
	}
	//HACK:MAKE SURE CASTING DOSE NOT FUCK THIS UP
	int x0 =(int) min((float)floor((x - dx) * h1), (float)this->numX - 1);
	float tx = ((x - dx) - x0 * h) * h1;
	int x1 = (int)min(x0 + 1.0f, (float)this->numX - 1);

	int y0 = (int)min((float)floor((y - dy) * h1), (float)this->numY - 1);
	float ty = ((y - dy) - y0 * h) * h1;
	int y1 = (int)min(y0 + 1.0f, (float)this->numY - 1);

	float sx = 1.0f - tx;
	float sy = 1.0f - ty;

	float val = sx * sy * f->at(x0 * n + y0) +
				tx * sy * f->at(x1 * n + y0) +
				tx * ty * f->at(x1 * n + y1) +
				sx * ty * f->at(x0 * n + y1);
	return val;
}

float Simulation::avgV(int i, int j)
{
	SPLICE;
	float val = (
		this->v->at((i - 1) * n + j) +
		this->v->at(i * n + j) +
		this->v->at((i - 1) * n + j+1) +
		this->v->at(i * n + j + 1)) * 0.25f;
	return val;
}


float Simulation::avgU(int i, int j)
{

	SPLICE;
	float val =(
		this->u->at(i * n + j - 1) +
		this->u->at(i * n + j) +
		this->u->at((i + 1) * n + j - 1) +
		this->u->at((i + 1) * n + j)) * 0.25f;
	return val;
}




void Simulation::advectVel(float dt)
{
	copyVector(this->u, this->newU);
	copyVector(this->v, this->newV);
	SPLICE;
	float h = this->h;
	float h2 = 0.5f * h;
	for (int i = 1; i < numX - 1; i++) // Adjusted loop conditions to match JavaScript
	{
		for (int j = 1; j < numY - 1; j++) // Adjusted loop conditions to match JavaScript
		{
			if (this->s->at(i * n + j) != 0.0 && this->s->at((i - 1) * n + j) != 0.0 && j < this->numY - 1)
			{
				float x = i * h;
				float y = j * h + h2;
				float u = this->u->at(i * n + j);
				float v = this->avgV(i, j);
				x = x - dt * u;
				y = y - dt * v;
				u = sampleField(x, y, 0);
				this->newU->at(i * n + j) = u;
			}

			if (this->s->at(i * n + j) != 0.0 && this->s->at(i * n + j - 1) != 0.0 && i < this->numX - 1)
			{
				float x = j * h + h2;
				float y = j * h;
				float u = this->avgU(i, j);
				float v = this->v->at(i * n + j);
				x = x - dt * u;
				y = y - dt * v;
				v = sampleField(x, y, 1);
				this->newV->at(i * n + j) = v;
			}
		}
	}
	copyVector(this->newU, this->u);
	copyVector(this->newV, this->v);
}

void Simulation::advectSmoke(float dt)
{
	SPLICE;
	float h = this->h;
	float h2 = 0.5f * h;
	copyVector(this->m, this->newM);
	for (int i = 1; i < numX - 1; i++) // Adjusted loop conditions to match JavaScript
	{
		for (int j = 1; j < numY - 1; j++) // Adjusted loop conditions to match JavaScript
		{
			if (this->s->at(i * n + j) != 0.0f)
			{
				float u = (this->u->at(i * n + j) + this->u->at((i + 1) * n + j)) * 0.5f;
				float v = (this->v->at(i * n + j) + this->v->at(i * n + j + 1)) * 0.5f;
				float x = j * h + h2 - dt * u;
				float y = j * h + h2 - dt * v;

				this->newM->at(i * n + j) = sampleField(x, y, 2);
			}
		}
	}

	copyVector(this->newM, this->m);
}

void Simulation::copyVector(std::vector<float>* source, std::vector<float>* destination)
{
	std::memcpy(destination->data(), source->data(), source->size() * sizeof(float));
}

