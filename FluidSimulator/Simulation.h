#pragma once
#include <vector>
#define SPLICE int n = this->numY
class Simulation
{
public:
	float denstiy;
	int numX;
	int numY;
	float h;
	std::vector<float>* u;
	std::vector<float>* v;
	std::vector<float>* newV;
	std::vector<float>* newU;
	std::vector<float>* p;
	std::vector<float>* s;
	std::vector<float>* m;
	std::vector<float>* newM;
	Simulation(float denstiy, int numX, int numY  , float h);
	int getSize();
	void integrate(float dt, float gravity);
	void solveIncompressibility(int numIters, float dt);
	void extrapolate();
	float sampleField(float x, float y,int FILED);
	float avgU(int i, int j);
	float avgV(int i, int j);
	void advectVel(float dt);
	void initializeS();
	void advectSmoke(float dt);
	void copyVector(std::vector<float>* source, std::vector<float>* destination);
};