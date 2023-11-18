#include "Grid.h"

Grid::Grid(int x, int y)
{
	this->vec = new std::vector<float>(x * y);
	this->x = x;
	this->y = y;
}

float Grid::getItem(int i, int j)
{
	return this->vec->at(j * y + i);
}
