#pragma once
#include <vector>
class Grid
{
private:
	std::vector<float>* vec;
	int x,y;

public:
	Grid(int x , int y );
	float getItem(int i, int j);
	~Grid();
	
};

