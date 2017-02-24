#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
#include "orbital_tools.h"
#include "rendering_tools.h"

using namespace sf;

int main()
{
	double time_val = 0;
	vector<body*> bodies = orbital_init(time_val);

	vector<vec_n> coordinates = update_all_and_convert(bodies, time_val);

	vec_n coordinate = coordinates[1];

	coordinate.x;
	coordinate.y;


}