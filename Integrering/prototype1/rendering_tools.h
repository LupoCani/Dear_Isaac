#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
#include "orbital_tools.h"

using namespace sf;

vec_n render() {

	Event input;

	double out_x = 0;
	double out_y = 0;

	while (window2.pollEvent(input)) {

		if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
			out_y++;
		}
		if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
			out_y--;
		}
		if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Right)) {
			out_x++;
		}
		if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Left)) {
			out_x--;
		}

	}
	vec_n out;
	out.x = out_x;
	out.y = out_y;

	return out;

}