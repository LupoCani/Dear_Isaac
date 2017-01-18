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

void main_render(vec_n cordinats[9], std::vector<CircleShape> planets, Sprite player) {

	Vector2f viewport_center;
	viewport_center.x = window2.getSize().x / 2;
	viewport_center.y = window2.getSize().y / 2;

	Vector2f modifi_cordinates;
	modifi_cordinates.x = cordinats[0].x-viewport_center.x;
	modifi_cordinates.y = cordinats[0].y - viewport_center.y;


	for (int i = 0; i < 9; i++) {
		cordinats[i].x -= modifi_cordinates.x;
		cordinats[i].y -= modifi_cordinates.y;
	}
	player.setPosition(Vector2f(cordinats[0].x, cordinats[0].y));
	for (int i = 0; i < 9; i++) {
		planets[i].setPosition(Vector2f(cordinats[i + 1].x, cordinats[i + 1].y));
	}

	window2.clear();
	
	window2.draw(player);

	for (int i = 0; i < 9;i++) {
		window2.draw(planets[i]);
	}

	window2.display();

}