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

void main_render(std::vector<vec_n> cordinats, std::vector<CircleShape> planets, Sprite player) {

	Vector2f viewport_center;
	viewport_center.x = window2.getSize().x / 2;
	viewport_center.y = window2.getSize().y / 2;

	Vector2f modifi_cordinates;
	modifi_cordinates.x = cordinats[cordinats.size()-1].x-viewport_center.x;
	modifi_cordinates.y = cordinats[cordinats.size()-1].y - viewport_center.y;


	for (int i = 0; i < cordinats.size(); i++) {
		cordinats[i].x -= modifi_cordinates.x;
		cordinats[i].y -= modifi_cordinates.y;
	}
	player.setPosition(Vector2f(cordinats[cordinats.size()-1].x, cordinats[cordinats.size()-1].y));
	for (int i = 0; i < cordinats.size()-1; i++) {
		planets[i].setPosition(Vector2f(cordinats[i].x, cordinats[i].y));
	}

	window2.clear();
	
	window2.draw(player);

	for (int i = 0; i < planets.size();i++) {
		window2.draw(planets[i]);
	}

	window2.display();

}

void planetss() { //paste into begining of main function

	std::vector<CircleShape>planets(9);
	Sprite player;

	int radius[9] = {40, 14, 28, 15, 10, 16, 25, 23}; 

	bool texture_loaded[8]; //error log for texture loading
	for (int i = 0; i < 9; i++) {
		texture_loaded[i] = 1;
	}

	Texture planet_textures[9];
	if (!planet_textures[0].loadFromFile("sun_texture")) {
		texture_loaded[0] = 0;
	}
	if (!planet_textures[1].loadFromFile("planet_texture1")) {
		texture_loaded[0] = 0;
	}
	if (!planet_textures[2].loadFromFile("planet_texture2.png")) {
		texture_loaded[1] = 0;
	}
	if (!planet_textures[3].loadFromFile("planet_texture3.png")) {
		texture_loaded[2] = 0;
	}
	if (!planet_textures[4].loadFromFile("planet_texture4.png")) {
		texture_loaded[3] = 0;
	}
	if (!planet_textures[5].loadFromFile("planet_texture5.png")) {
		texture_loaded[4] = 0;
	}
	if (!planet_textures[6].loadFromFile("planet_texture6.png")) {
		texture_loaded[5] = 0;
	}
	if (!planet_textures[7].loadFromFile("planet_texture7.png")) {
		texture_loaded[6] = 0;
	}
	if (!planet_textures[8].loadFromFile("planet_texture8.png")) {
		texture_loaded[7] = 0;
	}
	if (!planet_textures[9].loadFromFile("planet_texture9.png")) {
		texture_loaded[8] = 0;
	}

	Texture player_texture;
	if (!player_texture.loadFromFile("player_texture.png")) {
		//handle exception
	}

	player.setTexture(player_texture);
	player.setOrigin(32, 32);


	for (int i = 0; i < planets.size(); i++) {
		planets[i].setRadius(radius[i]);
		planets[i].setTexture(&planet_textures[i]);
		planets[i].setOrigin(radius[i], radius[i]);
	}

}

bool collision(std::vector<CircleShape>planets, vec_n cordinats) {

	int radius[9];
	for (int i = 0; i < 9;i++){
		radius[i] = planets[i].getRadius();
	}
	
	float precision = 3;



}