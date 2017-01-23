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

	Vector2f viewport_center; //cordinates at the center of the viewport/window
	viewport_center.x = window2.getSize().x / 2; 
	viewport_center.y = window2.getSize().y / 2;

	Vector2f modifi_cordinates; //value to modifi curent cordinates with on oder to put player att the center of the 
	modifi_cordinates.x = cordinats[cordinats.size()-1].x-viewport_center.x;
	modifi_cordinates.y = cordinats[cordinats.size()-1].y - viewport_center.y;


	for (int i = 0; i < cordinats.size(); i++) {
		cordinats[i].x -= modifi_cordinates.x;
		cordinats[i].y -= modifi_cordinates.y;
	}
	//set new position based on calculated values 
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

	std::vector<CircleShape>planets(9); //pass to render and  collision function
	Sprite player; //pass to render and  collision function

	int radius[9] = {40, 14, 28, 15, 10, 16, 25, 23}; //radius for the planets 

	bool texture_loaded[8]; //error log for texture loading (1 texture was loaded, 0 texture was not loaded)
	for (int i = 0; i < 9; i++) {
		texture_loaded[i] = 1;
	}

	Texture planet_textures[9];
	if (!planet_textures[0].loadFromFile("sun_texture")) {
		texture_loaded[0] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[1].loadFromFile("planet_texture1")) {
		texture_loaded[0] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[2].loadFromFile("planet_texture2.png")) {
		texture_loaded[1] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[3].loadFromFile("planet_texture3.png")) {
		texture_loaded[2] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[4].loadFromFile("planet_texture4.png")) {
		texture_loaded[3] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[5].loadFromFile("planet_texture5.png")) {
		texture_loaded[4] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[6].loadFromFile("planet_texture6.png")) {
		texture_loaded[5] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[7].loadFromFile("planet_texture7.png")) {
		texture_loaded[6] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[8].loadFromFile("planet_texture8.png")) {
		texture_loaded[7] = 0; //log texture not loaded in error log
	}
	if (!planet_textures[9].loadFromFile("planet_texture9.png")) {
		texture_loaded[8] = 0; //log texture not loaded in error log
	}

	Texture player_texture;
	if (!player_texture.loadFromFile("player_texture.png")) {
		//handle exception
	}

	player.setTexture(player_texture);
	player.setOrigin(32, 32); //center the origin of the player (half the with, half the height)
	float player_radius = player.getLocalBounds().width / 2; // radius of circle containing sprite; pass to collision function


	for (int i = 0; i < planets.size(); i++) { //set planet values
		planets[i].setRadius(radius[i]);
		planets[i].setTexture(&planet_textures[i]);
		planets[i].setOrigin(radius[i], radius[i]);
	}

}

bool collision(std::vector<CircleShape>planets, std::vector<vec_n> cordinats,Sprite player, float player_radius) { //function to detect if objects have colided 

	float player_x = cordinats[cordinats.size() - 1].x; //get player x position from last entry in cordinate vector
	float player_y = cordinats[cordinats.size() - 1].y; //get player y position from last entry in cordinate vector
	bool collision = 0; //0==no collision detected, 1==collision detected 
	
	float precision = 3;
	for (float i = cordinats[cordinats.size() - 1].x - player.getLocalBounds().width / 2; i < cordinats[cordinats.size() - 1].x + player.getLocalBounds().width / 2; i += 1) {

		float player_rout_collision = sqrt(player_radius*player_radius - (i - player_x)*(i - player_x)); //root of the circle equation for the circle containing the player (sqr(r^2-(x-x0)^2)
		float player_pos_collision = player_rout_collision + player_y; //value given by the positiv root
		float player_neg_collision = -player_rout_collision + player_y; //value given by the negativ root

		for (float k = 0; k < planets.size(); k++) { //circle trough the planets vector

			float radius = planets[k].getRadius();  
			float pos_x = cordinats[k].x; //get planet x cordinate from cordinats vector
			float pos_y = cordinats[k].y;  //get planet y cordinate from cordinats vector

			if (radius*radius - (i - pos_x)*(i - pos_x) >= 0) {

				float rout_collision = sqrt(radius*radius - (i - pos_x)*(i - pos_x)); //root of the circle equation for the Circleshape planets[k] (sqr(r^2-(x-x0)^2)
				float pos_collision = rout_collision + pos_y; //value given by the positiv root
				float neg_collision = -rout_collision + pos_y; //value given by the negativ root

				if (pos_collision >= player_pos_collision - precision && pos_collision <= player_pos_collision + precision) {
					collision = 1;
					break;
				}
				if (pos_collision >= player_neg_collision - precision && pos_collision <= player_neg_collision + precision) {
					collision = 1;
					break;
				}
				if (neg_collision >= player_pos_collision - precision && neg_collision <= player_pos_collision + precision) {
					collision = 1;
					break;
				}
				if (neg_collision >= player_neg_collision - precision && neg_collision <= player_neg_collision + precision) {
					collision = 1;
					break;
				}

			}

		}

		if (collision == 1) {
			break;
		}

	}

	return (collision);
}