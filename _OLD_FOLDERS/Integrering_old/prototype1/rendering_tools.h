#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
#include "orbital_tools.h"

using namespace sf;

sf::Font font;
Color yellow(225, 237, 7);
Color Red(165, 41, 13);
Color grey(177, 190, 198);
Color white(255, 255, 255);

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

	if (!font.loadFromFile("Alger.ttf")){
		//handle exception
	}

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

bool collision(std::vector<CircleShape>planets, std::vector<vec_n> cordinats, Sprite player, float player_radius) {

	Vector2f player_pos;
	player_pos.x = cordinats[0].x;
	player_pos.y = cordinats[0].y;
	bool collided = 0;

	for (int i = 1; i < cordinats.size; i++) {
		Vector2f planet_pos;
		planet_pos.x = cordinats[i].x;
		planet_pos.y = cordinats[i].y;

		float radius_compare = planets[i-1].getRadius() + player_radius;

		Vector2f diference;
		diference.x = player_pos.x - planet_pos.x;
		diference.y = player_pos.y - planet_pos.y;

		if (sqrt(diference.x*diference.x + diference.y*diference.y) <= radius_compare)
		{
			collided = 1;
			break;
		}
	}

	return (collided);
}

void option_menue() {

	int header_pos_y[2] = { 400,500 };
	Text option_header[2];
	option_header[0].setString("Sound");
	option_header[1].setString("BACK");
	for (int i = 0; i <= 1; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(40);
		option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
	bool run_option = 1;
	bool selected[2] = { 0,0 };

	while (run_option == 1) {
		Event input;

		for (int i = 0; i <= 1; i++) {
			float activ_header_poss = option_header[i].getPosition().y + option_header[i].getLocalBounds().height*1.3;
			if (underline.getPosition().y == activ_header_poss) {
				selected[i] = 1;
			}
			if (underline.getPosition().y != activ_header_poss) {
				selected[i] = 0;
			}
		}

		while (window2.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 1; i++) {
					if (selected[i] == 1) {
						if (i == 0) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
				for (int i = 0; i <= 1; i++) {
					if (selected[i] == 1) {
						if (i == 1) {
							underline.setSize(Vector2f(option_header[i - 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[1] + option_header[1].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code != Keyboard::Up) && (input.key.code != Keyboard::Down)) {
				if (selected[1] == 1) {
					run_option = 0;
					break;
				}
			}

		}

		window2.clear();
		for (int i = 0; i <= 1; i++) {
			window2.draw(option_header[i]);
		}
		window2.draw(underline);
		window2.display();

	}

}

void ingame_menue() {



	int header_pos_y[3] = { 340,440,550 };
	Text option_header[3];
	option_header[0].setString("Sound");
	option_header[1].setString("Quit to Main menue");
	option_header[2].setString("Back");
	for (int i = 0; i <= 2; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(50);
		option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
	bool run_option = 1;
	bool selected[3] = { 0,0,0 };

	int output;

	while (run_option == 1) {
		Event input;

		for (int i = 0; i <= 2; i++) {
			float activ_header_poss = option_header[i].getPosition().y + option_header[i].getLocalBounds().height*1.3;
			if (underline.getPosition().y == activ_header_poss) {
				selected[i] = 1;
			}
			if (underline.getPosition().y != activ_header_poss) {
				selected[i] = 0;
			}
		}

		while (window2.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i <= 1) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i >= 1) {
							underline.setSize(Vector2f(option_header[i - 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[2].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[2] + option_header[2].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code != Keyboard::Up) && (input.key.code != Keyboard::Down)) {
				if (selected[2] == 1) {
					run_option = 0;

					break;
				}
				if (selected[1] == 1) {
					run_option = 0;
					
					break;
				}
			}

		}

		window2.clear();
		for (int i = 0; i <= 2; i++) {
			window2.draw(option_header[i]);
		}
		window2.draw(underline);
		window2.display();

	}

	}

int start_menue() {
	Text title;
	title.setFont(font);
	title.setString("Dear Isaac");
	title.setCharacterSize(60);
	title.setFillColor(white);
	title.setPosition(Vector2f(window2.getSize().x / 2 - title.getLocalBounds().width*0.5, 200));

	int header_pos_y[3] = { 400,500,600 };

	

	Text option_header[3];
	option_header[0].setString("PLAY");
	option_header[1].setString("OPTION");
	option_header[2].setString("QUIT");
	for (int i = 0; i <= 2; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(40);
		option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));

	bool play = 0;
	bool selected[3] = { 0,0,0 };

	while (play == 0) {
		Event input;

		for (int i = 0; i <= 2; i++) {
			float activ_header_poss = option_header[i].getPosition().y + option_header[i].getLocalBounds().height*1.3;
			if (underline.getPosition().y == activ_header_poss) {
				selected[i] = 1;
			}
			if (underline.getPosition().y != activ_header_poss) {
				selected[i] = 0;

			}
		}

		while (window2.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i <= 1) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i >= 1) {
							underline.setSize(Vector2f(option_header[i - 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[2].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[2] + option_header[2].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code != Keyboard::Up) && (input.key.code != Keyboard::Down)) {
				if (selected[0] == 1) {
					play = 1;
				}
				if (selected[2] == 1) {
					window2.close();
					return 0;
					break;
				}
				if (selected[1] == 1) {
					option_menue();
				}
			}

		}
		window2.clear();
		for (int i = 0; i <= 2; i++) {
			window2.draw(option_header[i]);
		}
		window2.draw(title);
		window2.draw(underline);
		window2.display();


	}