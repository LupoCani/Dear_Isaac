#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
#define setFillColor2 setFillColor

namespace graph{

	using namespace phys;
	using namespace sf;

	sf::Font font;
	Color yellow(225, 237, 7);
	Color Red(165, 41, 13);
	Color grey(177, 190, 198);
	Color white(255, 255, 255);

	struct init_out
	{

	};

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

	vector<vec_n> handle_scale(vector<vec_n> list, int focus, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		vec_n origo = list[focus];
		origo.x *= scale;
		origo.y *= scale;

		for (int i = 0; i < list.size(); i++)
		{
			list[i].x *= scale;
			list[i].y *= scale;

			list[i].x = 500 + (list[i].x - origo.x);
			list[i].y = 500 - (list[i].y - origo.y);
		}

		return list;
	}

	std::vector<CircleShape>planets;
	std::vector<Texture> textures;
	Sprite player;
	Texture player_s;
	double player_radius;

	void main_render(std::vector<vec_n> coordinates, double zoom, int last_i = 0) {


		Vector2f viewport_center; //cordinates at the center of the viewport/window
		viewport_center.x = window2.getSize().x / 2;
		viewport_center.y = window2.getSize().y / 2;
		/*
		Vector2f modify_cordinates; //value to modify curent cordinates with on oder to put player att the center of the
		modify_cordinates.x = coordinates[coordinates.size() - 1].x - viewport_center.x;
		modify_cordinates.y = coordinates[coordinates.size() - 1].y - viewport_center.y;
		//modify_cordinates.x = coordinates[0].x;// - viewport_center.x;
		//modify_cordinates.y = coordinates[0].y;// - viewport_center.y;

		for (int i = 0; i < coordinates.size(); i++) {
		coordinates[i].x -= modify_cordinates.x;
		coordinates[i].y -= modify_cordinates.y;
		}
		//*/

		//set new position based on calculated values 
		vector<double> scales(0);

		coordinates = handle_scale(coordinates, last_i, zoom, viewport_center.x, viewport_center.y);

		player.setScale(0.5 *zoom, 0.5 *zoom);
		player.setPosition(coordinates[last_i]);
		for (int i = 0; i < 9; i++) {
			planets[i].setPosition(coordinates[i]);

			scales.push_back(planets[i].getRadius());
			planets[i].setRadius(scales[i] * zoom);
			planets[i].setOrigin(planets[i].getRadius() / 2, planets[i].getRadius() / 2);
		}

		window2.clear();

		for (int i = 0; i < planets.size(); i++) {
			window2.draw(planets[i]);
		}

		window2.draw(player);
		window2.display();

		for (int i = 0; i < 9; i++) {
			planets[i].setRadius(scales[i]);
		}
	}


	//Texture planet_textures[9];
	//*


	void render_init() { //paste into begining of main function


		if (!font.loadFromFile("ALGER.ttf")) {
			//handle exception
		}

		std::vector<CircleShape>planets_out(9); //pass to render and  collision function
		planets = planets_out;
		Sprite player; //pass to render and  collision function

		//int radius[9] = { 220, 68, 116, 70, 60, 72, 110, 86, 150 }; //radius for the planets 
		int radius[9] = { 20, 20, 20, 20, 20, 20, 20, 20, 20 }; //radius for the planets 

		vector<Texture> tx_out(9);
		textures = tx_out;

		if (!textures[0].loadFromFile("sun_texture.png")) {

		}
		if (!textures[1].loadFromFile("planet_texture2.png")) {

		}
		if (!textures[2].loadFromFile("planet_texture3.png")) {

		}
		if (!textures[3].loadFromFile("planet_texture4.png")) {

		}
		if (!textures[4].loadFromFile("planet_texture5.png")) {

		}
		if (!textures[5].loadFromFile("planet_texture6.png")) {

		}
		if (!textures[6].loadFromFile("planet_texture7.png")) {

		}
		if (!textures[7].loadFromFile("planet_texture8.png")) {

		}
		if (!textures[8].loadFromFile("planet_texture9.png")) {

		}

		Texture player_texture;
		if (!player_s.loadFromFile("Character_sprite.png")) {
			//handle exception
		}

		player.setTexture(player_s);
		player.setOrigin(32, 32); //center the origin of the player (half the with, half the height)
		float player_radius = player.getLocalBounds().width / 2.5; // radius of circle containing sprite; pass to collision function


		for (int i = 1; i < planets.size(); i++) { //set planet values
			planets[i].setRadius(radius[i]);
			planets[i].setTexture(&textures[i + 1]);
			planets[i].setOrigin(radius[i], radius[i]);
		}

		std::vector<vec_n> coordinates(10);

		planets[0].setTexture(&textures[0]);
		planets[0].setPosition(Vector2f(window2.getSize().x / 2, window2.getSize().y / 2));
		coordinates[0].x = planets[0].getPosition().x;
		coordinates[0].y = planets[0].getPosition().y;
		planets[0].setRadius(radius[0]);
		planets[0].setOrigin(radius[0], radius[0]);
	}
	//*/

	/*
	bool collision(std::vector<CircleShape>planets, std::vector<vec_n> coordinates, Sprite player, float player_radius) {

	Vector2f player_pos;
	player_pos.x = coordinates[0].x;
	player_pos.y = coordinates[0].y;
	bool collided = 0;

	for (int i = 1; i < coordinates.size; i++) {
	Vector2f planet_pos;
	planet_pos.x = coordinates[i].x;
	planet_pos.y = coordinates[i].y;

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
	}//*/


	void graph_init()
	{
		render_init();
	}

	void option_menue() {

		int header_pos_y[2] = { 400, 500 };
		Text option_header[2];
		option_header[0].setString("Sound");
		option_header[1].setString("BACK");
		for (int i = 0; i <= 1; i++) {
			option_header[i].setFont(font);
			option_header[i].setFillColor2(white);
			option_header[i].setCharacterSize(40);
			option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
		}


		RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
		underline.setFillColor(white);
		underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
		bool run_option = 1;
		bool selected[2] = { 0, 0 };

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


		int header_pos_y[3] = { 340, 440, 550 };
		Text option_header[3];
		option_header[0].setString("Sound");
		option_header[1].setString("Quit to Main menue");
		option_header[2].setString("Back");
		for (int i = 0; i <= 2; i++) {
			option_header[i].setFont(font);
			option_header[i].setFillColor2(white);
			option_header[i].setCharacterSize(50);
			option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
		}


		RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
		underline.setFillColor(white);
		underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
		bool run_option = 1;
		bool selected[3] = { 0, 0, 0 };

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
		title.setFillColor2(white);
		title.setPosition(Vector2f(window2.getSize().x / 2 - title.getLocalBounds().width*0.5, 200));

		int header_pos_y[3] = { 400, 500, 600 };



		Text option_header[3];
		option_header[0].setString("PLAY");
		option_header[1].setString("OPTION");
		option_header[2].setString("QUIT");
		for (int i = 0; i <= 2; i++) {
			option_header[i].setFont(font);
			option_header[i].setFillColor2(white);
			option_header[i].setCharacterSize(40);
			option_header[i].setPosition(Vector2f(window2.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
		}


		RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
		underline.setFillColor(white);
		underline.setPosition(Vector2f(window2.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));

		bool play = 0;
		bool selected[3] = { 0, 0, 0 };

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

	}

	void do_render()
	{
		main_render(phys::positions, 1, 1);
	}
}

#undef setFillColor