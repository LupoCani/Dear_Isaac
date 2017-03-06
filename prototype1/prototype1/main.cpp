#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
#include <ctime>

using namespace sf;



sf::RenderWindow window2(sf::VideoMode(810*2, 540*2), "Orbitals");
Color yellow(225, 237, 7);
Color Red(165, 41, 13);
Color grey(177, 190, 198);
Color white(255, 255, 255);
bool collided = 0;
bool restart = 0;

sf::Font font;

struct vec_n
{
	double x;
	double y;
	double z;
};

Vector2f generate_goal(std::vector<vec_n>cordinats, int sun_r, int radius) {

	Vector2f goal_cordinats;

	srand(time(0));

	for (;;) {

		goal_cordinats.x = rand() % (2 * radius) + cordinats[cordinats.size() - 1].x - radius;

		if (goal_cordinats.x < cordinats[cordinats.size()-1].x - sun_r) {
			break;
		}
		else if (goal_cordinats.x > cordinats[cordinats.size()-1].x + sun_r) {
			break;
		}

	}

	for (;;) {

		goal_cordinats.y = rand() % (2 * radius) + cordinats[cordinats.size() - 1].y - radius;

		if (goal_cordinats.y < cordinats[cordinats.size()-1].y - sun_r) {
			break;
		}
		else if (goal_cordinats.y > cordinats[cordinats.size()-1].y + sun_r) {
			break;
		}

	}
	std::cout << goal_cordinats.x << std::endl << goal_cordinats.y << std::endl;
	return(goal_cordinats);

}

bool goal_collision(CircleShape goal, Vector2f goal_cordinats, std::vector<vec_n>cordinats, float player_r) {

	bool goal_collided=0;

	float distance = sqrt((cordinats[cordinats.size() - 1].x - goal_cordinats.x)*(cordinats[cordinats.size() - 1].x - goal_cordinats.x) + (cordinats[cordinats.size() - 1].y - goal_cordinats.y)*(cordinats[cordinats.size() - 1].y - goal_cordinats.y));

	if (distance < player_r + goal.getRadius()) {
		goal_collided = 1;
	}

	return(goal_collided);

}

void option_menue() {

	int header_pos_y[2] = { window2.getSize().y/3, header_pos_y[0]+100};
	Text option_header[2];
	option_header[0].setString("-Sound+");
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



	int header_pos_y[3] = { window2.getSize().y / 3, header_pos_y[0] + 100,header_pos_y[0] + 200 };
	Text option_header[3];
	option_header[0].setString("-Sound+");
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
					restart = 1;
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

void main_render(std::vector<vec_n> cordinats, std::vector<CircleShape> planets, Sprite player, CircleShape goal, Vector2f goal_cordinats) {

	Vector2f viewport_center; //cordinates at the center of the viewport/window
	viewport_center.x = window2.getSize().x / 2;
	viewport_center.y = window2.getSize().y / 2;

	Vector2f modifi_cordinates; //value to modifi curent cordinates with on oder to put player att the center of the 
	modifi_cordinates.x = cordinats[cordinats.size()-1].x - viewport_center.x;
	modifi_cordinates.y = cordinats[cordinats.size() - 1].y - viewport_center.y;

	Vector2f goal_new = goal_cordinats;

	///*
	for (int i = 0; i < cordinats.size(); i++) {
		cordinats[i].x -= modifi_cordinates.x;
		cordinats[i].y -= modifi_cordinates.y;
	}

	goal_new.x -= modifi_cordinates.x;
	goal_new.y -= modifi_cordinates.y;

	//*/
	//set new position based on calculated values 
	player.setPosition(Vector2f(cordinats[cordinats.size() - 1].x, cordinats[cordinats.size() - 1].y));
	for (int i = 0; i < 9; i++) {
		planets[i].setPosition(Vector2f(cordinats[i].x, cordinats[i].y));
	}
	goal.setPosition(goal_new);
	window2.clear();

	for (int i = 0; i < planets.size(); i++) {
		window2.draw(planets[i]);
	}
	
	window2.draw(goal);

	window2.draw(player);

	window2.display();



}


bool collision(std::vector<CircleShape>planets, std::vector<vec_n> cordinats, Sprite player, float player_radius) {

	Vector2f player_pos;
	player_pos.x = cordinats[cordinats.size()-1].x;
	player_pos.y = cordinats[cordinats.size() - 1].y;
	bool collided = 0;

	for (int i = 0; i < cordinats.size()-1; i++) {
		Vector2f planet_pos;
		planet_pos.x = cordinats[i].x;
		planet_pos.y = cordinats[i].y;

		float radius_compare = planets[i].getRadius() + player_radius;

		Vector2f diference;
		diference.x = player_pos.x - planet_pos.x;
		diference.y = player_pos.y - planet_pos.y;

		if (sqrt(diference.x*diference.x + diference.y*diference.y) <= radius_compare)
		{
			collided = 1;
			restart = 1;
			break;
		}
	}

	return (collided);
}





int start_menue() {

	Text title;
	title.setFont(font);
	title.setString("Dear Isaac");
	title.setCharacterSize(60);
	title.setFillColor(white);
	title.setPosition(Vector2f(window2.getSize().x / 2 - title.getLocalBounds().width*0.5, window2.getSize().x/7));

	int header_pos_y[3] = { title.getPosition().y + 150, title.getPosition().y + 250, title.getPosition().y +350};

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

}

int main() {

	if (!font.loadFromFile("ALGER.ttf")) {
		//handle exception
	}
	float player_scale;
	std::vector<CircleShape>planets(9); //pass to render and  collision function
	Sprite player; //pass to render and  collision function
	int modifi_scale = 10;
	int radius[9] = { 20 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale, 10 * modifi_scale }; //radius for the planets 
	float scale = 0.5;
	CircleShape goal(30);
	goal.setOrigin(15, 15);

	Texture planet_textures[9];
	if (!planet_textures[0].loadFromFile("sun_texture.png")) {

	}
	if (!planet_textures[1].loadFromFile("planet_texture2.png")) {

	}
	if (!planet_textures[2].loadFromFile("planet_texture3.png")) {

	}
	if (!planet_textures[3].loadFromFile("planet_texture4.png")) {

	}
	if (!planet_textures[4].loadFromFile("planet_texture5.png")) {

	}
	if (!planet_textures[5].loadFromFile("planet_texture6.png")) {

	}
	if (!planet_textures[6].loadFromFile("planet_texture7.png")) {

	}
	if (!planet_textures[7].loadFromFile("planet_texture8.png")) {

	}
	if (!planet_textures[8].loadFromFile("planet_texture9.png")) {

	}


	Texture player_texture;
	if (!player_texture.loadFromFile("Character_sprite.png")) {
		//handle exception
	}
	player.setScale(1 * scale, 1 * scale);
	player.setTexture(player_texture);
	player.setOrigin(32, 32); //center the origin of the player (half the with, half the height)
	float player_radius = player.getLocalBounds().width / 2.5; // radius of circle containing sprite; pass to collision function


	for (int i = 1; i < planets.size(); i++) { //set planet values
		planets[i].setRadius(radius[i]);
		planets[i].setTexture(&planet_textures[i]);
		planets[i].setOrigin(radius[i], radius[i]);
	}

	std::vector<vec_n> cordinats(10);

	planets[0].setTexture(&planet_textures[0]);
	planets[0].setPosition(Vector2f(window2.getSize().x / 2, window2.getSize().y / 2));
	cordinats[0].x = planets[0].getPosition().x;
	cordinats[0].y = planets[0].getPosition().y;
	planets[0].setRadius(radius[0]);
	planets[0].setOrigin(radius[0], radius[0]);
	int scop = 2000;
	Vector2f goal_cordinats(generate_goal(cordinats, radius[8]*1.5, scop));
	goal.setPosition(goal_cordinats);

	while (window2.isOpen()) {
		restart = 0;
		
		if (start_menue() == 0) {
			return 0;
		}
		float k = 0;
		float r = 30;
		player.setRotation(0);
		collided = 0;
		while (restart == 0) {



			for (int i = 1; i < 9; i++) {
				cordinats[i].x = i * 50 * modifi_scale * cos(k / 50+10*i) + planets[0].getPosition().x;
				cordinats[i].y = i * 50 * modifi_scale* sin(k / 50+10 * i) + planets[0].getPosition().y;
			}



			cordinats[cordinats.size()-1].x = (r)* modifi_scale*cos(k / 50+35) + planets[0].getPosition().x;
			cordinats[cordinats.size() - 1].y = (r)* modifi_scale*sin(k / 50+35) + +planets[0].getPosition().y;
			player.setPosition(Vector2f(cordinats[cordinats.size() - 1].x, cordinats[cordinats.size() - 1].y));



			Event input;

			while (window2.pollEvent(input)) {

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Right)) {
					player.rotate(1);
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Left)) {
					player.rotate(-1);
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
					r += 10;
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
					r -= 10;
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Escape)) {
					ingame_menue();
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Space)) {
					std::cin >> scale;
					player.setScale(scale,scale);
					for (int i = 0; i < planets.size(); i++) {
						planets[i].setScale(scale, scale);
					}
				
				}

				if (input.type == Event::Closed) {
					window2.close();
					return 0;
				}
			}

			if (goal_collision(goal,goal_cordinats,cordinats,player_radius)==1){
				goal_cordinats = generate_goal(cordinats, radius[8] * 1.5, scop);
			}
		//collision(planets, cordinats, player, player_radius);
			main_render(cordinats, planets, player, goal, goal_cordinats);
			if (collided == 1)
			{
				break;
			}

			k += 0.01;

		}
		bool end = 1;
		if (restart == 0) {
			while (end == 1) {

				Event close;
				while (window2.pollEvent(close)) {

					if ((close.type == Event::KeyPressed) && (close.key.code == Keyboard::Space))
					{
						end = 0;
						break;
					}
					if (close.type == Event::Closed)
					{
						window2.close();
						return 0;
					}

				}

			}
		}

	}



}