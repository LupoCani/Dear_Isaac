#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>

using namespace sf;

sf::RenderWindow window(sf::VideoMode(1920, 1080), "Orbitals");
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

struct body
{
	double ax_a;    //Semimajor axis
	double ax_b;    //Semiminor axis
	double ecc;     //Orbital eccentrictity
	double t_0;     //Starting timestamp
	double t_p;
	double t_l;     //Timestamp at a given time
	double u;       //Gravitational parameter
	double normal;  //angle of the major axis
	double area;    //Orbital area
	double SOI;
	bool closed;
	bool inverse;
	vec_n vel;      //Velocity at a given time
	vec_n pos;      //Position at a given time

	double En;      //Orbital energy
	double Ar;      //Area swept per unit of time
	double Mn;      //Mean anomaly swept per unit of time

	body* parent;
};

void option_menue() {

	int header_pos_y[2] = { 400,500 };
	Text option_header[2];
	option_header[0].setString("Sound");
	option_header[1].setString("BACK");
	for (int i = 0; i <= 1; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(40);
		option_header[i].setPosition(Vector2f(window.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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

		while (window.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 1; i++) {
					if (selected[i] == 1) {
						if (i == 0) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[1] + option_header[1].getLocalBounds().height*1.3));
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

		window.clear();
		for (int i = 0; i <= 1; i++) {
			window.draw(option_header[i]);
		}
		window.draw(underline);
		window.display();

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
		option_header[i].setPosition(Vector2f(window.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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

		while (window.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i <= 1) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[2].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[2] + option_header[2].getLocalBounds().height*1.3));
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

		window.clear();
		for (int i = 0; i <= 2; i++) {
			window.draw(option_header[i]);
		}
		window.draw(underline);
		window.display();

	}

}

void viewport_render(Sprite player, CircleShape sun, CircleShape planet1, CircleShape planet2) {

	float modifi_x = player.getPosition().x - window.getSize().x / 2;
	float modifi_y = player.getPosition().y - window.getSize().y / 2;


	Vector2f player_newpos(player.getPosition().x - modifi_x, player.getPosition().y - modifi_y);
	Vector2f sun_newpos(sun.getPosition().x - modifi_x, sun.getPosition().y - modifi_y);
	Vector2f planet1_newpos(planet1.getPosition().x - modifi_x, planet1.getPosition().y - modifi_y);
	Vector2f planet2_newpos(planet2.getPosition().x - modifi_x, planet2.getPosition().y - modifi_y);

	player.setPosition(Vector2f(player_newpos.x, player_newpos.y));
	sun.setPosition(Vector2f(sun_newpos.x, sun_newpos.y));
	planet1.setPosition(Vector2f(planet1_newpos.x, planet1_newpos.y));
	planet2.setPosition(Vector2f(planet2_newpos.x, planet2_newpos.y));
	window.clear();
	window.draw(sun);
	window.draw(planet1);
	window.draw(planet2);
	window.draw(player);
	window.display();
}

void collision(Sprite player, float player_radius, CircleShape sun, CircleShape planet1, CircleShape planet2) {

	Vector2f player_center(player.getPosition().x, player.getPosition().y);
	float collision_box_begin = player_center.x - player_radius;
	float collision_box_end = collision_box_begin + player_radius * 2;

	float sun_radius = sun.getRadius();
	Vector2f sun_center(sun.getPosition().x + sun_radius, sun.getPosition().y + sun_radius);

	float planet1_radius = planet1.getRadius();
	Vector2f planet1_center(planet1.getPosition().x + planet1_radius, planet1.getPosition().y + planet1_radius);

	float planet2_radius = planet2.getRadius();
	Vector2f planet2_center(planet2.getPosition().x + planet2_radius, planet2.getPosition().y + planet2_radius);

	float precision = 1;
	for (int i = collision_box_begin; i <= collision_box_end; i += 1) {

		float rout_collision = sqrt(player_radius*player_radius - (i - player_center.x)*(i - player_center.x));
		float pos_collision = rout_collision + player_center.y;
		float neg_collision = -rout_collision + player_center.y;


		if (sun_radius*sun_radius - (i - sun_center.x)*(i - sun_center.x) >= 0) {

			float rout_sun_collision = sqrt(sun_radius*sun_radius - (i - sun_center.x)*(i - sun_center.x));
			float pos_sun_collision = rout_sun_collision + sun_center.y;
			float neg_sun_collision = -rout_sun_collision + sun_center.y;

			if (pos_sun_collision >= pos_collision - precision && pos_sun_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (pos_sun_collision >= neg_collision - precision && pos_sun_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_sun_collision >= pos_collision - precision && neg_sun_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_sun_collision >= neg_collision - precision && neg_sun_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}

		}

		if (planet1_radius*planet1_radius - (i - planet1_center.x)*(i - planet1_center.x) >= 0) {

			float rout_planet1_collision = sqrt(planet1_radius*planet1_radius - (i - planet1_center.x)*(i - planet1_center.x));
			float pos_planet1_collision = rout_planet1_collision + planet1_center.y;
			float neg_planet1_collision = -rout_planet1_collision + planet1_center.y;

			if (pos_planet1_collision >= pos_collision - precision && pos_planet1_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (pos_planet1_collision >= neg_collision - precision && pos_planet1_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_planet1_collision >= pos_collision - precision && neg_planet1_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_planet1_collision >= neg_collision - precision &&  neg_planet1_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}
		}

		if (planet2_radius*planet2_radius - (i - planet2_center.x)*(i - planet2_center.x) >= 0) {

			float rout_planet2_collision = sqrt(planet2_radius*planet2_radius - (i - planet2_center.x)*(i - planet2_center.x));
			float pos_planet2_collision = rout_planet2_collision + planet2_center.y;
			float neg_planet2_collision = -rout_planet2_collision + planet2_center.y;

			if (pos_planet2_collision >= pos_collision - precision && pos_planet2_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (pos_planet2_collision >= neg_collision - precision && pos_planet2_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_planet2_collision >= pos_collision - precision && neg_planet2_collision <= pos_collision + precision)
			{
				collided = 1;
				break;
			}
			if (neg_planet2_collision >= neg_collision - precision &&  neg_planet2_collision <= neg_collision + precision)
			{
				collided = 1;
				break;
			}
		}

	}

}

CircleShape B_to_C(body in)
{
	double pos_x = in.pos.x;
	double pos_y = in.pos.y;
}


int start_menue() {
	Text title;
	title.setFont(font);
	title.setString("Dear Isaac");
	title.setCharacterSize(60);
	title.setFillColor(white);
	title.setPosition(Vector2f(window.getSize().x / 2 - title.getLocalBounds().width*0.5, 200));

	int header_pos_y[3] = { 400,500,600 };

	Text option_header[3];
	option_header[0].setString("PLAY");
	option_header[1].setString("OPTION");
	option_header[2].setString("QUIT");
	for (int i = 0; i <= 2; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(40);
		option_header[i].setPosition(Vector2f(window.getSize().x / 2 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));

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

		while (window.pollEvent(input)) {

			if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
				for (int i = 0; i <= 2; i++) {
					if (selected[i] == 1) {
						if (i <= 1) {
							underline.setSize(Vector2f(option_header[i + 1].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[2].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(window.getSize().x / 2 - underline.getSize().x*0.5, header_pos_y[2] + option_header[2].getLocalBounds().height*1.3));
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
					window.close();
					return 0;
					break;
				}
				if (selected[1] == 1) {
					option_menue();
				}
			}

		}
		window.clear();
		for (int i = 0; i <= 2; i++) {
			window.draw(option_header[i]);
		}
		window.draw(title);
		window.draw(underline);
		window.display();


	}

}

int main() {



	if (!font.loadFromFile("Alger.ttf"))
	{

	}


	Texture planet_texture1;
	if (!planet_texture1.loadFromFile("planet texture.png")) {

	}

	Texture planet_texture2;
	if (!planet_texture2.loadFromFile("planet texture2.png")) {

	}

	Texture player_texture;
	if (!player_texture.loadFromFile("Character sprite.png")) {

	}
	Sprite player;
	player.setTexture(player_texture);
	player.setOrigin(32, 32);
	float player_radius = player.getLocalBounds().height / 2.5;



	CircleShape sun(75);
	sun.setFillColor(yellow);
	float sun_pos_y = 860 / 2 - sun.getRadius();
	float sun_pos_x = 1080 / 2 - sun.getRadius();
	sun.setPosition(Vector2f(sun_pos_x, sun_pos_y));

	std::vector<body*> bodies;

	CircleShape planet1(40);
	planet1.setTexture(&planet_texture1);

	CircleShape planet2(55);
	planet2.setTexture(&planet_texture2);


	while (window.isOpen()) {
		restart = 0;
		start_menue();
		float k = 0;
		float r = 1000;
		player.setRotation(0);
		collided = 0;
		while (restart == 0) {

			float planet1_pos_x = 300 * cos(k / 100) + sun_pos_x + sun.getRadius();
			float planet1_pos_y = 295 * sin(k / 100) + sun_pos_y + sun.getRadius();
			planet1.setPosition(Vector2f(planet1_pos_x, planet1_pos_y));

			float planet2_pos_x = 410 * cos(k / 50) + sun_pos_x + sun.getRadius();
			float planet2_pos_y = 409 * sin(k / 50) + sun_pos_y + sun.getRadius();
			planet2.setPosition(Vector2f(planet2_pos_x, planet2_pos_y));

			float player_pos_x = (r - 0.4*k)*cos(k / 50) + sun_pos_x;
			float player_pos_y = (r - 0.5*k)*sin(k / 50) + sun_pos_y;
			player.setPosition(Vector2f(player_pos_x, player_pos_y));

			Event input;

			while (window.pollEvent(input)) {

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Right)) {
					player.rotate(1);
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Left)) {
					player.rotate(-1);
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Up)) {
					r += 1;
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Down)) {
					r -= 1;
				}

				if ((input.type == Event::KeyPressed) && (input.key.code == Keyboard::Escape)) {
					ingame_menue();
				}

				if (input.type == Event::Closed) {
					window.close();
					return 0;
				}
			}

			viewport_render(player, sun, planet1, planet2);

			collision(player, player_radius, sun, planet1, planet2);

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
				while (window.pollEvent(close)) {

					if ((close.type == Event::KeyPressed) && (close.key.code == Keyboard::Space))
					{
						end = 0;
						break;
					}
					if (close.type == Event::Closed)
					{
						window.close();
						return 0;
					}

				}

			}
		}

	}



}