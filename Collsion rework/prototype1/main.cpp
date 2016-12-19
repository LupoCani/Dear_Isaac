#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>

using namespace sf;

sf::RenderWindow window(sf::VideoMode(1080, 860), "Orbitals");
Color yellow(225, 237, 7);
Color Red(165, 41, 13);
Color grey(177, 190, 198);

bool collided = 0;


void collision(Sprite player,float player_radius, CircleShape sun, CircleShape planet1) {

	Vector2f player_center(player.getPosition().x + player.getLocalBounds().width / 2, player.getPosition().y + player.getLocalBounds().height / 2);
	float collision_box_begin = player_center.x - player_radius;
	float collision_box_end = collision_box_begin + player_radius * 2;

	float sun_radius = sun.getRadius();
	Vector2f sun_center(sun.getPosition().x + sun_radius, sun.getPosition().y + sun_radius);

	float planet1_radius = planet1.getRadius();
	Vector2f planet1_center(planet1.getPosition().x + planet1_radius, planet1.getPosition().y + planet1_radius);

	float precision = 1;
	for (float i = collision_box_begin; i <= collision_box_end; i += 1) {

		float rout_collision = sqrt(player_radius*player_radius - (i - player_center.x)*(i - player_center.x));
		float pos_collision = rout_collision + player_center.y;
		float neg_collision = -rout_collision + player_center.y;


		if (sun_radius*sun_radius - (i - sun_center.x)*(i - sun_center.x) >= 0) {

			float rout_sun_collision = sqrt(sun_radius*sun_radius - (i - sun_center.x)*(i - sun_center.x));
			float pos_sun_collision = rout_sun_collision + sun_center.y;
			float neg_sun_collision = -rout_sun_collision + sun_center.y;

			if (pos_sun_collision>=pos_collision-precision && pos_sun_collision <= pos_collision + precision)
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

	}

}

int main() {

	Texture planet_texture1;
	if (!planet_texture1.loadFromFile("planet texture.png")) {
		std::cout << "error loading planet texture" << std::endl;
	}

	Texture player_texture;
	if (!player_texture.loadFromFile("Character sprite.png")) {
		std::cout << "error loading player texture" << std::endl;
	}
	Sprite player;
	player.setTexture(player_texture);
	float player_radius = player.getLocalBounds().height/2;
	

	CircleShape sun(50);
	sun.setPointCount(50);
	sun.setFillColor(yellow);
	float sun_pos_y = 860 / 2 - sun.getRadius();
	float sun_pos_x = 1080 / 2 - sun.getRadius();
	sun.setPosition(Vector2f(sun_pos_x, sun_pos_y));

	CircleShape planet1(20);
	planet1.setTexture(&planet_texture1);

	while (window.isOpen()) {

		float k = 0;

		while (true) {

			float planet1_pos_x = 200 * cos(k/100)+sun_pos_x+ sun.getRadius();
			float planet1_pos_y = 198 * sin(k/100)+sun_pos_y+ sun.getRadius();
			planet1.setPosition(Vector2f(planet1_pos_x, planet1_pos_y));

			float player_pos_x=k;
			float player_pos_y= 200 * sin(k/100) + 500;
			player.setPosition(Vector2f(player_pos_x, player_pos_y));

			window.clear();
			window.draw(sun);
			window.draw(planet1);
			window.draw(player);
			window.display();

			collision(player,player_radius, sun, planet1);

			if (collided==1)
			{
				break;
			}

			k += 0.05;
		
		}


		while (true) {

			Event close;
			while (window.pollEvent(close)) {

				if ((close.type == Event::KeyPressed) && (close.key.code == Keyboard::Space))
				{
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