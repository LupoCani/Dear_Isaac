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


void collision(Sprite player, float player_radius, CircleShape sun, CircleShape planet1) {

	Vector2f player_pos;
	player_pos.x = player.getPosition().x;
	player_pos.y = player.getPosition().y;

	Vector2f planet_pos;
	planet_pos.x = planet1.getPosition().x;
	planet_pos.y = planet1.getPosition().y;

	float radius_compare = planet1.getRadius() + player_radius;

	Vector2f diference;
	diference.x = player_pos.x - planet_pos.x;
	diference.y = player_pos.y - planet_pos.y;

	if (sqrt(diference.x*diference.x+diference.y*diference.y)<=radius_compare)
	{
		collided = 1;
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
	player.setOrigin(32, 32);
	float player_radius = player.getLocalBounds().height / 2;



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

			float planet1_pos_x = 300 * cos(k / 100) + sun_pos_x + sun.getRadius();
			float planet1_pos_y = 295 * sin(k / 100) + sun_pos_y + sun.getRadius();
			planet1.setPosition(Vector2f(planet1_pos_x, planet1_pos_y));

			float player_pos_x = k;
			float player_pos_y = 200 * -sin(k / 10) + 500;
			player.setPosition(Vector2f(player_pos_x, player_pos_y));
			player.setRotation(2 * k);

			window.clear();
			window.draw(sun);
			window.draw(planet1);
			window.draw(player);
			window.display();

			collision(player, player_radius, sun, planet1);

			if (collided == 1)
			{
				break;
			}

			k += 0.01;

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