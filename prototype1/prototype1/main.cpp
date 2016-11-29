#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>

using namespace sf;

sf::RenderWindow window(sf::VideoMode(1080, 860), "Orbitals");
Color yellow(255, 250, 7);
Color Red(165, 41, 13);
Color grey(177, 190, 198);

bool collision = 0;


void collide(CircleShape p1, CircleShape Sun, CircleShape planet1) {

	float Sun_r = Sun.getRadius();
	float P1_r = p1.getRadius();
	float planet1_r = planet1.getRadius();

	Vector2f Sun_center(Sun.getPosition());
	float sun_x = Sun_center.x + Sun_r;
	float sun_y = Sun_center.y + Sun_r;
	Vector2f p1_center(p1.getPosition());
	float p1_x = p1_center.x + P1_r;
	float p1_y = p1_center.y + P1_r;
	Vector2f planet1_center(planet1.getPosition());
	float planet1_x = planet1_center.x + planet1_r;
	float planet1_y = planet1_center.y + planet1_r;

	float precision = 2;

	for (float i = p1_center.x; i <= p1_x + P1_r; i += 1) {

		float p1_pos = sqrt(P1_r*P1_r - (i - p1_x)*(i - p1_x)) + p1_y;
		float p1_Neg = -sqrt(P1_r*P1_r - (i - p1_x)*(i - p1_x)) + p1_y;
		
		float pos;
		float Neg;

		if (Sun_r*Sun_r - (i - sun_x)*(i - sun_x) > 0) {

			 pos = sqrt(Sun_r*Sun_r - (i - sun_x)*(i - sun_x)) + sun_y;
			 Neg = -sqrt(Sun_r*Sun_r - (i - sun_x)*(i - sun_x)) + sun_y;

			if (p1_pos >= pos - precision && p1_pos <= pos + precision) {
				collision = 1;
				break;
			}

			if (p1_pos >= Neg - precision && p1_pos <= Neg + precision) {
				collision = 1;
				break;
			}

			if (p1_Neg >= pos - precision && p1_Neg <= pos + precision) {
				collision = 1;
				break;
			}

			if (p1_Neg >= Neg - precision && p1_Neg <= Neg + precision) {
				collision = 1;
				break;
			}

		}

		if (planet1_r*planet1_r - (i - planet1_x)*(i - planet1_x) >= 0) {

			 pos = sqrt(planet1_r*planet1_r - (i - planet1_x)*(i - planet1_x)) + planet1_y;
			 Neg = -sqrt(planet1_r*planet1_r - (i - planet1_x)*(i - planet1_x)) + planet1_y;

			 if (p1_pos >= pos - precision && p1_pos <= pos + precision) {
				 collision = 1;
				 break;
			 }

			 if (p1_pos >= Neg - precision && p1_pos <= Neg + precision) {
				 collision = 1;
				 break;
			 }

			 if (p1_Neg >= pos - precision && p1_Neg <= pos + precision) {
				 collision = 1;
				 break;
			 }

			 if (p1_Neg >= Neg - precision && p1_Neg <= Neg + precision) {
				 collision = 1;
				 break;
			 }

		}


	}

}

void main() {

	CircleShape p1(10);
	p1.setFillColor(grey);
	CircleShape planet1(20);
	planet1.setFillColor(Red);
	
	CircleShape Sun(35);
	Sun.setFillColor(yellow);
	Sun.setPosition(505, 395);

	while (window.isOpen()) {

		float k = 0;

		for (;;) {

			p1.setPosition(Vector2f(0.05*k, ( 0.5*k) * cos(k / 100) + 420));
			planet1.setPosition(Vector2f(300 * sin(k / 100) + 520, 300 * cos(k / 100) + 410));

			window.clear();
			window.draw(Sun);
			window.draw(p1);
			window.draw(planet1);
			window.display();
			collide(p1, Sun, planet1);

			Sleep(5);
			k += 0.1;

			if (collision == 1) {
				break;
			}

		}

		for (;;) {

		}

	}



}