#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>
using namespace sf;
bool colision = 0;


sf::RenderWindow window(sf::VideoMode(1080, 860), "Orbitals");
Color blue(10, 5, 100);
Color green(5, 100, 10);
Color red(100, 5, 10);

void rend(CircleShape p1, CircleShape planet1, CircleShape sun) {

	Vector2f pos1(p1.getPosition());
	Vector2f pos2(planet1.getPosition());
	Vector2f pos3(sun.getPosition());
	float x1 = pos1.x + 10;
	float y1 = pos1.y + 10;
	float xmax = x1 + 10;
	float xmin = x1 - 10;

	float x2 = pos2.x + 20;
	float y2 = pos2.y + 20;
	float x3 = pos3.x + 30;
	float y3 = pos3.y + 30;

	float precision = 2;

	for (float k = xmin; k <= xmax; k = k + 0.001) {
		float pointp1 = sqrt(100 - (k - x1)*(k - x1)) + y1;
		float negpoint1 = -sqrt(100 - (k - x1)*(k - x1)) + y1;
		
		if (20 * 20 - (x2 - k)*(x2 - k) >= 0) {
		
			float point = sqrt(20 * 20 - (k - x2)*(k - x2)) + y2;
			float negpoint = -sqrt(20 * 20 - (k - x2)*(k - x2)) + y2;

			if (pointp1 >= point - precision && pointp1 <= point + precision) {
				colision = 1;
				break;
			}
			if (negpoint1 >= point + precision && negpoint1 <= point + precision) {
				colision = 1;
				break;
			}
			if (pointp1 >= negpoint - precision && pointp1 <= negpoint + precision) {
				colision = 1;
				break;
			}
		}

		if (30 * 30 - (x3 - k)*(x3 - k) >= 0) {

			float point = sqrt(30 * 30 - (k - x3)*(k - x3)) + y3;
			float negpoint = -sqrt(30 * 30 - (k - x3)*(k - x3)) + y3;

			if (pointp1 >= point - precision && pointp1 <= point + precision) {
				colision = 1;
				break;
			}
			if (negpoint1 >= point + precision && negpoint1 <= point + precision) {
				colision = 1;
				break;
			}
			if (pointp1 >= negpoint - precision && pointp1 <= negpoint + precision) {
				colision = 1;
				break;
			}
		}


	}

}

void main() {

		CircleShape p1(10);
		p1.setFillColor(blue);
		CircleShape planet1(20);
		planet1.setFillColor(red);
		CircleShape reference(30);
		reference.setFillColor(green);
		reference.setPosition(Vector2f(510, 400));


	for (float i = 0; i <= 1080; i = i + 0.5) {

		p1.setPosition(Vector2f(i*i / (1000 - i),200*cos(i/30)+400));
		planet1.setPosition(Vector2f((200+i)* cos(i/300)+500,200*sin(i/500)+390));

		window.clear();
		window.draw(p1);
		window.draw(planet1);
		window.draw(reference);
		window.display();
		Sleep(5);
		rend(p1, planet1, reference);
		if (colision == 1) {
			break;
		}
}

	for (;;) {

	}
}