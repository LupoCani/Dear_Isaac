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
Color white(200, 200, 200);


int main() {

	sf::Font font;
	if (!font.loadFromFile("ALGER.ttf"))
	{
		std::cout << "font not loaded";
	}

	Text title;
	title.setFont(font);
	title.setString("Dear Isaac");
	title.setCharacterSize(60);
	title.setFillColor(white);
	title.setPosition(Vector2f(540 - title.getLocalBounds().width*0.5, 200));

	int header_pos_y[3] = { 400,500,600 };

	Text option_header[3];
	option_header[0].setString("PLAY");
	option_header[1].setString("OPTION");
	option_header[2].setString("QUIT");
	for (int i = 0; i <= 2; i++) {
		option_header[i].setFont(font);
		option_header[i].setFillColor(white);
		option_header[i].setCharacterSize(30);
		option_header[i].setPosition(Vector2f(540 - option_header[i].getLocalBounds().width*0.5, header_pos_y[i]));
	}


	RectangleShape underline(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
	underline.setFillColor(white);
	underline.setPosition(Vector2f(540 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));

	bool selected[3] = { 0,0,0 };
	std::cout << underline.getPosition().y << std::endl << option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3;
	while (window.isOpen()) {
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
							underline.setPosition(Vector2f(540 - underline.getSize().x*0.5, option_header[i + 1].getPosition().y + option_header[i + 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[0].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(540 - underline.getSize().x*0.5, option_header[0].getPosition().y + option_header[0].getLocalBounds().height*1.3));
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
							underline.setPosition(Vector2f(540 - underline.getSize().x*0.5, header_pos_y[i - 1] + option_header[i - 1].getLocalBounds().height*1.3));
							break;
						}
						else {
							underline.setSize(Vector2f(option_header[2].getLocalBounds().width * 1.5, 2));
							underline.setPosition(Vector2f(540 - underline.getSize().x*0.5, header_pos_y[2] + option_header[2].getLocalBounds().height*1.3));
							break;
						}
					}
				}
			}

			if ((input.type == Event::KeyPressed) && (input.key.code != Keyboard::Up) && (input.key.code != Keyboard::Down)) {
				if (selected[2] == 1) {
					return 0;
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

