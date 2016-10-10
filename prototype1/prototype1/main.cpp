#include <SFML/Graphics.hpp>

int main()
{


sf::RenderWindow window(sf::VideoMode(810, 540), "SFML works!"); //Render window


	sf::CircleShape shape(100.f);
	shape.setFillColor(sf::Color::Green);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed) {
				window.close();
				return 0;
			}
		}

		window.clear();
		window.draw(shape);
		window.display();
	}

	
}