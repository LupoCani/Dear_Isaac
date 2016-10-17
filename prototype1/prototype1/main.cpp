#include <SFML/Graphics.hpp>
#include <math.h>
#include <Windows.h>
using namespace sf;
RenderWindow main_window(VideoMode(960, 520), "main window");

void menue_Render() {

}

void main_Render(CircleShape test) {
	Sleep(50);
	main_window.clear();
	main_window.draw(test);
	main_window.display();

}

int main()
{
	CircleShape p1;
	p1.setRadius(15);
	p1.setFillColor(Color(25, 5, 100));

	for (;;) {

		for (;;) {
			for (float i = 0; i <= 500; i = i + 0.2) {
				float y = 200*sin(i);
				p1.setPosition(Vector2f(i, y));
				main_Render(p1);
			}
		}


	}
}