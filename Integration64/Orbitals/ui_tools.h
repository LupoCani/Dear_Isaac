#ifndef ORBITAL_TOOLS_LOADED
#include "orbital_tools.h"
#endif


/*
This is your house, so to speak. Any globals, "using namespace" and so forth are fine, so long as they're within the {} of namespace ui.

Any data you will need, you can get from shared::world_state::

If the data you need is not there, request it be put there.


Init will be called once, for starting the program, and declaring any variables you'll need to begin with.

Run will be called one hundred times per second, after the corresponding functions in render- and orbital tools.
*/

namespace ui
{
	using namespace sf;
	using namespace std;
	namespace ws = shared::world_state;

	void hp(double hp, double hp_max)
	{

	}

	void run_ui()
	{

	}

	void ui_init()
	{
		double hp = ws::health;
		double hp_max = ws::health_max;

		Font font;
		if (!font.loadFromFile("ALGER.ttf")) {
		}

		Text text("test", font);
		text.setCharacterSize(30);
		text.setStyle(Text::Bold);
		text.setColor(Color::Red);

		window draw.(text);
	}
}