#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <windows.h>
#include "orbital_tools.h"
#include "render_tools.h"
#include "ui_tools.h"
#define skip if (false)

using namespace std;


namespace basics
{
	using shared::window2;

	void begin()
	{
		window2.clear();
		input::run_input();
	}
	void done()
	{
		window2.display();
	}
}

void main()
{

	srand(time(0));

	phys::engine_init();
	graph::graph_init();
	ui::ui_init();

	shared::s_time = clock();

	long long i = 0;
	long long t_cap = 5;
	while (true)
	{
		using namespace shared;
		r_time = clock() - shared::s_time;

		if (l_time + cps * 0.00001 < r_time)
		{
			basics::begin();

			phys::run_engine();

			graph::do_render();
#ifdef RENDER_DEBUG_INSTALLED
			render_debug::render_all();
#endif
			ui::run_ui();

			basics::done();

			l_time = r_time;
			i++;
		}

		if ((r_time - s_time) / cps > t_cap)
		{
			cout << "Ticks per second: " << i / t_cap << endl;
			t_cap += 5;
		}
	}
	return;
}