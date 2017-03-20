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


void main()
{
	srand(time(0));

	phys::engine_init();
	graph::graph_init();

	shared::s_time = clock();

	long long i = 0;
	long long t_cap = 5;

	while (true)
	{
		using namespace shared;
		r_time = clock() - shared::s_time;

		if (l_time + cps * 0.00001 < r_time)
		{
			input::run_input();

			phys::run_engine();

#ifdef RENDER_DEBUG_INSTALLED
			render_debug::render_all();
#endif
			graph::do_render();

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

	/*
	//fill_tail();
	using namespace phys;
	int steps = 1000;
	int render_subdiv = 5;

	vec_n pos_temp, vel_temp;

	srand(time(0));

	cout << "Drawing. \n";

	for (int i2 = 0; i2 < bodies.size() - 1; i2++)
	{
		body &sat = *bodies[i2];

		for (int i = 0; i < int(steps / render_subdiv); i++)
		{
			double angle = to_rad(i, steps / render_subdiv);
			double radius = get_r(angle, sat);

			vec_n dot;
			dot.y = sin(angle) * radius;
			dot.x = cos(angle) * radius;

			dot = rot_vec(dot, sat.normal);
			//dot = t_to_screen(dot);

			draw_tail(dot.x, dot.y, false);
		}
	}

	draw_tail(    0,  200, true);
	draw_tail(    0, -200, true);
	draw_tail( -200,    0, true);
	draw_tail(  200,    0, true);

	draw_tail(0, 0, true);

	draw_tail(500, 500, true);
	graph::window2.display();

	cout << endl;
	long f_time_cache = 0;
	double zoom = 0.002;
	vec_n c_buff;
	vector<vec_n> positions;
	int focus = 0;

	bool eml = false;

	double e_time_last = 0;

	vector<int> tick_counts(0);
	int tick_count = 0;
	vector<double> tick_times(0);
	double tick_time = clock();

	do
	{
		double d_time = w_time;
		long f_time = clock();

		keys inp_in = render_cst();
		double fact = 1000;
		bool e_mode = inp_in.x + inp_in.y > 0.5;
		zoom *= pow(1.05, inp_in.z);
		if (inp_in.lc)
		{
			cin >> focus;
		}

		if (e_mode)
		{
			e_time_last = f_time;
			eml = true;
		}


		if ((e_time_last + 100 < f_time) && eml)
		{
			eml = false;
			cout << "Ecc:   " << plyr.ecc << endl;
			cout << "Ax_a:  " << plyr.ax_a << endl;	//Semimajor axis
			cout << "Ax_b:  " << plyr.ax_b << endl;	//Semiminor axis
			cout << "epoch: " << plyr.epoch << endl;		//Starting timestamp
			cout << "t_l:   " << plyr.t_l << endl;		//Timestamp at a given time
			cout << "u:     " << plyr.u << endl;		//Gravitational parameter
			cout << "norm:  " << plyr.normal << endl;	//angle of the major axis
			cout << "SOI:   " << yavin.SOI << endl;
			cout << "shape: " << plyr.shape << endl;
			cout << "inv:   " << plyr.inverse << endl;
			cout << "vel_x: " << plyr.vel.x << endl;		//Velocity at a given time
			cout << "vel_y: " << plyr.vel.y << endl;
			cout << "pos_x: " << plyr.pos.x << endl;
			cout << "pos_y: " << plyr.pos.y << endl;

			cout << "En:    " << plyr.En << endl;
			cout << "Ar:    " << plyr.Ar << endl;
			cout << "Mn:    " << plyr.Mn << endl;
			cout << "Name:  " << plyr.name << endl;
			cout << "Parent:" << (*(plyr.parent)).name << endl;

			cout << "___________________________\n__________________________\n";


			while (true)
			{
				if ((*(tail.end() - 1)).getFillColor() == sf::Color::Red)
				{
					tail.erase(tail.end() - 1);
					tail_coord.erase(tail_coord.end() - 1);

				}
					
				else break;
			}

			do_orbit(plyr, w_time);

			body sat = plyr;
			int dotcount = 200;
			for (int i = 0; i < dotcount; i++)
			{
				double angle = to_rad(i, dotcount);
				double radius = get_r(angle, sat);

				vec_n dot;
				dot.y = sin(angle) * radius;
				dot.x = cos(angle) * radius;

				dot = rot_vec(dot, sat.normal);
				//dot = t_to_screen(dot);

				draw_tail(dot.x, dot.y, true);
			}
		}

		positions = update_all_and_convert(bodies, w_time, e_mode, inp_in.x, inp_in.y, fact);

		graph::main_render();/*

		focus = focus % positions.size();

		vec_n p_pos_old = positions[focus];

		if (0 || (f_time_cache + 50) < f_time)
		{ 
			positions = graph::handle_scale(positions, focus, zoom, 500, 500);
			window2.clear();
			for (int i = 0; i < positions.size(); i++)
			{
				vec_n sat_pos = positions[i];
				c_buff = sat_pos;
				draw_circle(c_buff.x, c_buff.y, i, zoom);
			}

			draw_tail(0, 0, 0, true, 0, zoom, p_pos_old, (*(plyr.parent)).pos );


			sf::Text FPS;
			FPS.setFont(graph::font);
			FPS.setPosition(100, 100);
			FPS.setCharacterSize(25);
			//FPS.setString(to_string(countdown));
			window2.draw(FPS);

			sf::Text FPS2;
			FPS2.setFont(graph::font);
			FPS2.setPosition(100, 150);
			FPS2.setCharacterSize(25);
			//FPS2.setString(to_string(expire));
			window2.draw(FPS2);

			window2.display();
			f_time_cache = f_time;
		}//*/
		/*


		w_time += 0.0001;
		Sleep(1);

	} while (true);

	std::system("pause");
	std::system("pause"); */
}