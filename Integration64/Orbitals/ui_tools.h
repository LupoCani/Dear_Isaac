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
	namespace ws = shared::world_state;
	using shared::window2;

	using std::vector;
	using std::string;
	using shared::vec_n;
	Font font;

	vector<string> men_items;
	string men_head;
	short men_sel = 0;
	short men_l;

	vector<string> pause_items;
	string pause_head;
	short pause_sel = 0;
	short pause_l;


	vector<string> crep_items;
	vector<double> crep_dists(5);

	void hp(double hp, double hp_max)
	{
	}

	vec_n scale_single(vec_n pos, vec_n origo, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		pos -= origo;
		pos *= scale;
		pos.y *= -1;

		pos += vec_n(mid_x, mid_y);

		return pos;
	}

	vector<vec_n> handle_scale(vector<vec_n> list, vec_n origo, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		for (int i = 0; i < list.size(); i++)
			list[i] = scale_single(list[i], origo, scale, mid_x, mid_y);

		return list;
	}

	VertexArray render_dmg_ball(vec_n origo, double rad, int count = 30)
	{
		VertexArray lines(TrianglesFan, count + 2);

		lines[0].position = origo;
		lines[0].color = Color(255, 0, 0, 150);

		for (int i = 0; i < count + 1; i++)
		{
			phys::vec_r pos;

			pos.mag = rad;
			pos.ang = double(i) / count * phys::M_2PI;

			pos = pos + origo;

			lines[i + 1].position = vec_n(pos);
			lines[i + 1].color = Color::Transparent;
		}
		return lines;
	}

	void draw_kesslers(RenderWindow &window, vector<vec_n> pos, double size, vec_n origo, double scale = 1, vec_n center = vec_n(500, 500))
	{
		pos = handle_scale(pos, origo, scale);// , center.x, center.y);

		for (int i = 0; i < pos.size(); i++)
		{
			window.draw(render_dmg_ball(pos[i], scale * size));
		}
	}

	void draw_crep(RenderWindow &window, vector<double> dists, vector<string> items, vector<double> values)
	{
		for (int i = 0; i < 3; i++)
		{
			Text title;
			RectangleShape box;
			RectangleShape line;
			RectangleShape end;

			float mrg = dists[0];
			float x_col = mrg + dists[3];
			float x_end = x_col + dists[4];
			float y = mrg + dists[2] * i;

			title.setFont(font);
			title.setCharacterSize(dists[1]);
			title.setString(items[i]);

			box.setSize(vec_n(dists[4] * values[i], dists[2] * 0.9));
			line.setSize(vec_n(dists[4], dists[2] * 0.2));
			end.setSize(vec_n(2, dists[2] * 0.9));

			FloatRect text_size = title.getLocalBounds();
			RectangleShape back_rec;
			title.setOrigin(vec_n(0, (text_size.top + text_size.height) / 2));

			box.setOrigin(vec_n(0, box.getSize().y / 2.0));
			line.setOrigin(vec_n(0, line.getSize().y / 2.0));
			end.setOrigin(vec_n(0, end.getSize().y / 2.0));

			title.setPosition(vec_n(mrg, y + dists[2] / 2));
			box.setPosition(vec_n(x_col, y + dists[2] / 2));
			line.setPosition(vec_n(x_col, y + dists[2] / 2));
			end.setPosition(vec_n(x_end, y + dists[2] / 2));

			window.draw(title);
			window.draw(box);
			window.draw(line);
			window.draw(end);
		}
	}

	vec_n true_size(FloatRect rect)
	{
		vec_n out;
		out.x = rect.width + rect.left;
		out.y = rect.height + rect.top;
		return out;
	}

	string to_string(double in, short max)
	{
		string out;
		if (isnan(in))
			out = "N/A";
		else
			out = std::to_string(in);

		if (out.length() > max)
			out = out.substr(0, max);

		return out;
	}

	void draw_plnbox(RenderWindow &window, short kind, vec_n pos, string name, double vel, double dist = 0, double pr_dist = 0, double countdown = 0)
	{
		vec_n box_size(200, 300);
		vec_n box_pos;
		vec_n box_orig = box_size * -1;
		vec_n box_move = box_orig;
		short hl_count = 1;
		if (kind == 1 or kind == 0)
			hl_count = 2;

		if (pos.x < 0)
			pos.x += window.getSize().x;
		else
			box_move.x = 0;
		if (pos.y < 0)
			pos.y += window.getSize().y;
		else
			box_move.y = 0;

		box_pos = pos;

		RectangleShape box(box_size);
		box.setFillColor(Color(100, 100, 100));
		box_pos += box_move;
		box.setPosition(box_pos);
		window.draw(box);

		/*
		Box Types:
		0: target, not intercepting.
		1: target, intercepting.
		2: not target, intercepting
		*/

		vector<string> items;
		if (kind == 0 or kind == 1)
			items.push_back("Target");
		if (kind == 0)
			items.push_back("I");
		if (kind == 1 or kind == 2)
			items.push_back("INTERCEPTING");
		{
			items.push_back("Name:     " + name);
			items.push_back("Distance: " + to_string(dist, 10));
			items.push_back("Velocity: " + to_string(vel, 10));
		}
		if (kind == 0 || kind == 1 || kind == 2)
		{
			items.push_back("Closest:  " + to_string(pr_dist, 10));
			items.push_back("ETA:      " + to_string(countdown, 10));
		}

		double l_y = 0;
		for (int i = 0; i < items.size(); i++)
		{
			Text item;
			item.setString(items[i]);
			item.setFont(font);
			item.setPosition(box_pos);
			item.move(vec_n(5, l_y));
			item.setCharacterSize(15);
			double x_cent = box_orig.x * -0.5;

			if (i < hl_count)
			{
				item.setCharacterSize(30);
				item.setFillColor(Color(255, 255, 255));
				if (i >= 1)
				{
					item.setCharacterSize(20);
					item.setFillColor(Color(255, 127, 0));
					if (kind == 0)
						item.setFillColor(Color::Transparent);
				}
				if (kind == 2)
				{
					item.setCharacterSize(25);
					item.setFillColor(Color(255, 127, 0));
				}
				vec_n center(true_size(item.getLocalBounds()));
				center *= 0.5;
				center.y = 0;

				item.setOrigin(center);
				item.move(vec_n(x_cent, 0));
			}
			l_y += true_size(item.getLocalBounds()).y + 3;
			window.draw(item);

			if (not (i < hl_count))
			{
				VertexArray lines(sf::LinesStrip, 2);
				lines[0].position = vec_n(box_pos.x,              box_pos.y + l_y);
				lines[1].position = vec_n(box_pos.x + box_size.x, box_pos.y + l_y);

				lines[0].color = Color(50, 50, 50);
				lines[1].color = Color(50, 50, 50);

				window.draw(lines);
			}
		}
	}

	void run_ui()
	{
		using shared::game_state;
		Vector2u w_size_sfu = window2.getSize();
		vec_n w_size(w_size_sfu.x, w_size_sfu.y);
		input::flush_back::men_cmd = 0;

		if (input::win_inf::update)
		{
			shared::window2.setView(View(FloatRect(vec_n(), Vector2f(window2.getSize()))));
			input::win_inf::update = false;
		}

		if (game_state == 0)
		{
			vector<Text> lines;
			vec_n h_pos = w_size;

			{
				Text line;
				line.setFont(font);
				line.setString(men_head);
				line.setCharacterSize(60);

				h_pos.x *= 0.1;
				h_pos.y *= 0.3;
				line.setPosition(h_pos);

				lines.push_back(line);
			}


			for (int i = 0; i < men_l; i++)
			{
				Text line;
				line.setFont(font);
				line.setString(men_items[i]);
				line.setCharacterSize(30);

				vec_n l_pos = h_pos;
				l_pos.y += w_size.y * 0.1 * (i + 1);
				line.setPosition(l_pos);

				lines.push_back(line);
			}

			if (input::keyboard.wasPressed(input::key_state::keys::Up))
				men_sel--;
			if (input::keyboard.wasPressed(input::key_state::keys::Down))
				men_sel++;

			if (men_sel >= men_l)
				men_sel = 0;
			if (men_sel < 0)
				men_sel = men_l - 1;
			
			vec_n back_size = true_size(lines[1 + men_sel].getLocalBounds());
			RectangleShape back_rec;
			back_rec.setSize(back_size);

			back_rec.setFillColor(Color::White);
			back_rec.setPosition(lines[1 + men_sel].getPosition());
			lines[1 + men_sel].setFillColor(Color::Black);

			if (input::keyboard.wasPressed(input::key_state::keys::Return))
				input::flush_back::men_cmd = men_sel + 1;

			window2.draw(back_rec);
			for (int i = 0; i < lines.size(); i++)
				window2.draw(lines[i]);
		}

		if (game_state == 1)
		{
			double time_prc = ws::health / ws::health_max;
			double health_prc = ws::health / ws::health_max;
			double thrust_prc = ws::player_thrust / ws::player_thrust_max;

			draw_crep(window2, crep_dists, crep_items, { time_prc, health_prc, thrust_prc });

			int target = ws::target;
			if (target >= 0 && target < ws::names.size())
			{
				string name = ws::names[target];
				double vel = 0;
				double dist = 0;
				double min_dist = 0;
				double min_time = 0;
				if (ws::target_close.size())
					min_dist = phys::vmag(ws::target_close[0] - ws::player_close[0]);
				if (ws::target_time.size())
					min_time = ws::target_time[0];

				short box_mode = 0;
				dist = int(phys::vmag(ws::bodies.back() - ws::bodies[target]));
				if (ws::cepting and ws::target == ws::intercept)
				{
					box_mode = 1;
					min_time = ws::cept_time;
					min_dist = 0;
				}
				else
				{
					box_mode = 0;
					min_time = ws::target_time[0] - ws::world_time;
					min_dist = int(ws::target_min[0]);
					if (ws::cepting and ws::cept_time < ws::target_time[0])
					{
						min_time = NAN;
						min_dist = NAN;
					}
				}
				draw_plnbox(window2, box_mode, vec_n(-5, -5), name, vel, dist, min_dist, min_time);


				if (ws::cepting and ws::target != ws::intercept)
				{
					string name = ws::names[ws::intercept];
					double vel = 0;
					double dist = int(phys::vmag(ws::bodies.back() - ws::bodies[ws::intercept]));
					double min_dist = 0;
					double min_time = ws::cept_time - ws::world_time;
					short box_mode = 2;

					draw_plnbox(window2, box_mode, vec_n(-210, 5), name, vel, dist, min_dist, min_time);
				}
			}

			draw_kesslers(window2, ws::kesses, ws::sizes_kess.back(), ws::bodies[ws::focus], ws::zoom, vec_n(window2.getSize().x / 2, window2.getSize().y / 2));
		}

		if (game_state == 2)

		{
			vector<Text> lines;
			vec_n h_pos = w_size;

			{
				Text line;
				line.setFont(font);
				line.setString(pause_head);
				line.setCharacterSize(60);

				h_pos.x *= 0.1;
				h_pos.y *= 0.3;
				line.setPosition(h_pos);

				lines.push_back(line);
			}


			for (int i = 0; i < pause_l; i++)
			{
				Text line;
				line.setFont(font);
				line.setString(pause_items[i]);
				line.setCharacterSize(30);

				vec_n l_pos = h_pos;
				l_pos.y += w_size.y * 0.1 * (i + 1);
				line.setPosition(l_pos);

				lines.push_back(line);
			}

			if (input::keyboard.wasPressed(input::key_state::keys::Up))
				pause_sel--;
			if (input::keyboard.wasPressed(input::key_state::keys::Down))
				pause_sel++;

			if (pause_sel >= pause_l)
				pause_sel = 0;
			if (pause_sel < 0)
				pause_sel = pause_l - 1;

			vec_n back_size = true_size(lines[1 + pause_sel].getLocalBounds());
			RectangleShape back_rec;
			back_rec.setSize(back_size);

			back_rec.setFillColor(Color::White);
			back_rec.setPosition(lines[1 + pause_sel].getPosition());
			lines[1 + pause_sel].setFillColor(Color::Black);

			if (input::keyboard.wasPressed(input::key_state::keys::Return))
				input::flush_back::pause_cmd = pause_sel + 1;

			window2.draw(back_rec);
			for (int i = 0; i < lines.size(); i++)
				window2.draw(lines[i]);
		}
	}

	void ui_init()
	{
		font.loadFromFile("courier.ttf");
		{
			men_head = "Dear Isaac,";

			men_items.push_back("play,");
			men_items.push_back("options,");
			men_items.push_back("credits,");
			men_items.push_back("quit.");
			men_l = men_items.size();
		}
		{
			pause_head = "Paused,";
			pause_items.push_back("resume,");
			pause_items.push_back("controls,");
			pause_items.push_back("options,");
			pause_items.push_back("quit.");
		}
		{
			crep_items.push_back("Time:......");
			crep_items.push_back("Health:....");
			crep_items.push_back("Thrust:....");
		}
		{
			crep_dists[0] = 5;	//Margin
			crep_dists[1] = 15;	//Char size
			crep_dists[2] = 20;	//Line Height
			crep_dists[3] = 100;//Line Width
			crep_dists[4] = 200;//Box  Width
		}
	}
}