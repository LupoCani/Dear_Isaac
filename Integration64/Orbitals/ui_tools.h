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

	float crep_char;
	float crep_x;
	float crep_col_x;
	float crep_end_x;
	float crep_y;

	void hp(double hp, double hp_max)
	{
	}

	void run_ui()
	{
		using shared::game_state;
		Vector2u w_size_sfu = window2.getSize();
		vec_n w_size(w_size_sfu.x, w_size_sfu.y);
		input::flush_back::men_cmd = 0;

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
			if (men_sel > men_l)
				men_sel = 0;
			
			sf::FloatRect back_size = lines[1 + men_sel].getLocalBounds();
			RectangleShape back_rec(vec_n(back_size.width, back_size.height));
			back_rec.setFillColor(Color::White);
			back_rec.setPosition(lines[1 + men_sel].getPosition());

			if (input::keyboard.wasPressed(input::key_state::keys::Return))
				input::flush_back::men_cmd = men_sel + 1;

			window2.draw(back_rec);
			for (int i = 0; i < lines.size(); i++)
				window2.draw(lines[i]);
		}

		if (game_state == 1)
		{
			for (int i = 0; i < 3; i++)
			{
				Text title;
				RectangleShape box;
				RectangleShape line;
				RectangleShape end;
				float x = crep_x;
				float x_col = x + crep_col_x;
				float x_end = x_col + crep_end_x;
				float y = crep_x + crep_y * i;

				box.setPosition (vec_n(x, y));
				line.setPosition(vec_n(x_col, y + crep_y * 0.4 ));
				end.setPosition (vec_n(x_end, y));

				box.setSize(vec_n(crep_end_x, crep_y * 0.9));
				line.setSize(vec_n(crep_end_x, crep_y * 0.2));
				end.setSize(vec_n(2, crep_y * 0.9));

				string os = "";
				for (int i = 0; i < 10; i++)
				{
					os += char(phys::rand_part() * 200);
				}
				title.setFont(font);
				title.setString(os);
				
				window2.draw(title);
				window2.draw(box);
				window2.draw(line);
				window2.draw(end);
			}
		}
	}

	void ui_init()
	{
		men_head = "Dear Isaac,";

		men_items.push_back("play,");
		men_items.push_back("options,");
		men_items.push_back("credits,");
		men_items.push_back("quit.");
		men_l = men_items.size();

		font.loadFromFile("courier.ttf");
		crep_char = 50;
		crep_x = 50;
		crep_y = 50;

		crep_col_x = 500;
		crep_end_x = 500;
	}
}