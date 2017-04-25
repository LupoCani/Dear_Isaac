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

	vector<string> crep_items;
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

			if (men_sel >= men_l)
				men_sel = 0;
			if (men_sel < 0)
				men_sel = men_l - 1;
			
			sf::FloatRect back_size = lines[1 + men_sel].getLocalBounds();
			RectangleShape back_rec;
			back_rec.setSize(vec_n(back_size.left + back_size.width, back_size.top + back_size.height));

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

				title.setFont(font);
				title.setCharacterSize(crep_char);
				title.setString(crep_items[i]);

				box.setSize(vec_n(crep_end_x * phys::rand_part(), crep_y * 0.9));
				line.setSize(vec_n(crep_end_x, crep_y * 0.2));
				end.setSize(vec_n(2, crep_y * 0.9));

				sf::FloatRect text_size = title.getLocalBounds();
				RectangleShape back_rec;
				title.setOrigin(vec_n(0, (text_size.top + text_size.height)/2));

				box.setOrigin (vec_n(0, box.getSize().y  / 2.0));
				line.setOrigin(vec_n(0, line.getSize().y / 2.0));
				end.setOrigin (vec_n(0, end.getSize().y  / 2.0));

				title.setPosition(vec_n(x, y + crep_y / 2));
				box.setPosition(vec_n(x_col, y + crep_y / 2));
				line.setPosition(vec_n(x_col, y + crep_y / 2));
				end.setPosition(vec_n(x_end, y + crep_y / 2));
				
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

		crep_items.push_back("Time:......");
		crep_items.push_back("Health:....");
		crep_items.push_back("Thrust:....");

		font.loadFromFile("courier.ttf");
		crep_char = 15;
		crep_x = 5;
		crep_y = 20;

		crep_col_x = 100;
		crep_end_x = 200;
	}
}