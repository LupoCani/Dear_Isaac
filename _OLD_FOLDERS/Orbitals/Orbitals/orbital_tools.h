#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <windows.h>
#define skip if (false)


sf::RenderWindow window2(sf::VideoMode(1080, 860), "Orbitals");

namespace phys{
	using namespace std;

	double w_time = 0;

	clock_t cur_time = clock();
	vector <sf::RectangleShape> tail;
	int cps = CLOCKS_PER_SEC;
	const long double M_PI = acos(-1.0L);
	const long double M_2PI = M_PI * 2L;
	const long double M_PI2 = M_PI / 2L;
	const long double M_3PI2 = M_PI2 * 3L;

	const long double M_E = exp(1.0);

	struct vec_n
	{
		long double x = 0;
		long double y = 0;
		long double z = 0;

		operator sf::Vector2f() {
			return sf::Vector2f(x, y);
		}
	};
	vector <vec_n> tail_coord;


	sf::Color color_enum(int in)
	{
		switch (in)
		{
		case 0:
			return sf::Color::White;
		case 1:
			return sf::Color::Blue;
		case 2:
			return sf::Color::Yellow;
		case 3:
			return sf::Color::Green;
		case 4:
			return sf::Color::Cyan;
		case 5:
			return sf::Color::Magenta;
		default:
			return sf::Color::Red;
		}
	}

	void draw_circle(double x, double y, int beblue, double zoom) //Render either of the planets
	{

		sf::CircleShape planet;

		zoom *= 10;
		if (beblue > 5)
		{
			zoom *= 0.01;
		}
		if (zoom < 0.5)
			zoom = 0.5;

		planet.setPosition(x, y); // Subtract five to compensate for circle size
		planet.setOrigin(3 * zoom, 3 * zoom);
		planet.setFillColor(color_enum(beblue)); // Blue for orbiting body, red for central
		planet.setRadius(6 * zoom);


		window2.draw(planet);
	}

	vec_n handle_scale_single(vec_n list, vec_n origo, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		origo.x *= scale;
		origo.y *= scale;

		list.x *= scale;
		list.y *= scale;

		list.x = 500 + (list.x - origo.x);
		list.y = 500 - (list.y - origo.y);

		return list;
	}

	void draw_tail(double x, double y, int color, bool draw = false, bool add = true, double zoom = 1, vec_n origo = {}, vec_n parent = {}) //Render either of the tails
	{
		sf::RectangleShape holder; //Declare tail particle
		vec_n holder_vec;

		if (add)
		{
			holder_vec.x = x;
			holder_vec.y = y;
			holder.setPosition(x, y); //Apply properties
			holder.setSize(sf::Vector2f(2.6, 2.6)); //Radius
			holder.setOrigin(1.2, 1.3);
			holder.setFillColor(!color ? sf::Color::Cyan : sf::Color::Red); //Green for orbiting body, cyan for central

			tail.push_back(holder); //Store in main vector
			tail_coord.push_back(holder_vec);
		}

		/*
		if (draw) for (int i = 0; i < tail.size(); i++)
		{
		holder = tail[i];
		holder.setPosition(handle_scale_single(tail_coord[i], origo, zoom));
		window2.draw(holder); //Draw all tail particles
		}
		*/
		if (draw) for (int i = 0; i < tail.size(); i++)
		{
			holder = tail[i];
			vec_n pos = tail_coord[i];

			if (color)
			{
				pos.x += 10000;// parent.x;
				pos.y += 10000;// parent.y;
			}

			holder.setPosition(handle_scale_single(pos, origo, zoom));
			window2.draw(holder); //Draw all tail particles
		}
	}

	struct body
	{
		double ax_a;	//Semimajor axis
		double ax_b;	//Semiminor axis
		double ecc;		//Orbital eccentrictity
		double epoch;	//Timestamp of any passage through the periapsis
		double t_p;		//For elliptical orbits, the orbital period. Undefined for hyperbolic orbits.
		double t_l;		//Timestamp on the current values on time-dependent paramters
		double u;		//Gravitational parameter
		double normal;	//angle of the major axis
		double SOI;		//radius of the sphere of influence
		short shape;	//shape of the orbit. 0 denotes elliptical orbits, 1 denotes hyperbolic orbits
		bool inverse;	//true if the orbit runs counter-clockwise
		double expiry;	//The first point in time beyond which the orbit may enter another SOI

		vec_n vel;		//Velocity at a given time
		vec_n pos;		//Position at a given time

		double En;		//Orbital energy
		double Ar;		//Area swept per unit of time
		double Mn;		//Mean anomaly swept per unit of time

		bool isPlayer = false; //true if the body is the player
		bool isSun = false;	//true if the body is the sun

		body* parent;		//Pointer to parent body
		body* self = this;	//Pointer to self
		string name;		//Name of body
	};

	body sun, moon, planet, pluto, dune, yavin, plyr;
	vector<body*> bodies;



	struct world_state
	{
		vector<vec_n> bodies;
		vector<string> names;
		vector<vector<vec_n>> lines;

		double rotation;
		double throttle;
	};


	double clamp(double in)
	{
		return in > 1 ? 1 : (in < -1 ? -1 : in);
	}

	double sqr(double in)
	{
		return in*in;
	}

	double sqrt_a(double in)
	{
		return sqrt(abs(in));
	}

	vec_n vec_to_pos(double v, double r)
	{
		vec_n out;
		out.x = cos(v) * r;
		out.y = sin(v) * r;
		return out;
	}

	vec_n t_to_screen(vec_n in, double zoom = 1)
	{
		zoom = 100;
		vec_n out;
		out.y = 500 - in.y / zoom;
		out.x = 500 + in.x / zoom;

		return out;
	}

	double ang_wrap(double in, int quad = 4)
	{
		double ceil = double(quad) * M_PI2;
		double floor = ceil - M_2PI;

		if (in >= ceil) do
		{
			in -= M_2PI;
		} while (in >= ceil);
		else if (in < floor) do
		{
			in += M_2PI;
		} while (in < floor);
		return in;
	}

	vec_n get_rel(vec_n point, vec_n orig)
	{
		vec_n out;
		out.x = point.x - orig.x;
		out.y = point.y - orig.y;
		return out;
	}

	double dist(double x, double y)
	{
		return sqrt(sqr(x) + sqr(y));
	}

	double dist_v(vec_n in)
	{
		return dist(in.x, in.y);
	}

	double get_r(double ang, body sat)
	{
		return sat.ax_a * (1 - sqr(sat.ecc)) / (1 + sat.ecc * cos(ang));
	}

	double to_rad(double ang, double scale)
	{
		return M_2PI * ang / scale;
	}

	double get_V_r(double r, body sat)
	{
		return acos(clamp(sat.ax_a * (1 - sqr(sat.ecc)) / (sat.ecc * r) - 1 / sat.ecc));
	}

	void fill_tail() //Generate one thousand tail particles, fill main vector.
	{
		for (int i = 0; i < 10; i++)
		{
			sf::RectangleShape holder;
			holder.setPosition(0, 0);
			holder.setSize(sf::Vector2f(2.6, 2.6));
			holder.setFillColor(sf::Color::Green);
			tail.push_back(holder);
		}
	}



	struct dual_val
	{
		long double one = 0;
		long double two = 0;
	};

	double switch_co(double in)
	{
		return sqrt(1 - pow(in, 2));
	}

	double get_sinE(double tV, double ecc)
	{
		return (sqrt(1 - pow(ecc, 2)) * sin(tV)) / (1 + ecc * cos(tV));
	}
	double get_cosE(double tV, double ecc)
	{
		double cosV = cos(tV);
		return (ecc + cosV) / (1 + ecc * cosV);
	}
	double get_E(double tV, double ecc)
	{
		if (tV < M_PI)
			return acos(get_cosE(tV, ecc));
		else
			return M_2PI - acos(get_cosE(tV, ecc));
	}

	double get_V(double E, double ecc)
	{
		double cosE = cos(E);
		double V = acos((cosE - ecc) / (1 - ecc * cosE));

		if (E < M_PI)
			return V;
		else
			return M_2PI - V;
	}

	double get_M_H(double V, double ecc)
	{
		V = ang_wrap(V, 2);

		double cos_V = cos(V);
		double cosh_hE = (ecc + cos_V) / (1 + ecc * cos_V);
		double hE = acosh(cosh_hE);
		if (V < 0)
			hE = -hE;

		double sinh_hE = sinh(hE);

		return ecc * sinh_hE - hE;
	}


	double get_M(double tV, double ecc)
	{
		if (ecc > 1) return get_M_H(tV, ecc);
		double cosE = get_cosE(tV, ecc);
		double sinE = switch_co(cosE);

		if (tV < M_PI)
			return acos(cosE) - ecc * sinE;
		else
			return M_2PI - acos(cosE) - ecc * -sinE;
	}

	double get_M_E(double E, double ecc)
	{
		return E - ecc * sin(E);
	}

	double V_to_rinc(double V, body sat)
	{
		V = M_2PI - V;
		double E = get_E(V, sat.ecc);

		double rinc = V + atan(tan(M_PI2 - E) * sat.ax_b / sat.ax_a);
		if (V > M_PI) rinc += -M_PI;

		return ang_wrap(rinc - M_PI2, 2);
	}

	double rinc_to_V(double rinc, body sat, int precision, bool mid_range)
	{
		rinc = ang_wrap(rinc, 2);
		double ceil = M_2PI;
		double floor = 0;
		double mid = 0;

		if (mid_range)
		{
			ceil = get_V(M_3PI2, sat.ecc);
			floor = get_V(M_PI2, sat.ecc);
		}
		else
		{
			if (rinc > M_PI2)
				floor = get_V(M_3PI2, sat.ecc);
			else
				ceil = get_V(M_PI2, sat.ecc);
		}

		for (int i = 0; i < precision; i++)
		{
			mid = (ceil + floor) / 2;

			if (mid_range == V_to_rinc(mid, sat) > rinc)
				ceil = mid;
			else
				floor = mid;
		}

		return (ceil + floor) / 2;
	}

	vec_n rot_vec(vec_n in, double rot)
	{
		double l = dist_v(in);
		double angle = atan2(in.y, in.x);

		angle += rot;

		in.y = l * sin(angle);
		in.x = l * cos(angle);

		return in;
	}

	double get_esc_vel(body parent, double r)
	{
		return sqrt(2 * parent.u / r);
	}

	double do_orbit_precise_H(double part, body sat, int precision)
	{
		double ceil = acos(-1 / sat.ecc);
		double floor = -ceil;
		double mid = 0;

		for (int i = 0; i < precision; i++)
		{
			if (get_M(mid, sat.ecc) > part)
				ceil = mid;
			else
				floor = mid;

			mid = (ceil + floor) / 2;
		}

		return (ceil + floor) / 2;
	}

	double do_orbit_precise(double part, body sat, int precision)
	{
		if (sat.shape) return do_orbit_precise_H(part, sat, precision);

		double ceil = M_2PI;
		double floor = 0;
		double mid = M_PI;

		for (int i = 0; i < precision; i++)
		{
			if (get_M(mid, sat.ecc) > part)
				ceil = mid;
			else
				floor = mid;

			mid = (ceil + floor) / 2;
		}

		return mid;
	}

	void do_orbit_phys(body &sat, double d_time, double px, double py)
	{
		double m_time = d_time - sat.t_l;
		sat.t_l = d_time;

		body &parent = *sat.parent;
		vec_n pos_rel;
		vec_n force_vec;
		pos_rel = get_rel(sat.pos, parent.pos);

		double angle = atan2(pos_rel.y, pos_rel.x);
		double pos_mag = dist_v(pos_rel);
		double F = -parent.u / sqr(pos_mag);

		force_vec.x = F * cos(angle);
		force_vec.y = F * sin(angle);

		sat.vel.x += force_vec.x * m_time *px;
		sat.vel.y += force_vec.y * m_time *py;

		sat.pos.x += sat.vel.x * m_time;
		sat.pos.y += sat.vel.y * m_time;
	}


	double time_to_M(body sat, double w_time)
	{
		double d_time = w_time - sat.epoch;

		double M = d_time * sat.Mn;
		if (sat.shape == 0)
			return ang_wrap(M);

		return M;
	}

	double M_to_time(body sat, double M, double w_time_min = -1)
	{
		double d_time = M / sat.Mn;
		double w_time = d_time + sat.epoch;

		if (sat.shape == 0 && w_time_min >= 0)
		{
			double r_time = M_2PI / sat.Mn;
			while (w_time < w_time_min)
				w_time += r_time;
		}

		return w_time;
	}

	void do_orbit(body &sat, double w_time)
	{
		sat.t_l = w_time;
		double d_time = w_time - sat.epoch;

		body &parent = *sat.parent;

		double part = d_time * sat.Mn;
		double V;
		double radius;

		if (sat.shape == 0)
			part = ang_wrap(part);

		V = do_orbit_precise(part, sat, 100);
		if (sat.inverse)
			V = M_2PI - V;

		radius = get_r(V, sat);

		sat.pos = vec_to_pos(-V, radius);
		sat.pos = rot_vec(sat.pos, sat.normal);

		sat.pos.x += parent.pos.x;
		sat.pos.y += parent.pos.y;

		double vel_mag = sqrt(2 * (sat.En + parent.u / radius));
		double vel_ang = acos(sat.Ar / (radius * vel_mag));
		if (V < M_PI)
			vel_ang = -vel_ang;

		if (sat.inverse)
			vel_mag = -vel_mag;

		vel_ang = ang_wrap(-V - vel_ang + sat.normal + M_PI2);

		sat.vel.x = -vel_mag * cos(vel_ang);
		sat.vel.y = -vel_mag * sin(vel_ang);

		sat.vel.x += parent.vel.x;
		sat.vel.y += parent.vel.y;
	}

	double Mn_from_V(double ecc, double r, double v)
	{
		double dV = sqrt(1 - sqr(ecc)) * (1 - ecc) / (1 + ecc);
		double Mn = dV * v / r;

		return Mn;
	}

	double Mn_from_r(double ax_a, double u)
	{
		return sqrt_a(u / pow(ax_a, 3));
	}

	vector<dual_val> get_extremes(double vel, double ang, double r, double u)
	{
		dual_val out_v, out_r;
		vector<dual_val> out(2);

		double part_div = u / (cos(ang) * vel * r);
		double part_root = sqrt(pow(-part_div, 2) - 2 * u / r + pow(vel, 2));
		double part_area = cos(ang) * vel * r;

		out_v.one = part_div + part_root;
		out_v.two = part_div - part_root;

		out_r.one = part_area / out_v.one;
		out_r.two = part_area / out_v.two;

		out[0] = out_r;
		out[1] = out_v;

		return out;
	}


	void get_expiry(body sat_in, vector<body*> list_in, double w_time)
	{
		body &sat = *(sat_in.self);
		vector<body*> list(0);
		double detail = 1000;

		for (int i = 0; i < list_in.size(); i++)
		{
			body sat_2 = (*(list_in[i]));
			if (sat_2.parent == sat.parent && sat_2.self != sat.self)
				list.push_back(sat.parent);
		}

		double w_time_end = HUGE_VAL;
		double dist_end = HUGE_VAL;
		bool expire_true = false;

		for (int i = 0; i < list.size(); i++)
		{
			body sat_2 = (*(list[i]));
			double SOI = sat_2.SOI;
			double dist_min = HUGE_VAL;
			double w_time_close;
			bool expire = false;

			for (int i2 = 0; i2 < detail; i2++)
			{
				double V = to_rad(i2, detail);
				double r = get_r(V, sat);
				vec_n pos = vec_to_pos(sat.normal - V, r);

				double M = get_M(V, sat.ecc);
				double w_time_pot = M_to_time(sat, M, w_time);
				double M_2 = time_to_M(sat_2, w_time_pot);

				double V_2 = do_orbit_precise(M_2, sat_2, log2(detail));
				double r_2 = get_r(V_2, sat_2);
				vec_n pos_2 = vec_to_pos(sat_2.normal - V_2, r_2);

				double dist = dist_v(get_rel(pos, pos_2));

				if (expire_true && w_time_pot > w_time_end)
					break;

				if (dist < dist_min)
				{
					dist_min = dist;
					w_time_close = w_time_pot;
				}
				if (dist < SOI)
				{
					expire = true;
					break;
				}
			}

			if (expire)
			{
				w_time_end = w_time_close;
				dist_end = dist_min;
				expire_true = true;
			}
		}

	}

	body* get_parent(body sat, vector<body*> list)
	{
		body *parent_presumed = sat.parent;
		if (!sat.isPlayer)
			throw;

		for (vector<body*>::iterator i = list.begin(); i != list.end(); i++)
		{
			if (*i == sat.self)
				list.erase(i);

			if ((**i).isSun)
				parent_presumed = *i;
		}

		for (int i = 0; i < list.size(); i++)
		{
			body &parent_pot = *(list[i]);
			body &parent_cur = *parent_presumed;

			double dist = dist_v(get_rel(sat.pos, parent_pot.pos));

			if (dist < parent_pot.SOI)
				if (parent_pot.SOI < parent_cur.SOI)
					parent_presumed = list[i];
		}

		return parent_presumed;
	}

	void make_orbit(body &sat_in, double w_time)
	{
		body &sat = *sat_in.self;
		body &parent = *(sat.parent);

		vec_n pos, vel;
		double pos_mag, pos_ang, vel_mag, vel_ang, rel_ang, u;

		pos = get_rel(sat.pos, parent.pos);
		vel = get_rel(sat.vel, parent.vel);
		sat.inverse = false;
		u = parent.u;

		bool repeat;
		do
		{
			repeat = false;
			pos_mag = dist_v(pos);
			pos_ang = atan2(pos.y, pos.x);

			vel_mag = dist_v(vel);
			vel_ang = atan2(vel.y, vel.x);

			rel_ang = M_PI2 + vel_ang - pos_ang;
			rel_ang = ang_wrap(rel_ang, 2);

			if (abs(rel_ang) > M_PI2)
			{
				vel.x = -vel.x;
				vel.y = -vel.y;
				sat.inverse = true;
				repeat = true;
			}
			if (vel_mag == get_esc_vel(parent, pos_mag))
			{
				vel.x = vel.x * 1.00001;
				vel.y = vel.y * 1.00001;
				repeat = true;
			}

		} while (repeat);

		vector<dual_val> r_vec_vector = get_extremes(vel_mag, rel_ang, pos_mag, u);

		dual_val r_vec = r_vec_vector[0];
		dual_val r_vec_vel = r_vec_vector[1];

		sat.ax_a = (r_vec.two + r_vec.one) / 2;
		sat.ecc = (r_vec.two - r_vec.one) / (r_vec.two + r_vec.one);
		sat.ax_b = sqrt_a(1 - sqr(sat.ecc)) * sat.ax_a;
		sat.shape = sat.ecc > 1;
		sat.SOI = sat.ax_a * pow(sat.u / u, 0.4);


		double tV = get_V_r(pos_mag, sat);
		if (rel_ang < 0)
			tV = ang_wrap(-tV);

		sat.normal = pos_ang + tV;

		if (sat.inverse)
			tV = ang_wrap(-tV);
		sat.En = sqr(vel_mag) / 2 - u / pos_mag;			//Store orbital Energy, velocital energy plus gravitational energy
		sat.Ar = vel_mag * pos_mag * cos(rel_ang);			//Store area swept per unit of time, altitude times velocity times cosine-of-the-angle
		sat.Mn = Mn_from_r(sat.ax_a, u);					//Store mean angular motion of the body, as per yet another equation dug up from obscure wikipedia stubs
		sat.epoch = w_time - get_M(tV, sat.ecc) / sat.Mn;

		if (sat.shape || r_vec.two > parent.SOI)
		{
			double V_max = get_V_r(parent.SOI, sat);
			double M_max = get_M(V_max, sat.ecc);
			sat.expiry = M_max / sat.Mn + sat.epoch;
		}
	}

	struct keys{
		double y = 0;
		double x = 0;
		double z = 0;
		double lc = 0;
		double pause;


		vector<int> pressed;
		struct list : sf::Keyboard{
			using Keyboard::Key;
		};
	};

	keys render_cst() {

		sf::Event input;
		keys out;

		while (window2.pollEvent(input)) {

			if ((input.type == sf::Event::KeyPressed) && (input.key.code == keys::list::Up)) {
				out.y += 1;
			}
			if ((input.type == sf::Event::KeyPressed) && (input.key.code == sf::Keyboard::Down)) {
				out.y += -1;
			}
			if ((input.type == sf::Event::KeyPressed) && (input.key.code == sf::Keyboard::Right)) {
				out.x += 1;
			}
			if ((input.type == sf::Event::KeyPressed) && (input.key.code == sf::Keyboard::Left)) {
				out.x += -1;
			}
			if (input.type == sf::Event::MouseWheelMoved) {
				out.z = input.mouseWheel.delta;
			}
			if ((input.type == sf::Event::KeyPressed) && (input.key.code == sf::Keyboard::RControl)) {
				out.lc = 1;
			}
			if ((input.type == sf::Event::KeyPressed) && (input.key.code == sf::Keyboard::RShift)) {
				out.pause = 1;
			}
		}
		return out;
	}


	vector<vec_n> update_all_and_convert(vector<body*> bodies, double w_time, bool e_mode, double vel_x, double vel_y, double fact)
	{
		vector<vec_n> out(1);
		out[0] = (*bodies[0]).pos;
		bool emodelast = false;


		out[0].x = 500;
		out[0].y = 500;
		for (int i = 1; i < bodies.size(); i++)
		{
			body &sat = *((*bodies[i]).self);

			if (sat.isPlayer)
			{
				body* parent_new = get_parent(sat, bodies);
				if (parent_new != sat.parent)
				{
					sat.parent = parent_new;
					make_orbit(sat, w_time);
				}
			}

			if (!sat.isPlayer || !e_mode)
			{
				do_orbit(sat, w_time);
			}
			else
			{
				do_orbit_phys(sat, w_time, vel_x * fact, vel_y*fact);
				make_orbit(sat, w_time);
				//do_orbit(newsat, w_time);
			}
			//change_parent(sat, bodies);
			emodelast = e_mode;
			vec_n pos_temp = sat.pos;
			//pos_temp = t_to_screen(pos_temp);
			out.push_back(pos_temp);
		}

		//get_expiry(*bodies[bodies.size() - 1], bodies, w_time);
		return out;
	}

	vector<vec_n> positions;

	void run_engine()
	{
		positions = update_all_and_convert(bodies, shared::r_time, 0, 0, 0, 0);
	}


	void phys_init()
	{
		sun.pos.x = 0;
		sun.pos.y = 0;
		sun.vel.x = 0;
		sun.vel.y = 0;
		sun.u = pow(10, 14);
		sun.SOI = HUGE_VAL;
		sun.isSun = true;
		sun.name = "Sun";

		moon.pos.x = -25100;
		moon.pos.y = 0;
		moon.vel.x = 0;
		moon.vel.y = 40000;
		moon.u = pow(10, 8);
		moon.parent = &sun;
		moon.name = "moon";

		planet.pos.x = -25000;
		planet.pos.y = 0;
		planet.vel.x = 0;
		planet.vel.y = 39900;
		planet.u = pow(10, 10);
		planet.parent = &sun;
		planet.name = "planet";

		pluto.pos.x = -30000;
		pluto.pos.y = 0;
		pluto.vel.x = 0;
		pluto.vel.y = 40000;
		pluto.u = pow(10, 10);
		pluto.parent = &sun;
		pluto.name = "pluto";

		dune.pos.x = -35000;
		dune.pos.y = 0;
		dune.vel.x = 0;
		dune.vel.y = 40000;
		dune.u = pow(10, 10);
		dune.parent = &sun;
		dune.name = "dune";

		yavin.pos.x = -40000;
		yavin.pos.y = 0;
		yavin.vel.x = 0;
		yavin.vel.y = 40000;
		yavin.u = pow(10, 10);
		yavin.parent = &sun;
		yavin.name = "yavin";

		plyr.pos.x = -40200;
		plyr.pos.y = 0;
		plyr.vel.x = 0;
		plyr.vel.y = 43000;
		plyr.u = pow(10, 1);
		plyr.parent = &sun;
		plyr.name = "plyr";

		bodies.push_back(&sun);		//White
		bodies.push_back(&moon);	//Blue
		bodies.push_back(&planet);	//Yellow
		bodies.push_back(&pluto);	//Green
		bodies.push_back(&dune);	//Cyan
		bodies.push_back(&yavin);	//Magenta
		bodies.push_back(&plyr);	//Red

		make_orbit(moon, w_time);
		make_orbit(planet, w_time);
		make_orbit(pluto, w_time);
		make_orbit(dune, w_time);
		make_orbit(yavin, w_time);
		make_orbit(plyr, w_time);
	}
}

