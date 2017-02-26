#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <windows.h>
#define skip if (false)
#define ORBITAL_TOOLS_LOADED true


namespace shared 					//Declare basic shared parameters
{	
	sf::RenderWindow window2;

	clock_t r_time;
	clock_t l_time;
	clock_t s_time;
	const double cps = CLOCKS_PER_SEC;

	struct world_state;
}

namespace input					//Declare the input system. In a namespace becuse putting sf:: in phys feels unclean
{
	using namespace sf;

	struct key_state {

		short scroll = 0;
		std::vector<int> pressed;
		struct keys : sf::Keyboard {
			using Keyboard::Key;
		};

		bool isPressed(int key)
		{
			for (int i = 0; i < pressed.size(); i++)
				if (key == pressed[i])
					return true;

			return false;
		}
	};

	key_state keyboard;

	void run_input()
	{
		Event input;
		keyboard.scroll = 0;

		while (shared::window2.pollEvent(input)) {

			if (input.type == Event::KeyPressed)
			{
				bool isAdded = false;

				for (int i = 0; i < keyboard.pressed.size(); i++)
					if (input.key.code == keyboard.pressed[i])
						isAdded = true;

				if (!isAdded)
					keyboard.pressed.push_back(input.key.code);
			}

			if (input.type == Event::KeyReleased)
				for (int i = 0; i < keyboard.pressed.size(); i++)
					if (input.key.code == keyboard.pressed[i])
						keyboard.pressed.erase(keyboard.pressed.begin() + i);

			if (input.type == Event::MouseWheelMoved)
				keyboard.scroll = input.mouseWheel.delta;
		}
	}
}

namespace phys					//Declare various classes and functions
{
	struct vec_n;

	double clamp(double in);
	double sqr(double in);
	double sqrt_a(double in);
	vec_n vec_to_pos(double v, double r);
	double atan2(vec_n pos);
	vec_n t_to_screen(vec_n in, double zoom = 1);
	double ang_wrap(double in, int quad = 4);
	double vmag(double x, double y);
	double vmag(vec_n in);
	vec_n rot_vec(vec_n in, double rot);
	double to_rad(double ang, double scale);
	using std::atan2;

	struct vec_r
	{
		double ang = 0;
		double mag = 0;
	};

	struct vec_n
	{
		double x = 0;
		double y = 0;

		operator sf::Vector2f() {
			return sf::Vector2f(x, y);
		}

		operator phys::vec_r() {
			phys::vec_r lhs;
			lhs.ang = phys::atan2(*this);
			lhs.mag = phys::vmag(*this);
			return lhs;
		}

		vec_n& operator +=(const vec_n& rhs)
		{
			(*this).x = (*this).x + rhs.x;
			(*this).y = (*this).y + rhs.y;
			return *this;
		}

		vec_n& operator-=(const vec_n& rhs)
		{
			(*this).x = (*this).x - rhs.x;
			(*this).y = (*this).y - rhs.y;
			return *this;
		}
		vec_n& operator*=(const vec_n& rhs)
		{
			(*this).x = (*this).x * rhs.x;
			(*this).y = (*this).y * rhs.y;
			return *this;
		}
		vec_n& operator*=(const double& rhs)
		{
			(*this).x = (*this).x * rhs;
			(*this).y = (*this).y * rhs;
			return *this;
		}
	};

	inline vec_n operator+(vec_n lhs, const vec_n& rhs)
	{
		lhs += rhs;
		return lhs;
	}

	inline vec_n operator-(vec_n lhs, const vec_n& rhs)
	{
		lhs -= rhs;
		return lhs;
	}

	inline vec_n operator*(vec_n lhs, const vec_n& rhs)
	{
		lhs *= rhs;
		return lhs;
	}
	inline vec_n operator*(vec_n lhs, const double& rhs)
	{
		lhs *= rhs;
		return lhs;
	}

}

namespace shared				//Declare the ever-so-important vec_n. In shared because render_tools needs to use it.
{
	using phys::vec_n;

	struct world_state
	{
		std::vector<vec_n> bodies;
		std::vector<std::string> names;
		std::vector<std::vector<vec_n>> paths;

		double zoom = 1;
		int focus = 0;

		double player_rotation; //Angle in radians, sweeping counter-clockwise.
	};

	world_state screen_state;
}

namespace phys
{
	using std::vector;

	double w_time = 0;

	clock_t cur_time = clock();
	vector <sf::RectangleShape> tail;
	int cps = CLOCKS_PER_SEC;
	const long double M_PI = acos(-1.0L);
	const long double M_2PI = M_PI * 2L;
	const long double M_PI2 = M_PI / 2L;
	const long double M_3PI2 = M_PI2 * 3L;

	const long double M_E = exp(1.0);

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


		shared::window2.draw(planet);
	}

	vec_n handle_scale_single(vec_n list, vec_n origo, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		origo *= scale;

		list *= scale;

		list -= origo;

		list.x = 500 + list.x;
		list.y = 500 - list.y;

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
			shared::window2.draw(holder); //Draw all tail particles
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
		std::string name;		//Name of body
	};

	body sun, moon, planet, pluto, dune, yavin, plyr;
	vector<body*> bodies;


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

	double atan2(vec_n pos)
	{
		return std::atan2(pos.y, pos.x);
	}

	vec_n t_to_screen(vec_n in, double zoom)
	{
		zoom = 100;
		vec_n out;
		out.y = 500 - in.y / zoom;
		out.x = 500 + in.x / zoom;

		return out;
	}

	double ang_wrap(double in, int quad)
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

	double vmag(double x, double y)
	{
		return sqrt(sqr(x) + sqr(y));
	}

	double vmag(vec_n in)
	{
		return vmag(in.x, in.y);
	}

	vec_n rot_vec(vec_n in, double rot)
	{
		double l = vmag(in);
		double angle = atan2(in);

		angle += rot;

		in.y = l * sin(angle);
		in.x = l * cos(angle);

		return in;
	}

	double to_rad(double ang, double scale)
	{
		return M_2PI * ang / scale;
	}

	double get_r(double ang, body sat)
	{
		return sat.ax_a * (1 - sqr(sat.ecc)) / (1 + sat.ecc * cos(ang));
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

	double switch_co(double in) //Changes cos x to sin x and vice versa. No guarantees as to the sign of the output.
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

	void do_orbit_phys(body &sat, double w_time, vec_n force_add)
	{
		double d_time = w_time - sat.t_l;
		body &parent = *sat.parent;
		vec_r pos_r = sat.pos - parent.pos;
		vec_n force_vec;

		double F = -parent.u / sqr(pos_r.mag);

		force_vec = vec_to_pos(pos_r.ang, F);

		sat.vel += force_vec * d_time;
		sat.vel += force_add * d_time;

		sat.pos += sat.vel * d_time;
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

	vector<vec_n> do_orbit(body &sat, double w_time)
	{
		double d_time = w_time - sat.epoch;

		body &parent = *sat.parent;

		double part = d_time * sat.Mn;
		double V;
		double radius;
		vector<vec_n> out(2);
		if (sat.isSun)
			return out;

		if (sat.shape == 0)
			part = ang_wrap(part);

		V = do_orbit_precise(part, sat, log2(pow(10, 10)));
		if (sat.inverse)
			V = M_2PI - V;

		radius = get_r(V, sat);

		out[0] = vec_to_pos(sat.normal - V, radius);

		double vel_mag = sqrt(2 * (sat.En + parent.u / radius));
		double vel_ang = acos(sat.Ar / (radius * vel_mag));
		if (V < M_PI)
			vel_ang = -vel_ang;

		if (sat.inverse)
			vel_mag = -vel_mag;

		vel_ang = ang_wrap(-V - vel_ang + sat.normal + M_PI2);

		out[1].x = -vel_mag * cos(vel_ang);
		out[1].y = -vel_mag * sin(vel_ang);

		return out;
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
		double vmag_end = HUGE_VAL;
		bool expire_true = false;

		for (int i = 0; i < list.size(); i++)
		{
			body sat_2 = (*(list[i]));
			double SOI = sat_2.SOI;
			double vmag_min = HUGE_VAL;
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

				double dist = vmag(pos - pos_2);

				if (expire_true && w_time_pot > w_time_end)
					break;

				if (dist < vmag_min)
				{
					vmag_min = dist;
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
				vmag_end = vmag_min;
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

			double r = vmag(sat.pos - parent_pot.pos);

			if (r < parent_pot.SOI)
				if (parent_pot.SOI < parent_cur.SOI)
					parent_presumed = list[i];
		}

		return parent_presumed;
	}

	void make_orbit(body &sat_in, double w_time)
	{
		body &sat = *sat_in.self;
		body &parent = *(sat.parent);

		double rel_ang, u;

		vec_r pos_r = sat.pos - parent.pos;
		vec_r vel_r = sat.vel - parent.vel;
		sat.inverse = false;
		u = parent.u;

		bool repeat;
		do
		{
			repeat = false;

			rel_ang = M_PI2 + vel_r.ang - pos_r.ang;
			rel_ang = ang_wrap(rel_ang, 2);

			if (abs(rel_ang) > M_PI2)
			{
				vel_r.mag *= -1;
				sat.inverse = true;
				repeat = true;
			}
			if (vel_r.mag == get_esc_vel(parent, pos_r.mag))
			{
				vel_r.mag = nextafter(vel_r.mag, 0);
				repeat = true;
			}

		} while (repeat);

		vector<dual_val> r_vec_vector = get_extremes(vel_r.mag, rel_ang, pos_r.mag, u);

		dual_val r_vec = r_vec_vector[0];
		dual_val r_vec_vel = r_vec_vector[1];

		sat.ax_a = (r_vec.two + r_vec.one) / 2;
		sat.ecc = (r_vec.two - r_vec.one) / (r_vec.two + r_vec.one);
		sat.ax_b = sqrt_a(1 - sqr(sat.ecc)) * sat.ax_a;
		sat.shape = sat.ecc > 1;
		sat.SOI = sat.ax_a * pow(sat.u / u, 0.4);


		double tV = get_V_r(pos_r.mag, sat);
		if (rel_ang < 0)
			tV = ang_wrap(-tV);

		sat.normal = pos_r.ang + tV;

		if (sat.inverse)
			tV = ang_wrap(-tV);
		sat.En = sqr(vel_r.mag) / 2 - u / pos_r.mag;			//Store orbital Energy, velocital energy plus gravitational energy
		sat.Ar = vel_r.mag * pos_r.mag * cos(rel_ang);			//Store area swept per unit of time, altitude times velocity times cosine-of-the-angle
		sat.Mn = Mn_from_r(sat.ax_a, u);					//Store mean angular motion of the body, as per yet another equation dug up from obscure wikipedia stubs
		sat.epoch = w_time - get_M(tV, sat.ecc) / sat.Mn;

		if (sat.shape || r_vec.two > parent.SOI)
		{
			double V_max = get_V_r(parent.SOI, sat);
			double M_max = get_M(V_max, sat.ecc);
			sat.expiry = M_max / sat.Mn + sat.epoch;
		}
	}



	vector<vec_n> do_phys_tick(vector<body*> bodies, double w_time, bool phys_mode = 0)
	{
		vector<vec_n> out(1);
		vector<vec_n> out2;
		bool emodelast = false;


		for (int i = 0; i < bodies.size(); i++)
		{
			body &sat = *bodies[i];
			vector<vec_n> pos_vel = do_orbit(sat, w_time);

			sat.pos = pos_vel[0] + (*sat.parent).pos;
			sat.vel = pos_vel[0] + (*sat.parent).vel;
			sat.t_l = w_time;

			out2.push_back(sat.pos);
		}

		//get_expiry(*bodies[bodies.size() - 1], bodies, w_time);
		return out2;
	}


	double zoom_mem = 1;

	void run_engine()
	{
		using shared::screen_state;

		zoom_mem *= pow(1.1, input::keyboard.scroll);

		screen_state.zoom = zoom_mem;
		screen_state.bodies = do_phys_tick(bodies, shared::r_time * 0.0001, 0);
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
		sun.parent = &sun;
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

	void engine_init()
	{
		phys_init();
	}

}