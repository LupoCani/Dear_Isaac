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
#define RENDER_DEBUG_INSTALLED true


namespace shared 					//Declare basic shared parameters
{
	sf::RenderWindow window2;
	bool window_is_clear = false;

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
		std::vector<int>  pressed;
		std::vector<int>  pressed_last;
		std::vector<bool> buttons;
		std::vector<bool> buttons_last;
		struct keys : sf::Keyboard {
			using Keyboard::Key;
		};
		struct btns : sf::Mouse {
			using Mouse::Button;
		};

		int isAt(int key, std::vector<int> list)
		{
			for (int i = 0; i < list.size(); i++)
				if (key == list[i])
					return i;

			return -1;
		}

		bool isPressed(int key)
		{
			return isAt(key, pressed) >= 0;
		}

		int pressedAt(int key)
		{
			return isAt(key, pressed);
		}

		bool wasPressed(int key)
		{
			return isAt(key, pressed) >= 0 && isAt(key, pressed_last) < 0;
		}

		bool msDown(short button)
		{
			return buttons[button] && !buttons_last[button];
		}
		bool msPressed(short button)
		{
			return buttons[button];
		}
	};

	key_state keyboard;

	void run_input()
	{
		Event input;
		keyboard.scroll = 0;
		keyboard.pressed_last = keyboard.pressed;
		keyboard.buttons_last = keyboard.buttons;
		keyboard.buttons = std::vector<bool>(Mouse::Button::ButtonCount);

		while (shared::window2.pollEvent(input)) {

			if (input.type == Event::KeyPressed)
				if (!keyboard.isPressed(input.key.code))
					keyboard.pressed.push_back(input.key.code);

			if (input.type == Event::KeyReleased)
				keyboard.pressed.erase(keyboard.pressed.begin() + keyboard.pressedAt(input.key.code));

			keyboard.buttons[Mouse::Button::Left] = Mouse::isButtonPressed(Mouse::Button::Left);
			keyboard.buttons[Mouse::Button::Right] = Mouse::isButtonPressed(Mouse::Button::Right);
			keyboard.buttons[Mouse::Button::Middle] = Mouse::isButtonPressed(Mouse::Button::Middle);
			keyboard.buttons[Mouse::Button::XButton1] = Mouse::isButtonPressed(Mouse::Button::XButton1);
			keyboard.buttons[Mouse::Button::XButton2] = Mouse::isButtonPressed(Mouse::Button::XButton2);

			if (input.type == Event::MouseWheelMoved)
				keyboard.scroll = input.mouseWheel.delta;
		}
	}
}

namespace phys					//Declare various classes and functions
{
	struct vec_n;
	struct vec_r;

	const long double M_PI = acos(-1.0L);		//Equal to π
	const long double M_2PI = 2L * M_PI;		//Equal to 2π
	const long double M_PI2 = M_PI / 2L;		//Equal to π/2
	const long double M_3PI2 = 3L * M_PI2;		//Equal to .75π
	const long double M_E = exp(1.0);

	double clamp(double in);
	double sqr(double in);
	double sqrt_a(double in);
	vec_n vec_to_pos(double v, double r);
	vec_n vec_to_pos(vec_r in);
	double atan2(vec_n pos);
	vec_n t_to_screen(vec_n in, double zoom = 1);
	double ang_wrap(double in, int quad = 4);
	double vmag(double x, double y);
	double vmag(vec_n in);
	vec_n rot_vec(vec_n in, double rot);
	double to_rad(double ang, double scale, double sc_out = M_2PI);
	using std::atan2;

	struct vec_r
	{
		double ang = 0;
		double mag = 0;
		operator vec_n();

		void normalize();
	};

	struct vec_n
	{
		double x = 0;
		double y = 0;

		vec_n(double in_x = 0, double in_y = 0)
		{
			x = in_x;
			y = in_y;
		}

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

	vec_r::operator vec_n()
	{
		return vec_to_pos(*this);
	}
	void vec_r::normalize()
	{
		vec_n norm = *this;
		*this = norm;
	}

}

namespace shared
{
	//struct vec_n : phys::vec_n {};//Copy the ever-so-important vec_n into shared
	using phys::vec_n;

	struct world_state
	{
		std::vector<vec_n> bodies;
		std::vector<std::string> names;
		std::vector<std::vector<vec_n>> paths;

		double zoom = 1;
		int focus = 0;
		int target = -1;
		double target_time;
		vec_n target_close;
		vec_n player_close;

		double player_rotation; //Angle in radians, sweeping counter-clockwise.
	};

	world_state screen_state;
}


namespace phys
{
	using std::vector;

	struct body
	{
		double ax_a;	//The length of the semimajor axis
		double ax_b;	//The length of the semiminor axis. So far, there's been no use for it, meaning it's sorta-undefined
		double ecc;		//The orbital eccentricity of the body
		double epoch;	//Timestamp of the latest passage through the periapsis, according to the current shape of the orbit
		double t_p;		//For elliptical orbits, the orbital period. Undefined for hyperbolic orbits.
		double t_l;		//Timestamp on the time-dependent paramters' current values
		double u;		//The body's SGP
		double normal;	//The angle of the orbit's periapsis. Counter-clockwise from +X
		double SOI;		//Radius of the sphere of influence
		short shape;	//Shape of the orbit. 0 denotes elliptical orbits, 1 denotes hyperbolic orbits
		bool inverse;	//Whether or not the orbit is counter-clockwise, with all the special cases that entails

		short expire;	//Whether or not there is a predicted expiry. 0 = no confirmed expiry, 1 = expiry by leaving SOI, 2 = expiry by entering SOI
		double expiry;	//Before this time, it is guaranteed no type-one expiry will occur. The time of expiry if expire > 0;
		double safe;	//Before this time, it is guaranteed no expiry will occur.
		double V_exp;	//The true anomaly at the predicted point of expiry.

		vec_n vel;		//Velocity at t_l
		vec_n pos;		//Position at t_l

		double En;		//Orbital energy
		double Ar;		//Area swept per unit of time. Unreliable for calculating Mn
		double Mn;		//Mean angular motion, aka mean anomaly anomaly swept per unit of time

		bool isPlayer = false;	//True if the body is the player
		bool isSun = false;		//True if the body is the sun

		body* parent;		//Pointer to parent body
		body* self = this;	//Pointer to self
		std::string name;	//Name of body
	};

	struct dual_val
	{
		long double one = 0;
		long double two = 0;
	};

	struct vars_struct
	{
		vector<body*> bodies;
		vector<vector<vec_n>> tails;
		double zoom_mem = 2;
		double last_predict = 0;
		double w_time = 0;
		double w_time_last = 0;
		double w_time_diff = 0;
		double d_time_fact = 1;

		int target = 0;
		bool isTarget = false;
		double min_dist;
		vec_n target_pos;
		vec_n player_pos;
	};

	struct plyr_struct
	{
		double spin = 0;
		double gyro = 0;
		double rot = 0;

		double thrust = 0;
		double fuel = 0;

		double mass = 0;
		double f_p_t = 1;
		double eng_F = 1;
	};

	plyr_struct plr;
	vars_struct gen;

	//Begin generic vector math functions
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

	vec_n vec_to_pos(vec_r in)
	{
		return vec_to_pos(in.ang, in.mag);
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

	double to_rad(double ang, double scale, double sc_out)
	{
		return sc_out * ang / scale;
	}

	double switch_co(double in) //Changes cos x to sin x and vice versa. No guarantees as to the sign of the output.
	{
		return sqrt(1 - pow(in, 2));
	}
	//End generic vector math functions

	//Begin orbital equation functions
	double get_r(double ang, body sat)
	{
		return sat.ax_a * (1 - sqr(sat.ecc)) / (1 + sat.ecc * cos(ang));
	}

	double get_V_r(double r, body sat)
	{
		return acos(clamp(sat.ax_a * (1 - sqr(sat.ecc)) / (sat.ecc * r) - 1 / sat.ecc));
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

		if (!(sat.shape || w_time_min < 0))
		{
			while (w_time < w_time_min)
				w_time += sat.t_p;
		}

		return w_time;
	}

	double Mn_from_V(double ecc, double r, double v)
	{
		double dV = sqrt(1 - sqr(ecc)) * (1 - ecc) / (1 + ecc);
		return dV * v / r;
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
		double part_root = sqrt(sqr(part_div) - 2 * u / r + pow(vel, 2));
		double part_area = cos(ang) * vel * r;

		out_v.one = part_div + part_root;
		out_v.two = part_div - part_root;

		out_r.one = part_area / out_v.one;
		out_r.two = part_area / out_v.two;

		out[0] = out_r;
		out[1] = out_v;

		return out;
	}

	vec_n get_pos_ang(double V, body sat)
	{
		vec_r point;
		point.ang = -V;
		point.mag = get_r(point.ang, sat);
		point.ang += (sat).normal;

		return point;
	}

	double get_V_phys(body sat, vec_n pos = vec_n())		//Get the current true anomaly from the physical position
	{
		if (vmag(pos) > 0)
			return sat.normal - atan2(pos);
		else
			return sat.normal - atan2(sat.pos - (*sat.parent).pos);
	}

	double V_bound_H(double ecc)
	{
		return acos(-1 / ecc);
	}

	//End orbital equation functions

	//Begin orbital algorithms
	double do_orbit_precise_H(double part, body sat, int precision)
	{
		double ceil = V_bound_H(sat.ecc);
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
			if (get_M_E(mid, sat.ecc) > part)
				ceil = mid;
			else
				floor = mid;

			mid = (ceil + floor) / 2;
		}

		return get_V(mid, sat.ecc);
	}

	void do_orbit_phys(body &sat, double w_time, vec_n force_add, vec_n &pos_out, vec_n &vel_out)
	{
		double d_time = w_time - sat.t_l;
		body &parent = *sat.parent;
		vec_r pos_r = sat.pos - parent.pos;
		vec_r vel_r = sat.vel - parent.vel;
		vec_n force_vec;

		double F = -parent.u / sqr(pos_r.mag);

		force_vec = vec_to_pos(pos_r.ang, F);

		vel_out = vel_r + force_vec * d_time;

		vel_out += force_add * d_time;

		pos_out = pos_r + vel_out * d_time;
	}

	void do_orbit(body &sat, double w_time, vec_n &pos_out, vec_n &vel_out = vec_n())
	{
		double d_time = w_time - sat.epoch;

		body &parent = *sat.parent;

		double part = time_to_M(sat, w_time);
		double V;
		double radius;
		if (sat.isSun)
			return;

		if (sat.shape == 0)
			part = ang_wrap(part);

		V = do_orbit_precise(part, sat, log2(pow(10, 20)));
		if (sat.inverse)
			V = ang_wrap(V);

		radius = get_r(V, sat);

		pos_out = get_pos_ang(V, sat);

		vec_r vel;

		vel.mag = sqrt(2 * (sat.En + parent.u / radius));
		vel.ang = acos(sat.Ar / (radius * vel.mag));
		if (V < M_PI)
			vel.ang = -vel.ang;

		if (sat.inverse)
			vel.mag = -vel.mag;

		vel.ang = -V -vel.ang + sat.normal + M_PI2;
		vel.ang = ang_wrap(vel.ang + M_PI);

		vel_out = vel;
	}

	void get_expiry(body &sat, vector<body*> list, vector<body*> &list_out, vector<vector<double>> &mins_out, bool &will_expire, int &new_parent)
		//This one's gonna be tough.
	{
		vector<body*> list_in = list;
		vector<double> starts;			//Current V-values for bodies
		vector<double> ends;			//V-values for bodies at the end of the relevant time span
		vector<int>	its;				//The iterators for the bodies
		list.clear();
		int precision;
		int precision_base = 1000;

		vector<vector<double>> pairs_dist;
		vector<vector<double>> pairs_time;

		vector<double> min_dists;
		vector<double> min_times;

		/*
		We'll try to find if there's any body, orbiting the same parent as the player, into whose SOI the player will enter, and when.
		We could just check a lot of moments in time with do_orbit(), but that function is expensive as hell.

		Thus, we increment 
		*/


		for (int i = 0; i < list_in.size(); i++)
		{
			body &pln = *list_in[i];
			if (pln.parent == sat.parent && !(pln.isPlayer || pln.isSun))
			{
				vec_n cur_pos = pln.pos - (*pln.parent).pos;;
				if (sat.safe - sat.t_l > 0)
					do_orbit(pln, sat.t_l, cur_pos), std::cout << "hi!\n";

				vec_n end_pos;
				do_orbit(pln, sat.expiry, end_pos);

				list.push_back(pln.self);
				starts.push_back(get_V_phys(pln, cur_pos));		//Get the current true anomaly of the body. 
				ends.push_back(get_V_phys(pln, end_pos));
				its.push_back(0);

				pairs_dist.push_back(vector<double>());
				pairs_time.push_back(vector<double>());
				min_times.push_back(DBL_MAX);
				min_dists.push_back(DBL_MAX);
			}
		}

		double start, end;
		if (list.size()) {
			vec_n sat_pos;
			do_orbit(sat, sat.t_l, sat_pos);
			start = get_V_phys(sat, sat_pos);

			do_orbit(sat, sat.expiry, sat_pos);
			end = get_V_phys(sat, sat_pos);
		}
		//double end = sat.V_exp;
		if (list.size())
			std::cout << M_to_time(sat, get_M(end, sat.ecc), sat.t_l-0.01) -M_to_time(*list[0], get_M(ends[0], (*list[0]).ecc), 0 ) << std::endl, system("pause");


		vec_n sat_pos_old;
		double sat_tim_old = 0;

		precision = precision_base * ang_wrap(end - start) / M_2PI;
		if (precision < precision_base / 20)
			precision = precision_base / 20;

		if (list.size())
		{
			std::cout << list.size() << ": \n";

			for (int i2 = 0; i2 < precision; i2++)
			{
				int i_sat = 0;
				int i_pln = 0;
				double Vs = to_rad(i2, precision, end - start);

				vec_n sat_pos = get_pos_ang(start + Vs, sat);
				double sat_tim = M_to_time(sat, get_M(Vs, sat.ecc), sat.t_l);

				vec_n sat_pos_diff = sat_pos - sat_pos_old;
				double sat_tim_diff = sat_tim - sat_tim_old;

				std::cout << i2 << ": " << std::endl;

				for (int i = 0; i < list.size(); i++)
				{
					body &pln = *list[i];
					double Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
					double pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);
					vec_n pln_pos = get_pos_ang(Vp, pln);

					std::cout << "    " << i << std::endl;

					while (pln_tim < sat_tim_old)
					{
						its[i]++;
						Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
						pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);
					}
					std::cout << "    " << "Midway!" << std::endl;
					while (pln_tim < sat_tim)
					{
						double diff_part = (pln_tim - sat_tim_old) / sat_tim_diff;	//What percentage of the current time interval has passed in the moment we're observing.
						vec_n sat_pos_est = sat_pos_old + sat_pos_diff * diff_part;

						double dist_est = vmag(sat_pos_est - pln_pos);
						pairs_dist[i].push_back(dist_est);
						pairs_time[i].push_back(pln_tim);

						its[i]++;
						Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
						pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);

						std::cout << "    " << Vp << std::endl;
						std::cout << "    " << get_M(Vp, pln.ecc) << std::endl;
						std::cout << "    " << pln_tim << " - " << sat_tim << std::endl;
						std::cout << "    " << " ----- " << std::endl;
						Sleep(20);

					}
				}
				sat_pos_old = sat_pos;
				sat_tim_old = sat_tim;
			}
		}
			

		body *parent_next_p = (*sat.parent).parent;
		double end_body = -1;
		double end_time = sat.expiry;
		bool ending = false;


		for (int i = 0; i < list.size(); i++)
		{
			for (int i2 = 0; i2 < pairs_time.size(); i2++)
			{
				std::cout << i << " " << i2 << std::endl;
				if (pairs_dist[i][i2] 
					< min_dists[i])
				{
					min_dists[i] = pairs_dist[i][i2];
					min_times[i] = pairs_time[i][i2];
				}
				if (pairs_dist[i][i2] < (*list[i]).SOI)
					if (pairs_time[i][i2] < end_time)
					{
						ending = true;
						end_body = i;
						break;
					}
			}
		}

		list_out = list;
		will_expire = ending;
		new_parent = end_body;

		mins_out.push_back(min_dists);
		mins_out.push_back(min_times);
	}

	body* get_parent(body &sat, vector<body*> list)
	{
		body *parent_presumed = list[0];
		if (sat.isSun)
			return sat.self;

		for (int i = 0; i < list.size(); i++)
		{
			body &parent_pot = *list[i];
			body &parent_cur = *parent_presumed;

			double r = vmag(sat.pos - parent_pot.pos);

			if (list[i] == sat.self || (*list[i]).u <= sat.u)
				continue;

			if (r < parent_pot.SOI)
				if (parent_pot.SOI < parent_cur.SOI)
					parent_presumed = list[i];
		}
		return parent_presumed;
	}

	void make_orbit(body &sat, double w_time)
	{
		body &parent = *sat.parent;
		if (sat.isSun)
			return;

		double rel_ang;							//The angle between the velocity and the perpendicular of the radius
		double u;								//The SGP of the parent body

		vec_r pos_r = sat.pos - parent.pos;		//The radial vector position of the sat, relative the parent
		vec_r vel_r = sat.vel - parent.vel;		//The radial vector velocity of the sat, relative the parent
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
				vel_r.normalize();
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

		double tV = get_V_r(pos_r.mag, sat);	//The current true anomaly. Sweeps counter-clockwise
		if (rel_ang < 0)
			tV = ang_wrap(-tV);

		sat.normal = pos_r.ang + tV;

		if (sat.inverse)
			tV = ang_wrap(-tV);

		sat.En = sqr(vel_r.mag) / 2 - u / pos_r.mag;			//Store orbital Energy, velocital energy plus gravitational energy
		sat.Ar = vel_r.mag * pos_r.mag * cos(rel_ang);			//Store area swept per unit of time, altitude times velocity times cosine-of-the-angle
		sat.Mn = Mn_from_r(sat.ax_a, u);						//Store mean angular motion of the body, as per yet another equation dug up from obscure wikipedia stubs
		sat.epoch = w_time - get_M(tV, sat.ecc) / sat.Mn;
		sat.t_p = M_2PI / sat.Mn;
		sat.expiry = w_time + sat.t_p;
		sat.safe = w_time;
		sat.expire = 0;

		if (sat.shape || r_vec.two > parent.SOI)
		{
			double V_max = get_V_r(parent.SOI, sat);
			double M_max = get_M(V_max, sat.ecc);
			sat.expiry = M_to_time(sat, M_max, sat.epoch);
			sat.expire = 1;
			sat.V_exp = V_max;
		}
	}

	//End orbital algorithms

	//Begin graph algorithms

	vector<vec_n> make_tail(body sat, int subdiv)
	{
		vector<vec_n> out;

		double V_start = ang_wrap(get_V_phys(sat), 2);
		double V_end;

		if (sat.expire)
			V_end = sat.V_exp;
		else
			V_end = V_start + M_2PI;

		V_end = nextafter(V_end, 0);

		double V_span = V_end - V_start;

		for (int i2 = 0; i2 <= subdiv; i2++)
			out.push_back(get_pos_ang(V_start + to_rad(i2, subdiv, V_span), sat));

		return out;
	}

	//End graph algorithms

	vector<vec_n> do_phys_tick(vector<body*> bodies, double w_time, bool phys_mode, vec_n thrust = vec_n())
	{
		vector<vec_n> out;
		bool expired = false;

		vector<vec_n> pos_buffer;
		vector<vec_n> vel_buffer;

		for (int i = 0; i < bodies.size(); i++)
		{
			body &sat = *bodies[i];
			vec_n pos;
			vec_n vel;
			if (!phys_mode || !sat.isPlayer)
				do_orbit(sat, w_time, pos, vel);
			else
				do_orbit_phys(sat, w_time, thrust, pos, vel);

			pos_buffer.push_back(pos);
			vel_buffer.push_back(vel);

		}
		for (int i = 0; i < bodies.size(); i++)
		{
			body &sat = *bodies[i];

			sat.pos = pos_buffer[i] + (*sat.parent).pos;
			sat.vel = vel_buffer[i] + (*sat.parent).vel;
			sat.t_l = w_time;


			if (!sat.expire)
				sat.expiry = sat.t_l + sat.t_p;
			else
			{
				body &new_parent = *get_parent(sat, bodies);
				if (new_parent.self != sat.parent)
					expired = true;
				sat.parent = &new_parent;
			}
			out.push_back(sat.pos);
		}

		body &plyr = *bodies.back();

		if (phys_mode || expired || plyr.expire)
		{
			if (expired || phys_mode)
				make_orbit(plyr, w_time);

			if (expired || phys_mode || plyr.safe < plyr.t_l)
				plyr.safe = plyr.t_l;

			gen.tails[gen.tails.size() - 1] = make_tail(plyr, 1000);
		}


		//get_expiry(*bodies[bodies.size() - 1], bodies, w_time);
		return out;
	}

	vec_n do_game_tick()
	{
		using namespace input;

		if (!keyboard.isPressed(key_state::keys::RShift))
			gen.zoom_mem *= pow(1.1, keyboard.scroll);
		else
			gen.d_time_fact *= pow(1.1, keyboard.scroll);

		plr.spin -= plr.gyro * gen.w_time_diff * keyboard.isPressed(key_state::keys::Left);
		plr.spin += plr.gyro * gen.w_time_diff * keyboard.isPressed(key_state::keys::Right);

		plr.rot += gen.w_time_diff * plr.spin;

		plr.rot = ang_wrap(plr.rot, 2);

		plr.thrust += keyboard.isPressed(key_state::keys::Up);
		plr.thrust -= keyboard.isPressed(key_state::keys::Down);

		if (keyboard.isPressed(key_state::keys::Numpad1))
		{
			plr.spin = 0;
		}
		if (keyboard.isPressed(key_state::keys::Numpad3))
		{
			gen.d_time_fact = 100;
		}
		if (keyboard.isPressed(key_state::keys::Numpad6))
		{
			gen.d_time_fact = 1;
		}
		if (keyboard.isPressed(key_state::keys::Numpad2))
		{
			plr.thrust = 0;
		}
		if (keyboard.isPressed(key_state::keys::Numpad0))
		{
			std::cin >> gen.target;
		}
		if (keyboard.isPressed(key_state::keys::Numpad4))
		{
			std::cin >> shared::screen_state.focus;
		}

		double thrust_actual = 0;
		vec_r accel;

		if (plr.thrust < 0)
			plr.thrust = 0;

		if (plr.fuel > 0)
			thrust_actual = plr.thrust * plr.eng_F;

		plr.fuel -= thrust_actual * plr.f_p_t;
		if (plr.fuel < 0)
			plr.fuel = 0;

		accel.mag = thrust_actual / (plr.mass + plr.fuel);
		accel.ang = plr.rot;

		using shared::screen_state;
		screen_state.focus = 6;
		screen_state.zoom = gen.zoom_mem;
		screen_state.player_rotation = plr.rot;

		return vec_to_pos(plr.rot, -thrust_actual / (plr.mass + plr.fuel) );
	}

#ifdef RENDER_DEBUG_INSTALLED
	bool emode = 0;
	vec_n thrust_debug;

#endif

	void run_engine()
	{
		using shared::screen_state;

		gen.w_time_last = gen.w_time;
		double diff_time = (shared::r_time - shared::l_time) / 100.0 / shared::cps / gen.d_time_fact;
		if (diff_time < 0.01)
			gen.w_time += diff_time;
		else
			gen.w_time += 0.01;
		gen.w_time_diff = gen.w_time - gen.w_time_last;

		vec_n eng_thrust = do_game_tick();
		bool eng_mode = vmag(eng_thrust);

		screen_state.bodies = do_phys_tick(gen.bodies, gen.w_time, eng_mode, eng_thrust);

		if (gen.last_predict + 0.1 / shared::cps < shared::r_time)
		{
			vector<vector<double>> pairs;
			body &plyr = *gen.bodies.back();
			vector<body*> co_sats;
			bool expiring;
			int new_parent;
			double expire_time;

			get_expiry(plyr, gen.bodies, co_sats, pairs, expiring, new_parent);

			if (expiring)
			{
				plyr.expire = 2;
				plyr.expiry = pairs[0][new_parent];
			}
			else
			{
				plyr.safe = plyr.expiry;
			}

			//std::cout << co_sats.size() << "found one!\n";

			for (int i = 0; i < co_sats.size(); i++)
			{
				if (co_sats[i] == gen.bodies[gen.target])
				{
					gen.min_dist = pairs[0][i];
				}
			}


			gen.last_predict = shared::r_time;
		}

		vector<vector<vec_n>> tails_out = gen.tails;

		for (int i = 0; i < tails_out.size(); i++)
		{
			vec_n par_pos = (*(*gen.bodies[i]).parent).pos;
			for (int i2 = 0; i2 < tails_out[i].size(); i2++)
			{
				tails_out[i][i2] += par_pos;
			}
		}

		screen_state.paths = tails_out;

#ifdef RENDER_DEBUG_INSTALLED
		emode = eng_mode;
		thrust_debug = eng_thrust;
#endif
	}

	vector<body*> sort_bodies(vector<body*> unsorted)
	{
		vector<body*> sorted;
		for (int i = 0; i < unsorted.size(); i++)
			(*unsorted[i]).self = unsorted[i];

		{
			while (unsorted.size() > 0)
			{
				int   body_i = 0;
				body *body_p = unsorted[body_i];

				for (int i = 0; i < unsorted.size(); i++)
					if ((*unsorted[i]).u >(*body_p).u)
					{
						body_p = unsorted[i];
						body_i = i;
					}
				sorted.push_back(body_p);
				unsorted.erase(unsorted.begin() + body_i);
			}

			//Name the largest body as 'sun';
			(*sorted[0]).isSun = true;
			(*sorted[0]).SOI = DBL_MAX;
			//To begin with, assume all sorted orbit the sun.
			for (int i = 0; i < sorted.size(); i++)
				(*sorted[i]).parent = sorted[0];
		}

		return sorted;
	}

	vector<vector<vec_n>> get_tails_basic(vector<body*> list, int subdiv = 100)
	{
		vector<vector<vec_n>> paths;

		for (int i = 0; i < list.size(); i++)
		{
			paths.push_back(make_tail(*list[i], subdiv));
		}
		return paths;
	}

	void phys_init()
	{
		body& sun = *new body;
		sun.pos.x = 0;
		sun.pos.y = 0;
		sun.vel.x = 0;
		sun.vel.y = 0;
		sun.u = pow(10, 14);
		sun.name = "Sun";

		body& moon = *new body;
		moon.pos.x = -25300;
		moon.pos.y = 0;
		moon.vel.x = 0;
		moon.vel.y = 46000;
		moon.u = pow(10, 8);
		moon.name = "moon";

		body& planet = *new body;
		planet.pos.x = -25000;
		planet.pos.y = 0;
		planet.vel.x = 0;
		planet.vel.y = 40000;
		planet.u = pow(10, 10);
		planet.name = "planet";

		body& pluto = *new body;
		pluto.pos.x = -30000;
		pluto.pos.y = -0;
		pluto.vel.x = 0;
		pluto.vel.y = 40000;
		pluto.u = pow(10, 10);
		pluto.name = "pluto";

		body& dune = *new body;
		dune.pos.x = -35000;
		dune.pos.y = 0;
		dune.vel.x = 0;
		dune.vel.y = 40000;
		dune.u = pow(10, 10);
		dune.name = "dune";

		body& yavin = *new body;
		yavin.pos.x = -40000;
		yavin.pos.y = 0;
		yavin.vel.x = 0;
		yavin.vel.y = 40000;
		yavin.u = pow(10, 10);
		yavin.name = "yavin";

		body& plyr = *new body;
		plyr.pos.x = -40200;
		plyr.pos.y = 0;
		plyr.vel.x = 0;
		plyr.vel.y = 50000;
		plyr.u = 1;
		plyr.isPlayer = true;
		plyr.name = "plyr";

		vector<body*> unsorted;

		unsorted.push_back(&sun);		//White
		unsorted.push_back(&moon);		//Blue
		unsorted.push_back(&planet);	//Yellow
		unsorted.push_back(&pluto);		//Green
		unsorted.push_back(&dune);		//Cyan
		unsorted.push_back(&yavin);		//Magenta
		unsorted.push_back(&plyr);		//Red

		gen.bodies = sort_bodies(unsorted);

		for (int i2 = 0; i2 < 5; i2++)
			for (int i = 0; i < gen.bodies.size(); i++)
			{
				make_orbit(*gen.bodies[i], 0);								//Get everything's orbit around the sun
				(*gen.bodies[i]).parent = get_parent(*gen.bodies[i], gen.bodies);	//Check if the newly-determined SOIs make them into moons
			}

		gen.tails = get_tails_basic(gen.bodies);

		plr.fuel = 1000;
		plr.mass = 1000000000;
		plr.f_p_t = 0.0;
		plr.gyro = 100000000;
		plr.eng_F = 1000000000000;
	}

	void engine_init()
	{
		phys_init();
	}
}

#ifdef RENDER_DEBUG_INSTALLED
namespace render_debug			//To be removed once the neccesary render_tools functions are implemented
{
	using sf::Vertex;
	using sf::Color;
	using sf::CircleShape;
	using sf::Vector2f;
	using std::vector;
	using shared::vec_n;
	using shared::window2;
	sf::Font debug_font;
	bool nonsense = debug_font.loadFromFile("arial.ttf");

	bool window_is_clear = false;

	vector<vec_n> handle_scale(vector<vec_n> list, vec_n origo, double scale = 1, double mid_x = 500, double mid_y = 500)
	{
		for (int i = 0; i < list.size(); i++)
		{
			list[i] -= origo;
			list[i] *= scale;
			list[i].y *= -1;

			list[i] += vec_n(mid_x, mid_y);
		}

		return list;
	}

	void render_line(vector<vec_n> in, vec_n origo, double zoom)
	{
		sf::VertexArray lines(sf::LinesStrip, in.size());

		in = handle_scale(in, origo, zoom);

		for (int i = 0; i < in.size(); i++)
		{
			lines[i].position = in[i];
			lines[i].color = Color::Cyan;
		}

		window2.draw(lines);
	}
	void render_lines(vector<vector<vec_n>> in, vec_n origo, double zoom)
	{
		for (int i = 0; i < in.size(); i++)
		{
			render_line(in[i], origo, zoom);
		}
	}

	void render_text(std::string text, vec_n coords)
	{
		sf::Text FPS2;
		FPS2.setFont(debug_font);
		FPS2.setPosition(coords);
		FPS2.setCharacterSize(25);
		FPS2.setString(text);
		window2.draw(FPS2);
	}

	void render_texts(std::vector<std::string> texts)
	{
		for (int i = 0; i < texts.size(); i++)
			render_text(texts[i], vec_n(0, i * 25));
	}

	void render_player(double rot, vec_n pos, vec_n origo, double zoom)
	{
		vector<vec_n> temp;
		temp.push_back(pos);
		pos = handle_scale(temp, origo, zoom)[0];

		CircleShape plyr(10, 3);
		plyr.setFillColor(sf::Color::Green);
		plyr.setOrigin(5, 5);
		plyr.setPosition(pos);
		plyr.setScale(vec_n(1, 2));
		plyr.setRotation(rot / phys::M_2PI * 360.0 + 90.0);

		window2.draw(plyr);
	}

	void render_all(shared::world_state in)
	{
		window2.clear();
		window_is_clear = true;
		render_lines(in.paths, in.bodies[in.focus], in.zoom);

		std::vector<std::string> texts;
		texts.push_back(std::to_string(phys::plr.thrust));
		texts.push_back(std::to_string(phys::to_rad(phys::plr.rot, phys::M_2PI, 4)));
		texts.push_back(std::to_string(phys::to_rad(phys::atan2(phys::thrust_debug), phys::M_2PI, 4)));

		texts.push_back(std::to_string((*phys::gen.bodies.back()).ecc));

		texts.push_back(std::to_string(phys::gen.min_dist));

		render_texts(texts);
		render_player(in.player_rotation, in.bodies.back(), in.bodies[in.focus], in.zoom);

		//window2.display();
	}
}
#endif // RENDER_DEBUG_INSTALLED
