#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <windows.h>
#include <ciso646>
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

	//namespace world_state;
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
		bool isPressed_any()
		{
			return pressed.size();
		}
		bool wasPressed_any()
		{
			return pressed.size() && !pressed_last.size();
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

	namespace win_inf
	{
		bool update = 0;
		bool close_h = 1;
	}

	void run_input()
	{
		Event input;
		keyboard.scroll = 0;
		keyboard.pressed_last = keyboard.pressed;
		keyboard.buttons_last = keyboard.buttons;
		keyboard.buttons = std::vector<bool>(Mouse::Button::ButtonCount);

		keyboard.buttons[Mouse::Button::Left] = Mouse::isButtonPressed(Mouse::Button::Left);
		keyboard.buttons[Mouse::Button::Right] = Mouse::isButtonPressed(Mouse::Button::Right);
		keyboard.buttons[Mouse::Button::Middle] = Mouse::isButtonPressed(Mouse::Button::Middle);
		keyboard.buttons[Mouse::Button::XButton1] = Mouse::isButtonPressed(Mouse::Button::XButton1);
		keyboard.buttons[Mouse::Button::XButton2] = Mouse::isButtonPressed(Mouse::Button::XButton2);

		while (shared::window2.pollEvent(input)) {

			if (input.type == Event::KeyPressed)
			{
				if (!keyboard.isPressed(input.key.code))
					keyboard.pressed.push_back(input.key.code);
			}

			else if (input.type == Event::KeyReleased)
				keyboard.pressed.erase(keyboard.pressed.begin() + keyboard.pressedAt(input.key.code));

			else if (input.type == Event::MouseWheelMoved)
				keyboard.scroll = input.mouseWheel.delta;

			else if (input.type == Event::Resized)
				win_inf::update = true;
			else if (input.type == Event::Closed)
				win_inf::close_h = true;
		}
	}

	namespace flush_back
	{
		bool play = false;
		int men_cmd = 0;
		int pause_cmd = 0;
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

	double clamp(double in, double max = 1.0, double min = -1.0);
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
	double to_rad(double ang, double scale);
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
		*this = vec_n(*this);
	}

}

namespace shared
{
	using phys::vec_n;

	short game_state = 0;
	/*
	0: start menu
	1: playing
	2: paused
	4: game over
	5: closing
	*/

	namespace world_state		//Data for rendering/ui functions to use in rendering
	{
		double world_time = 0;

		//Physical variables:

		/**/std::vector<vec_n> bodies;				//List of bodies' positions. The first value is the sun, the last is the player.
		std::vector<double> sizes;				//List of bodies' sizes. Specifically, unzoomed radius.
		std::vector<vec_n> kesses;				//List of kessler fields
		std::vector<double> sizes_kess;			//List of kessler fields' sizes. Specifically, the unzoomed radius
		
		std::vector<std::string> names;			//List of the names for every corresponding body
		/**/std::vector<std::vector<vec_n>> paths;	//List of paths. Every path is a vector of points to be rendered as connected lines.

		//Screen variables:

		double zoom = 1;	//The current zoom level. Larger values mean more close zoom.
		/**/int focus = 0;		//The body the view focuses on.
		/**/int target = -1;	//The body selected by the player to target.

		int intercept = -1;		//The body the player is currently intercepting
		bool cepting = false;	//Whether or not the player is set to intercept a body
		double cept_time = 0;	//Time until projected intercept

		//Target approach variables:

		std::vector<double> target_time(1);	//The time of the two closest approaches the player will make
		std::vector<double> target_min(1);
		std::vector<vec_n> target_close(1);	//The target's position in the two closest approaches
		std::vector<vec_n> player_close(1);	//The player's position in the two closest approaches
		
		//Assist variables:

		std::vector<vec_n> min_pos;	//A list of minimum altitudes for existing bodies.
		std::vector<vec_n> max_pos;	//A list of maximum altitudes for existing bodies. 

		//Player status variables:

		double player_rotation;		//Angle in radians, sweeping counter-clockwise.
		double player_thrust;		//The current amount of engine thrust.
		double player_thrust_max;	//The maximum amount of engine thrust.

		//Game variables:

		/**/vec_n  goal_pos;	//The center of the target zone.
		/**/double goal_rad;	//The radius of the target zone.

		/**/int goal_parent;	//The id of the body around which the target is placed.
		/**/int goal_parent_2;	//If the above is a moon, this is the id of the planet the moon orbits.

		//Score varibles:

		/**/double eng_secs = 0;	//Abstract amount of time engine has been on
		/**/double goal_count = -1;	//Amount of goals that have been collected
		/**/double score = 0;		//Amount of points, represented by the score

		/**/double health;			//Current health
		/**/double health_max;		//Size of health bar
	};
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

		bool exiting;	//Whether or not there is a predicted expiry by leaving the SOI.
		double exit;	//The time of such an exit.
		double V_exit;

		bool entering;	//Whether or not there is a predicted expiry by entering an SOI.
		double entry;	//The time of such an entry.
		double safe;	//Before this time, it is guaranteed no expiry will occur.
		double V_ent;	//The true anomaly at the 

		short expire()	//The predicted point of expiry
		{
			return entering ? 2 : exiting;
		}
		double expiry()
		{
			return entering ? entry : (exiting ? exit : t_l + t_p);
		}
		double V_exp()	//The true anomaly at the predicted point of expiry.
		{
			return entering ? V_ent : V_exit;
		}

		vec_n vel;		//Velocity at t_l
		vec_n pos;		//Position at t_l

		double En;		//Orbital energy
		double Ar;		//Area swept per unit of time. Unreliable for calculating Mn
		double Mn;		//Mean angular motion, aka mean anomaly anomaly swept per unit of time

		bool isPlayer = false;	//True if the body is the player
		bool isSun = false;		//True if the body is the sun
		bool isKess = false;

		body* parent;		//Pointer to parent body
		body* self = this;	//Pointer to self
		std::string name;	//Name of body

		void split()
		{
			self = this;
		}

		double size = 1;
	};

	struct dual_val
	{
		long double one = 0;
		long double two = 0;
	};

	struct pred_p				//A struct for storing local extremes of planetary distances, be they maximal or minimal
	{
		double time = DBL_MAX;	//The time of the distance
		double dist = DBL_MAX;	//The distance itself
		short dlt = 1;			//The nature of the extreme point. 1 for minimals, -1 for maximals.

		pred_p(double dist_in = DBL_MAX, double time_in = DBL_MAX, int dlt_in = 1)
		{
			dist = dist_in;
			time = time_in;
			dlt = dlt_in;
		}
	};

	namespace gen
	{
		double w_time = 0;
		double w_time_last = 0;
		double w_time_diff = 0;
		double d_time_fact = 1;

		double last_predict = 0;

		vector<body*> bodies;
		vector<body*> kesses;
		vector<double> bodies_size;
		vector<double> kesses_size;
		vector<vector<vec_n>> tails;
		vector<vector<vec_n>> tails_future;
		vector<vector<vec_n>> tails_out;
		vector<vec_n> bodies_pos;
	}

	namespace game
	{
		int target = 0;
		bool isTarget = false;
		double min_dist;
		double min_time;
		double cur_dist;
		vec_n target_pos;
		vec_n player_pos;

		double zoom_mem = 2;
		int focus = 1;

		double health;
		double health_max;
		double eng_secs;

		double dmg_rad;
		double dmg_int;

		vec_n goal_coords;
		double goal_size;
		bool goal_active;
		int goal_id;
		int goal_parent;
		int goal_count = -1;

		vector<vec_n> kesses_pos;
	}

	namespace rock
	{
		double spin = 0;
		double gyro = 0;
		double rot = 0;

		double thrust = 0;
		double fuel = 0;

		double mass = 0;
		double f_p_t = 1;
		double eng_F = 1;

		vec_n eng_thrust;
		bool eng_mode;
	}

	namespace pred
	{
		double begin;
		double end;

		bool invalid = false;

		short expire_count = 0;

		int next_parent;	//The id of the predicted parent

		vec_n enter_plyr_pos;
		vec_n enter_pln_pos;
	}

	template<class vec_cont>
	void push_any(vector<vec_cont> &in, int ind, vec_cont element)
	{
		vector<vec_cont>::iterator ref = in.begin();
		if (ind < 0)
			ref = in.end();

		in.insert(ref + ind, element);
	}

	template<class vec_cont>
	void pop_any(vector<vec_cont> &in, int ind)
	{
		vector<vec_cont>::iterator ref = in.begin();
		if (ind < 0)
			ref = in.end();

		in.erase(ref + ind);
	}

	template<class vec_cont>
	int find_in(vector<vec_cont> &in, vec_cont element)
	{
		for (int i = 0; i < in.size(); i++)
			if (in[i] == element)
				return i;

		return -1;
	}

	double rand_part(double scale = 1)
	{
		return rand() % 1000 / 1000.0 * scale;
	}

	//Begin generic vector math functions
	double clamp(double in, double max, double min)
	{
		return in > max ? max : (in < min ? min : in);
	}

	double ang_scale(double scale, double ang, double scope)
	{
		return scale * ang / scope;
	}

	double sqr(double in)	//Because writing pow(num, 2) all the time got annoying.
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

	double ang_wrap_old(double in, int quad)
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
	double ang_wrap(double in, int quad)
	{
		double interval = M_2PI;
		double max = double(quad) * M_PI2;
		double min = max - interval;

		if (in > max)
		{
			in -= ceil((in - max) / interval) * interval;
		}
		if (in < min)
		{
			in += ceil((min - in) / interval) * interval;
		}


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

	vec_r rot_vec(vec_r in, double rot)
	{
		in.ang += rot;
		return in;
	}

	double to_rad(double ang, double scale)
	{
		return M_2PI * ang / scale;
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
		return acos(clamp(
			sat.ax_a * (1 - sqr(sat.ecc)) / (sat.ecc * r) - 1 / sat.ecc
		));
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

		tV = ang_wrap(tV);

		double cosE = get_cosE(tV, ecc);
		double sinE = switch_co(cosE);

		if (tV < M_PI)
			return acos(cosE) - ecc * sinE;
		else
			return M_2PI - acos(cosE) - ecc * -sinE;
	}

	double get_dM_E(double E, double ecc)
	{
		return 1 - ecc * cos(E);
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

	double M_to_time(body sat, double M, double w_time_min = NAN, double digits = 16)
	{
		double d_time = M / sat.Mn;
		double w_time = d_time + sat.epoch;
		double diff = pow(10, -digits);

		if (!(sat.shape || isnan(w_time_min) ))
		{
			if (w_time_min > (w_time + diff))
			{
				w_time += ceil((w_time_min - (w_time + diff)) / sat.t_p) * sat.t_p;
			}

		}

		return w_time;
	}

	double V_to_time(double V, body sat, double w_time_min = NAN)
	{
		return M_to_time(sat, get_M(V, sat.ecc), w_time_min);
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

	vector<double> get_extremes(double vel, double ang, double r, double u)
	{
		vector<double> out_v(2), out_r(2);

		double part_div = u / (cos(ang) * vel * r);
		double part_root = sqrt(sqr(part_div) - 2 * u / r + sqr(vel));
		double part_area = cos(ang) * vel * r;

		out_v[0] = part_div + part_root;
		out_v[1] = part_div - part_root;

		out_r[0] = part_area / out_v[0];
		out_r[1] = part_area / out_v[1];

		return out_r;
	}

	vec_n get_pos_ang(double V, body sat)
	{
		vec_r point;
		point.ang = V;
		point.mag = get_r(point.ang, sat);
		point.ang += sat.normal;

		return point;
	}

	double get_V_phys(body sat, vec_n pos = vec_n(NAN, NAN))		//Get the current true anomaly from the physical position
	{
		if (!isnan(pos.x))
			return atan2(pos) - sat.normal;
		else
			return  atan2(sat.pos - (*sat.parent).pos) - sat.normal;
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
			if (get_M_H(mid, sat.ecc) > part)
				ceil = mid;
			else
				floor = mid;

			mid = (ceil + floor) / 2;
		}

		return ang_wrap(mid);
	}

	double do_orbit_precise_N(double part, body sat, int precision)
	{
		if (sat.shape) return do_orbit_precise_H(part, sat, precision);

		double E_pot = part;
		double M_pot = 0;
		double diff = pow(2, -precision);

		for (int i = 0; i < precision *2; i++)
		{
			M_pot = get_M_E(E_pot, sat.ecc);
			E_pot = E_pot - M_pot / get_dM_E(E_pot, sat.ecc);

			if (abs(M_pot - part) < diff)
				break;
		}
		return get_V(E_pot, sat.ecc);
	}

	double do_orbit_precise(double part, body sat, int precision)
	{
		if (sat.shape) return do_orbit_precise_H(part, sat, precision);

		double ceil = M_2PI;
		double floor = 0;
		double mid = M_PI;
		double M_pot;

		for (int i = 0; i < precision; i++)
		{
			M_pot = get_M_E(mid, sat.ecc);	//Try a potential value of E

			if (M_pot > part)				//If the corresponding M is too large
				ceil = mid;					//halv the search space by lowering the upper bound.
			else if (M_pot < part)			//If it's too small,
				floor = mid;				//half the search space by increasing the lower bound.
			else
				break;						//In the unlikely event that the value should correspond exactly, return it.

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

	void do_orbit(body &sat, double w_time, int digits, vec_n &pos_out, vec_n &vel_out = vec_n(NAN, NAN))
	{
		double d_time = w_time - sat.epoch;

		body &parent = *sat.parent;

		double part = time_to_M(sat, w_time);
		double V;
		double radius;
		if (sat.isSun)
			return;

		if (!sat.shape)
			part = ang_wrap(part);

		V = do_orbit_precise(part, sat, log2(pow(10, digits)));

		radius = get_r(V, sat);

		pos_out = get_pos_ang(V, sat);

		if (isnan(vel_out.x)) return;

		vec_r vel;

		vel.mag = sqrt(2 * (sat.En + parent.u / radius));
		vel.ang = acos(sat.Ar / (radius * vel.mag));
		if (V > M_PI)
			vel.ang = -vel.ang;

		if (!sat.inverse)
			vel.mag = -vel.mag;

		vel.ang = V -vel.ang + sat.normal + M_PI2;
		vel.ang = ang_wrap(vel.ang + M_PI);

		vel_out = vel;
	}

	struct get_expiry
	{
		//With this many return variables, it's easier to make it a single-function class
		vector<body*> list_out;
		vector<vector<pred_p>> mins_out;
		bool will_expire;
		double V_exp;
		double expire_time;
		int new_parent;

		get_expiry(body &sat, vector<body*> list)
			//This one's gonna be tough.
		{
			vector<body*> list_in = list;
			vector<double> starts;			//Current V-values for bodies
			vector<double> ends;			//V-values for bodies at the end of the relevant time span
			vector<double> diffs;
			list.clear();
			int precision;
			int precision_base = 1000;

			/*
			We'll try to find if there's any body, orbiting the same parent as the player, into whose SOI the player will enter, and when.
			We could just check a lot of moments in time with do_orbit(), but that function is expensive as hell.
			Instead, we use the inverse process. We check an interval of values of V for the player
			*/

			double start, end, diff, start_t, end_t;

			start_t = sat.safe;
			end_t = sat.expiry();

			double fidelty = 20;
			int t_fid = 15;

			vec_n sat_pos, sat_pos_old;
			double sat_tim, sat_tim_old;
			{
				vec_n cur_pos, end_pos;

				do_orbit(sat, start_t, fidelty, cur_pos);
				do_orbit(sat, end_t, fidelty, end_pos);

				start = get_V_phys(sat, cur_pos);
				end = get_V_phys(sat, end_pos);

				start = ang_wrap(start);

				double laps = 0;
				if (sat.ecc < 0.9999999999)
					laps = floor((end_t - start_t) / sat.t_p) * M_2PI;

				if (end == start + laps)
					end += M_2PI;

				if (end < start + laps)
					end += ceil(((start + laps) - end) / M_2PI) * M_2PI ;

				diff = end - start;
				sat_pos = cur_pos;
				sat_tim = start_t;
			}


			for (int i = 0; i < list_in.size(); i++)
			{
				body &pln = *list_in[i];
				if (pln.parent == sat.parent && !(pln.isPlayer || pln.isSun))
				{
					vec_n cur_pos, end_pos;
					do_orbit(pln, start_t, fidelty, cur_pos);
					do_orbit(pln, end_t, fidelty, end_pos);

					double pln_start = get_V_phys(pln, cur_pos);
					double pln_end = get_V_phys(pln, end_pos);

					pln_start = ang_wrap(pln_start);

					double laps = floor((end_t - start_t) / pln.t_p) * M_2PI;

					if (pln_end < pln_start + laps)
						pln_end += ceil(((pln_start + laps) - pln_end) / M_2PI) * M_2PI;

					starts.push_back(pln_start);			//Get the current true anomaly of the body. 
					ends.push_back(pln_end);

					list.push_back(pln.self);
					diffs.push_back(pln_end - pln_start);
				}
			}

			vector<vector<pred_p>> mins(list.size());

			vector<double> last_times(list.size(), start_t - pow(10.0, -10.0));
			vector<int>	its(list.size(), 0);
			vector<short> dlts(list.size(), 1);
			vector<unsigned char> past_peak(list.size(), false);

			pred_p basic_queue_point(0, 0);
			vector<pred_p> basic_queue(10, basic_queue_point);

			vector<vector<pred_p>> logs_td(list.size(), basic_queue);
			basic_queue.clear();
			vector<vector<pred_p>> pred_lists(list.size(), basic_queue);

			int end_body = -1;
			double end_time = DBL_MAX;
			double end_V = NAN;
			bool ending = false;

			if (list.size())
			{
				precision = precision_base * ((end)-start) / M_2PI;
				if (precision < precision_base / 20)
					precision = precision_base / 20;


				for (int i2 = 0; i2 < precision; i2++)
				{
					double Vs = start + ang_scale(diff, i2, precision);

					sat_pos_old = sat_pos;
					sat_tim_old = sat_tim;

					sat_pos = get_pos_ang(Vs, sat);
					sat_tim = M_to_time(sat, get_M(Vs, sat.ecc), sat_tim_old, t_fid);

					vec_n sat_pos_diff = sat_pos - sat_pos_old;
					double sat_tim_diff = sat_tim - sat_tim_old;


					for (int i = 0; i < list.size(); i++)
					{
						body &pln = *list[i];
						double Vp = starts[i] + ang_scale(diffs[i], its[i], precision);
						double pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), last_times[i], t_fid);
						last_times[i] = pln_tim;
						vec_n pln_pos = get_pos_ang(Vp, pln);

						while (pln_tim < sat_tim_old && its[i] < precision)
						{
							its[i]++;
							Vp = starts[i] + ang_scale(diffs[i], its[i], precision);
							pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), last_times[i], t_fid);
							last_times[i] = pln_tim;
						}

						while (pln_tim < sat_tim && its[i] < precision)
						{
							double diff_part = (pln_tim - sat_tim_old) / sat_tim_diff;	//What percentage of the current time interval has passed in the moment we're observing.
							vec_n sat_pos_est = sat_pos_old + sat_pos_diff * diff_part;
							vec_n pln_pos = get_pos_ang(Vp, pln);
							double dist_est = vmag(sat_pos_est - pln_pos);
							short dlt_new;

							pred_p pln_td(dist_est, pln_tim);

							logs_td[i].push_back(pln_td);
							logs_td[i].erase(logs_td[i].begin());

							if (logs_td[i][0].dist > logs_td[i].back().dist)
								dlt_new = -1;
							else
								dlt_new = 1;

							if (dlts[i] != dlt_new)
							{
								pred_p queue_ext = logs_td[i][0];
								for (int i3 = 0; i3 < logs_td[i].size(); i3++)
									if (queue_ext.dist < logs_td[i][i3].dist == dlt_new > 0)
										queue_ext = logs_td[i][i3];

								queue_ext.dlt = dlt_new;

								pred_lists[i].push_back(queue_ext);

								dlts[i] = dlt_new;
							}
							if (pln_td.dist < pln.SOI && (pln_tim < end_time || !ending))
							{
								ending = true;
								end_body = i;
								end_time = pln_tim;
								end_V = Vs;
								break;
							}

							if (pln_td.dist < logs_td[i][0].dist)
								past_peak[i] = true;

							its[i]++;
							Vp = starts[i] + ang_scale(diffs[i], its[i], precision);
							pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), last_times[i], t_fid);
							last_times[i] = pln_tim;
						}
					}
				}
			}

			for (int i = 0; i < list.size(); i++)
				for (int i2 = 0; i2 < pred_lists[i].size(); i2++)
					if (pred_lists[i][i2].dlt > 0)
						mins[i].push_back(pred_lists[i][i2]);


			list_out = list;
			will_expire = ending;
			new_parent = end_body;
			expire_time = end_time;
			V_exp = end_V;

			mins_out = mins;
		}
	};

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

	void reset_expiry(body &sat)
	{
		sat.exit = false;
		sat.V_exit = 0;
	}

	void set_expiry_regular(body &sat, double apo = NAN)
	{
		if (isnan(apo))
			apo = get_r(M_PI, sat);

		body &parent = *sat.parent;
		if (sat.shape || apo > parent.SOI)
		{
			double V_max = get_V_r(parent.SOI, sat);
			double M_max = get_M(V_max, sat.ecc);
			sat.exit = M_to_time(sat, M_max, sat.epoch);
			sat.exiting = 1;
			sat.V_exit = V_max;

		}

		if (sat.inverse)
			sat.V_exit *= -1;
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

			rel_ang = pos_r.ang + M_PI2 - vel_r.ang;
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

		vector<double> r_vec = get_extremes(vel_r.mag, rel_ang, pos_r.mag, u);

		sat.ax_a = (r_vec[1] + r_vec[0]) / 2;
		sat.ecc = (r_vec[1] - r_vec[0]) / (r_vec[1] + r_vec[0]);
		sat.ax_b = sqrt_a(1 - sqr(sat.ecc)) * sat.ax_a;
		sat.shape = sat.ecc > 1;
		sat.SOI = sat.ax_a * pow(sat.u / u, 0.4);

		double tV = get_V_r(pos_r.mag, sat);	//The current true anomaly.
		if (rel_ang < 0)
			tV = ang_wrap(-tV);

		sat.normal = ang_wrap(pos_r.ang - tV, 2);
		sat.Mn = Mn_from_r(sat.ax_a, u);						//Store mean angular motion of the body, as per yet another equation dug up from obscure wikipedia stubs
		sat.En = sqr(vel_r.mag) / 2 - u / pos_r.mag;			//Store orbital Energy, velocital energy plus gravitational energy
		sat.Ar = vel_r.mag * pos_r.mag * cos(rel_ang);			//Store area swept per unit of time, altitude times velocity times cosine-of-the-angle

		if (sat.inverse)
			sat.Mn *= -1;

		sat.epoch = w_time - get_M(tV, sat.ecc) / sat.Mn;
		sat.t_p = M_2PI / abs(sat.Mn);
		sat.safe = w_time;
		sat.exiting = 0;

		set_expiry_regular(sat, r_vec[1]);
	}

	//End orbital algorithms

	//Begin graph algorithms

	vector<vec_n> make_tail(body sat, int subdiv, bool debug = false)
	{
		vector<vec_n> out;

		double V_start = ang_wrap(get_V_phys(sat), 2);
		double V_end;

		if (sat.expire())
			V_end = sat.V_exp();
		else
			V_end = V_start + M_2PI;

		V_end = nextafter(V_end, 0);

		double V_span = V_end - V_start;

		for (int i = 0; i <= subdiv; i++)
			out.push_back(get_pos_ang(V_start + ang_scale(V_span, i, subdiv), sat));

		return out;
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

	void set_tails_future(vector<body*> list, int subdiv = 100)
	{
		body plyr = *list.back();
		gen::tails_future.clear();

		body& pln_fut = *new body;
		body& plyr_fut = *new body;
		body& gramp = *new body;

		if (plyr.entering)
		{
			body pln = *list[pred::next_parent];

			double exp_time = plyr.entry;

			vec_n fut_plyr_pos, fut_pln_pos, fut_rel_pos;
			vec_n fut_plyr_vel, fut_pln_vel, fut_rel_vel;

			do_orbit(plyr, exp_time, 10, fut_plyr_pos, fut_plyr_vel);
			do_orbit(pln, exp_time, 10, fut_pln_pos, fut_pln_vel);

			game::target_pos = fut_pln_pos;
			game::player_pos = fut_plyr_pos;

			{
				pln_fut = pln;
				pln_fut.split();
				pln_fut.pos = fut_pln_pos;
				pln_fut.vel = fut_pln_vel;
			}

			{
				plyr_fut.pos = fut_plyr_pos;
				plyr_fut.vel = fut_plyr_vel;
				plyr_fut.u = 1;
				plyr_fut.parent = pln_fut.self;
			}

			make_orbit(plyr_fut, exp_time);
			set_expiry_regular(plyr_fut);

			gen::tails_future.push_back(make_tail(plyr_fut, subdiv));

			gramp = *plyr.parent;
			gramp.split();
		}
		else
		{
			plyr_fut = plyr;
			pln_fut = *plyr.parent;
			gramp = *pln_fut.parent;

			gramp.split();
			plyr_fut.split();
			pln_fut.split();
		}

		gramp.pos = vec_n();
		gramp.vel = vec_n();

		if (plyr_fut.exiting)
		{
			double exit = plyr_fut.exit;

			do_orbit(plyr_fut,	exit, 10, plyr_fut.pos,	plyr_fut.vel);
			do_orbit(pln_fut,	exit, 10, pln_fut.pos,	pln_fut.vel);

			plyr_fut.pos += pln_fut.pos;
			plyr_fut.vel += pln_fut.vel;

			plyr_fut.parent = gramp.self;

			make_orbit(plyr_fut, exit);
			set_expiry_regular(plyr_fut);

			gen::tails_future.push_back(make_tail(plyr_fut, subdiv*10));
		}

		delete plyr_fut.self;
		delete pln_fut.self;
	}

	//End graph algorithms

	//Begin game functions

	double calc_dmg(double dist, double radius, double intensity)
	{
		if (dist > radius)
			return 0;

		double dist_part = dist / radius;
		double dmg_base = 1 - pow(dist_part, 0.35);

		if (isnan(dmg_base))
			dmg_base = 0;

		return dmg_base * intensity;
	}

	//End game functions

	void generate_goal(double radius, int par_id = 0)
	{
		double sat_r = (*gen::bodies[par_id]).size;
		vec_r goal_coords;

		goal_coords.mag = rand_part() * (radius - sat_r) + sat_r;
		goal_coords.ang = to_rad(rand_part(), 1);

		game::goal_coords = goal_coords;
		game::goal_parent = par_id;
		game::goal_size = goal_coords.mag / 10;
	}

	bool compare_goal()
	{
		vec_n goal_pos = game::goal_coords + (*gen::bodies[game::goal_parent]).pos;

		body &plyr = *gen::bodies.back();

		double goal_dist = vmag(goal_pos - plyr.pos);

		if (goal_dist < game::goal_size)
			return true;

		return false;
	}

	void run_health()
	{
		using gen::kesses;
		double dmg_tot = 0;
		body plyr = *gen::bodies.back();

		for (int i = 0; i < kesses.size(); i++)
		{
			body kess = *kesses[i];
			double dist = vmag(kess.pos - plyr.pos);

			double dmg = calc_dmg(dist, game::dmg_rad, game::dmg_int);

			dmg_tot += dmg;
		}

		{
			body parent = *plyr.parent;

			double dist = vmag(plyr.pos - parent.pos) - parent.size;
			double atmosphere = parent.size * 0.05;

			if (dist < 0)
				game::health = 0;

			double dmg = calc_dmg(dist, atmosphere, game::dmg_int);
			dmg_tot += dmg;
		}

		game::health -= dmg_tot;
		if (game::health < 0)
			game::health = 0;
	}

	void add_kesslers()
	{
		body &plyr = *gen::bodies.back();
		body &parent = *plyr.parent;

		for (int i = 0; gen::kesses.size() < 5; i++)
		{
			body &kess = *new body;

			vec_r kess_pos;
			vec_r kess_vel;
			kess_pos.mag = parent.SOI * rand_part();
			kess_vel.mag = get_esc_vel(parent, kess_pos.mag) * rand_part() * 0.9;

			kess_pos.ang = rand_part() * M_2PI;
			kess_vel.ang = kess_pos.ang + M_PI2;

			kess.pos = kess_pos + parent.pos;
			kess.vel = kess_vel + parent.vel;
			kess.parent = parent.self;
			kess.isKess = true;
			kess.u = 0.01;
			kess.size = game::dmg_rad;

			make_orbit(kess, plyr.epoch);

			if (parent.size < get_r(0, kess) and get_r(M_PI, kess) < parent.SOI)
			{
				gen::kesses.push_back(&kess);
				gen::kesses_size.push_back(kess.size);
			}
			else
				delete kess.self;
		}
	}

	void kill_kesslers()
	{
		using gen::kesses;
		gen::kesses_size.clear();

		while (kesses.size())
		{
			body* sat_p = kesses.back();
			kesses.pop_back();
			delete sat_p;
		}

		std::cout << "Kill\n";
	}

	void do_kess_tick(vector<body*> list, vector<body*> list_pln, double w_time)
	{
		vector<vec_n> out;
		bool expired = false;

		vector<vec_n> pos_buffer;
		vector<vec_n> vel_buffer;

		if (!list.size())
			return;

		body& parent = *(*list_pln.back()).parent;

		for (int i = 0; i < list.size(); i++)
		{
			body &sat = *list[i];
			vec_n pos;
			vec_n vel;

			do_orbit(sat, w_time, 20, pos, vel);

			sat.pos = pos + parent.pos;
			sat.vel = vel + parent.vel;
			sat.t_l = w_time;

			out.push_back(sat.pos);
		}

		game::kesses_pos = out;
	}

	void do_phys_tick(vector<body*> bodies, double w_time, bool phys_mode, vec_n thrust = vec_n())
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
				do_orbit(sat, w_time, 10, pos, vel);
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

			out.push_back(sat.pos);
		}

		body &plyr = *bodies.back();


		{
			body &new_parent = *get_parent(plyr, bodies);
			if (new_parent.self != plyr.parent)
				expired = true;
			plyr.parent = &new_parent;
		}

		if (phys_mode || expired || plyr.expire())
		{
			if (expired || phys_mode)
			{
				make_orbit(plyr, w_time);
				pred::invalid = true;
			}

			if (expired || phys_mode || plyr.safe < plyr.t_l)
				plyr.safe = plyr.t_l;

			if (expired)
			{
				plyr.exiting = false;
				plyr.entering = false;
				set_expiry_regular(plyr);

				
			}
				
		}

		if (expired)
		{
			kill_kesslers();

			if (!(*plyr.parent).isSun)
				add_kesslers();
		}

		gen::bodies_pos = out;
	}

	void do_game_tick()
	{
		using namespace input;

		if (!keyboard.isPressed(key_state::keys::RShift))
			game::zoom_mem *= pow(1.1, keyboard.scroll);
		else
			gen::d_time_fact *= pow(1.1, keyboard.scroll);

		rock::spin += rock::gyro * gen::w_time_diff * keyboard.isPressed(key_state::keys::Left);
		rock::spin -= rock::gyro * gen::w_time_diff * keyboard.isPressed(key_state::keys::Right);

		rock::rot += gen::w_time_diff * rock::spin;

		rock::rot = ang_wrap(rock::rot, 2);

		rock::thrust += keyboard.isPressed(key_state::keys::Up);
		rock::thrust -= keyboard.isPressed(key_state::keys::Down);

		if (keyboard.isPressed(key_state::keys::Numpad1))
		{
			rock::spin = 0;
		}
		if (keyboard.isPressed(key_state::keys::Numpad3))
		{
			gen::d_time_fact = 100;
		}
		if (keyboard.isPressed(key_state::keys::Numpad6))
		{
			gen::d_time_fact = 1;
		}
		if (keyboard.isPressed(key_state::keys::Numpad2))
		{
			rock::thrust = 0;
		}
		if (keyboard.isPressed(key_state::keys::Numpad0))
		{
			std::cin >> game::target;
		}
		if (keyboard.isPressed(key_state::keys::Numpad4))
		{
			std::cin >> game::focus;
		}

		double thrust_actual = 0;
		vec_r accel;

		if (rock::thrust < 0)
			rock::thrust = 0;

		if (rock::fuel > 0)
			thrust_actual = rock::thrust * rock::eng_F;

		rock::fuel -= thrust_actual * rock::f_p_t;
		if (rock::fuel < 0)
			rock::fuel = 0;

		accel.mag = thrust_actual / (rock::mass + rock::fuel);
		accel.ang = rock::rot;

		namespace ws = shared::world_state;
		ws::focus = 6;
		ws::zoom = game::zoom_mem;
		ws::player_rotation = rock::rot;

		rock::eng_thrust = vec_to_pos(rock::rot, thrust_actual / (rock::mass + rock::fuel) );
		rock::eng_mode = vmag(rock::eng_thrust);
	}

	void do_predict()
	{
		using gen::bodies;
		using namespace pred;

		body &plyr = *bodies.back();

		if (shared::r_time < gen::last_predict + shared::cps * 0.1) return;
		if (plyr.inverse) return;

		get_expiry data(plyr, gen::bodies);

		vector<body*> co_sats = data.list_out;

		vector<vector<pred_p>> pairs = data.mins_out;

		vector<int> cnv;	//Table for translating the prediction function's indexes to regular bodies indexes.

		for (int i = 0; i < pairs.size(); i++)
			for (int i2 = 0; i2 < bodies.size(); i2++)
				if (bodies[i2] == co_sats[i])
					cnv.push_back(i2);

		bool expiring = data.will_expire;
		int new_parent = data.new_parent;
		double expire_time = data.expire_time;

		if (expiring)
			std::cout << "New parent: " << (*co_sats[new_parent]).name << "\n";

		for (int i = 0; i < pairs.size(); i++)
		{
			int ip = cnv[i];
			body &pln = *bodies[ip];

			if (expiring && i == new_parent) continue;

			for (int i2 = 0; i2 < pairs[i].size(); i2++)
			{
				vector<double> dists(3);

				double t_part = (pln.t_p / 1000);

				for (int i3 = 0; i3 < dists.size(); i3++)
				{
					vec_n true_pos_plyr, true_pos_pln;

					double tdiff = t_part * (i3 - 1);

					do_orbit(plyr, pairs[i][i2].time + tdiff, 10, true_pos_plyr);
					do_orbit(pln, pairs[i][i2].time + tdiff, 10, true_pos_pln);

					dists[i3] = vmag(true_pos_plyr - true_pos_pln);
				}

				if (dists[1] > dists[0] && dists[1] > dists[2])	//If the point is literally a local maximum, be rid of it.
				{
					pop_any(pairs[i], i2);
					i2--;
					continue;
				}

				int orbits_done = 0;
				int budget = 100;

				for (int i3 = 2; dists[1] > dists[2]; i3++)
				{
					vec_n true_pos_plyr, true_pos_pln;

					double tdiff = t_part * i3;

					do_orbit(plyr, pairs[i][i2].time + tdiff, 15, true_pos_plyr);
					do_orbit(pln, pairs[i][i2].time + tdiff, 15, true_pos_pln);

					double true_dist = vmag(true_pos_plyr - true_pos_pln);

					pop_any(dists, 0);
					dists.push_back(true_dist);

					orbits_done++;

					if (orbits_done > budget)
						break;
				}

				for (int i3 = 2; dists[1] > dists[0]; i3++)
				{
					vec_n true_pos_plyr, true_pos_pln;

					double tdiff = t_part * -i3;

					do_orbit(plyr, pairs[i][i2].time + tdiff, 15, true_pos_plyr);
					do_orbit(pln, pairs[i][i2].time + tdiff, 15, true_pos_pln);

					double true_dist = vmag(true_pos_plyr - true_pos_pln);

					push_any(dists, 0, true_dist);
					dists.pop_back();

					orbits_done++;

					if (orbits_done > budget)
						break;
				}

				pairs[i][i2].dist = dists[1];

				if (dists[1] > dists[0] || dists[1] > dists[2])	//If the point is still not a local minimum, be rid of it.
				{
					pop_any(pairs[i], i2);
					i2--;
				}

			}
		}

		plyr.safe = plyr.t_l;

		if (expiring)
		{
			vec_n true_pos_plyr, true_pos_pln;
			body &pln = *bodies[cnv[new_parent]];

			do_orbit(plyr, expire_time, 10, true_pos_plyr);
			do_orbit(pln, expire_time, 10, true_pos_pln);

			double true_dist = vmag(true_pos_plyr - true_pos_pln);

			if (true_dist > pln.SOI)
				expiring = false;
		}

		if (expiring)
		{
			plyr.entering = true;

			plyr.V_ent = data.V_exp;

			plyr.entry = expire_time;

			if (plyr.shape ? expire_time <= plyr.expire() : expire_time <= plyr.t_l + 2 * plyr.t_p)
				plyr.entry = expire_time;

			pred::next_parent = cnv[data.new_parent];
			pred::expire_count = 50;
			pred::end = plyr.entry;
		}
		else if (pred::expire_count > 0)
		{
			pred::expire_count--;
		}
		else
		{
			plyr.entering = false;
		}

		for (int i = 0; i < co_sats.size(); i++)
		{
			if (co_sats[i] == gen::bodies[game::target])
			{
				if (!pairs[i].size()) break;

				if (pairs[i][0].dist == DBL_MAX) break;

				game::min_dist = pairs[i][0].dist;
				game::min_time = pairs[i][0].time;

				break;
			}
		}
		gen::last_predict = shared::r_time;
	}

	void run_goal()
	{
		if (compare_goal() || game::goal_count < 0)
		{
			int par_id = rand() % gen::bodies.size();
			double radius = (*gen::bodies[par_id]).SOI * 0.8;

			if (par_id == 0)
			{
				double radius = 0;

				for (int i = 0; i < gen::bodies.size(); i++)
				{
					double r_max = get_r(M_PI, *gen::bodies[i]);

					if (r_max > radius)
						radius = r_max;
				}
				radius *= 0.8;
			}

			generate_goal(radius, par_id);

			game::goal_id = par_id;
			game::goal_count++;
		}
	}

	void handle_flushback(int flushback)
	{
		if (shared::game_state == 0)
		{
			if (flushback == 1)			//Play game
				shared::game_state = 1;
			if (flushback == 2)			//View options
				shared::game_state = 3;
			if (flushback == 3)			//View credits
				shared::game_state = 5;
			if (flushback == 4)			//End game
				shared::game_state = -1;
		}
	}

	void send_world_state()
	{
		using namespace shared::world_state;

		world_time = gen::w_time;

		health = game::health;
		health_max = game::health_max;
		target = game::target;

		score = pow(10, 9) * sqr(game::goal_count) / game::eng_secs;
		goal_count = game::goal_count;
		eng_secs = game::eng_secs;

		bodies = gen::bodies_pos;
		focus = game::focus;
		paths = gen::tails_out;

		goal_pos = game::goal_coords;
		goal_rad = game::goal_size;
		goal_parent = game::goal_id;
		goal_parent_2 = find_in(gen::bodies, (*gen::bodies[game::goal_id]).parent);

		sizes = gen::bodies_size;
		kesses = game::kesses_pos;
		sizes_kess = gen::kesses_size;

		cepting = (*gen::bodies.back()).entering;
		intercept = pred::next_parent;
		cept_time = pred::end;

		target_close[0] = game::target_pos;
		player_close[0] = game::player_pos;
		target_min[0] = game::min_dist;
		target_time[0] = game::min_time;
	}

#ifdef RENDER_DEBUG_INSTALLED
	bool emode = 0;
	vec_n thrust_debug;
#endif

	void run_engine()
	{
		handle_flushback(input::flush_back::men_cmd);

		gen::w_time_last = gen::w_time;
		double diff_time = (shared::r_time - shared::l_time) / 100.0 / shared::cps / gen::d_time_fact;
		if (diff_time < 0.01)
			gen::w_time_diff = diff_time;
		else
			gen::w_time_diff = 0.01;

		if (shared::game_state == 1)
			gen::w_time += gen::w_time_diff;	
		else
			return;

#ifdef RENDER_DEBUG_INSTALLED
		thrust_debug = (*gen::bodies.back()).vel;
#endif

		do_game_tick();
		run_health();

		do_phys_tick(gen::bodies, gen::w_time, rock::eng_mode, rock::eng_thrust);
		do_kess_tick(gen::kesses, gen::bodies, gen::w_time);

		body &plyr = *gen::bodies.back();
		gen::tails[gen::tails.size() - 1] = make_tail(plyr, 1000);

		do_predict();

		run_goal();

		set_tails_future(gen::bodies, 100);

		game::cur_dist = vmag((*gen::bodies.back()).pos - (*gen::bodies[game::target]).pos);
		vector<vector<vec_n>> tails_out = gen::tails;

		if (gen::tails_future.size())
			tails_out.push_back(gen::tails_future[0]);

		for (int i = 0; i < tails_out.size(); i++)
		{
			vec_n par_pos = (*gen::bodies[pred::next_parent]).pos;

			if (i < gen::bodies.size())
				par_pos = (*(*gen::bodies[i]).parent).pos;

			for (int i2 = 0; i2 < tails_out[i].size(); i2++)
				tails_out[i][i2] += par_pos;
		}

		gen::tails_out = tails_out;

#ifdef RENDER_DEBUG_INSTALLED
		thrust_debug = (*gen::bodies.back()).vel - thrust_debug;
#endif
		send_world_state();
	}

	vector<body*> sort_bodies(vector<body*> unsorted)	//Bubble-ish sorting of the bodies by mass, largest to smallest.
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
					if ((*unsorted[i]).u > (*body_p).u)
					{
						body_p = unsorted[i];
						body_i = i;
					}
				sorted.push_back(body_p);
				pop_any(unsorted, body_i);
			}

			//Name the largest body as 'sun';
			(*sorted[0]).isSun = true;
			(*sorted[0]).SOI = DBL_MAX;
			//To begin with, assume all sorted bodies orbit the sun.
			for (int i = 0; i < sorted.size(); i++)
				(*sorted[i]).parent = sorted[0];

			for (int i = 0; i < sorted.size(); i++)
			{
				(*sorted[i]).pos -= (*sorted[0]).pos;
				(*sorted[i]).vel -= (*sorted[0]).vel;
			}
		}

		return sorted;
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
		moon.vel.y = -46000;
		moon.u = pow(10, 9);
		moon.name = "moon";

		body& planet = *new body;
		planet.pos.x = -25000;
		planet.pos.y = 0;
		planet.vel.x = 0;
		planet.vel.y = -40000;
		planet.u = pow(10, 10);
		planet.name = "planet";

		body& pluto = *new body;
		pluto.pos.x = -30000;
		pluto.pos.y = -0;
		pluto.vel.x = 0;
		pluto.vel.y = -40000;
		pluto.u = pow(10, 10);
		pluto.name = "pluto";

		body& dune = *new body;
		dune.pos.x = -35000;
		dune.pos.y = 0;
		dune.vel.x = 0;
		dune.vel.y = -40000;
		dune.u = pow(10, 10);
		dune.name = "dune";

		body& laythe = *new body;
		dune.pos.x = -35000;
		dune.pos.y = 5000;
		dune.vel.x = 0;
		dune.vel.y = -40000;
		dune.u = pow(10, 10);
		dune.name = "layhte";

		body& yavin = *new body;
		yavin.pos.x = -40000;
		yavin.pos.y = 0;
		yavin.vel.x = 0;
		yavin.vel.y = -40000;
		yavin.u = pow(10, 10);
		yavin.name = "yavin";

		body& plyr = *new body;
		plyr.pos.x = -40200;
		plyr.pos.y = 0;
		plyr.vel.x = 0;
		plyr.vel.y = -50000;
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

		gen::bodies = sort_bodies(unsorted);

		for (int i2 = 0; i2 < 5; i2++)
			for (int i = 0; i < gen::bodies.size(); i++)
			{
				make_orbit(*gen::bodies[i], 0);								//Get everything's orbit around the sun
				(*gen::bodies[i]).parent = get_parent(*gen::bodies[i], gen::bodies);	//Check if the newly-determined SOIs make them into moons
			}

		for (int i = 0; i < gen::bodies.size(); i++)
		{
			body &sat = *gen::bodies[i];

			sat.size = cbrt(sat.u) / 40;
			gen::bodies_size.push_back(sat.size);
		}

		for (int i = 0; i < gen::bodies.size(); i++)
			shared::world_state::names.push_back((*gen::bodies[i]).name);

		gen::tails = get_tails_basic(gen::bodies, 1000);

		rock::fuel = 1000;
		rock::mass = 1000000000;
		rock::f_p_t = 0.0;
		rock::gyro = 100000000;
		rock::eng_F = 1000000000000;

		game::health = 100;
		game::health_max = game::health;

		game::dmg_rad = 100;
		game::dmg_int = 0.001;

		game::focus = gen::bodies.size() - 1;
	}

	void engine_init()
	{
		phys_init();
		add_kesslers();
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
	bool nonsense = debug_font.loadFromFile("courier.ttf");

	bool window_is_clear = false;

	short count = 0;


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
		FPS2.setFillColor(Color(150, 150, 150));
		if (nonsense)
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
		pos = scale_single(pos, origo, zoom);

		CircleShape plyr(10, 3);
		plyr.setFillColor(sf::Color::Green);
		plyr.setOrigin(5, 5);
		plyr.setPosition(pos);
		plyr.setScale(vec_n(1, 2));
		plyr.setRotation(-rot / phys::M_2PI * 360.0 + 90.0);

		window2.draw(plyr);
	}

	void render_kesslers(vector<vec_n> list, vec_n origo, double zoom)
	{
		list = handle_scale(list, origo, zoom);

		for (int i = 0; i < list.size(); i++)
		{
			CircleShape kess(2, 4);
			kess.setFillColor(sf::Color::Red);
			kess.setOrigin(1, 1);
			kess.setPosition(list[i]);
			
			
			window2.draw(kess);
		}
	}

	void render_goal(vec_n pos, vec_n origo, double zoom)
	{
		pos += (*phys::gen::bodies[phys::game::goal_parent]).pos;

		pos = scale_single(pos, origo, zoom);

		double size = phys::game::goal_size * zoom;

		CircleShape kess(size);
		kess.setFillColor(sf::Color::Magenta);
		kess.setOrigin(size / 2, size / 2);
		kess.setPosition(pos);

		window2.draw(kess);
	}

	void render_coll(vec_n origo, double zoom)
	{
		
		{
			double size = 6;
			shared::vec_n pos = phys::game::target_pos;

			pos = scale_single(pos, origo, zoom);

			CircleShape kess(size);
			kess.setFillColor(sf::Color::Blue);
			kess.setOrigin(size / 2, size / 2);
			kess.setPosition(pos);

			window2.draw(kess);
		}
		{
			double size = 3;
			shared::vec_n pos = phys::game::player_pos;

			pos = scale_single(pos, origo, zoom);

			CircleShape kess(size);
			kess.setFillColor(sf::Color::Red);
			kess.setOrigin(size / 2, size / 2);
			kess.setPosition(pos);

			window2.draw(kess);
		}
	}

	void render_all()
	{
		if (shared::game_state != 1)
		{
			return;
		}
		else if (count < 10)
		{
			count++;
			return;
		}


		namespace ws = shared::world_state;

		window_is_clear = true;
		render_lines(ws::paths, ws::bodies[ws::focus], ws::zoom);

		double gm_rot = phys::rock::rot;// phys::ang_scale(4, phys::rock::rot, phys::M_2PI);
		double ph_rot = phys::atan2(phys::thrust_debug);// phys::ang_scale(4, phys::atan2(phys::thrust_debug), phys::M_2PI);
		double df_rot = gm_rot - ph_rot;

		while (df_rot < 0)
			df_rot += 4;
		while (df_rot > 4)
			df_rot -= 4;

		std::vector<std::string> texts;
		texts.push_back("Engine: " + std::to_string(phys::rock::thrust));
		texts.push_back("GM_rot: " + std::to_string(gm_rot));
		texts.push_back("PH_rot: " + std::to_string(ph_rot));
		texts.push_back("DF_rot: " + std::to_string(df_rot));

		texts.push_back("Eccent: " + std::to_string((*phys::gen::bodies.back()).ecc));

		texts.push_back("PrDist: " + std::to_string(phys::game::min_dist));
		texts.push_back("PrTime: " + std::to_string(phys::game::min_time - phys::gen::w_time));
		texts.push_back("AcDist: " + std::to_string(phys::game::cur_dist));
		texts.push_back("Expire: " + std::to_string((*phys::gen::bodies.back()).expire()));

		render_texts(texts);
		render_kesslers(shared::world_state::kesses, ws::bodies[ws::focus], ws::zoom);
		render_player(ws::player_rotation, ws::bodies.back(), ws::bodies[ws::focus], ws::zoom);

		render_goal(phys::game::goal_coords, ws::bodies[ws::focus], ws::zoom);
		render_coll(ws::bodies[ws::focus], ws::zoom);

		//window2.display();
	}
}
#endif // RENDER_DEBUG_INSTALLED