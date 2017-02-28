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
	};

	struct vec_n
	{
		double x = 0;
		double y = 0;

		vec_n(double in_x = 0, double in_y =0)
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

		double player_rotation; //Angle in radians, sweeping counter-clockwise.
	};

	world_state screen_state;
}

#define RENDER_DEBUG_INSTALLED true
#ifdef RENDER_DEBUG_INSTALLED
namespace render_debug			//To be removed once the neccesary render_tools functions are implemented
{
	using sf::Vertex;
	using sf::Color;
	using sf::Vector2f;
	using std::vector;
	using shared::vec_n;
	using shared::window2;

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

	void render_all(shared::world_state in)
	{
		window2.clear();
		window_is_clear = true;
		render_lines(in.paths, in.bodies[in.focus], in.zoom);
	}
}
#endif // RENDER_DEBUG_INSTALLED



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
		double expiry;	//The last point in time where the current orbital parameters are guaranteed accurate
		short expire;	//Whether or not there is a predicted expiry. 0 = no confirmed expiry, 1 = expiry by leaving SOI, 2 = expiry by entering SOI

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

	body sun, moon, planet, pluto, yavin, dune, plyr;
	vector<body*> bodies;
	vector<body> bodies_actual;
	vector<vector<vec_n>> tails;
	double zoom_mem = 0.01;

	struct dual_val
	{
		long double one = 0;
		long double two = 0;
	};

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

	double M_to_time(body sat, double M, double w_time_min = NAN)
	{
		double d_time = M / sat.Mn;
		double w_time = d_time + sat.epoch;

		if (sat.shape == 0 && !isnan(w_time_min))
		{
			double r_time = M_2PI / sat.Mn;
			while (w_time < w_time_min)
				w_time += r_time;
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

	vec_n get_pos_ang(double V, body sat)
	{
		vec_r point;
		point.ang = -V;
		point.mag = get_r(point.ang, sat);
		point.ang += (sat).normal;

		return point;
	}

	//End orbital equation functions

	//Begin orbital algorithms
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
			if (get_M_E(mid, sat.ecc) > part)
				ceil = mid;
			else
				floor = mid;

			mid = (ceil + floor) / 2;
		}

		return get_V(mid, sat.ecc);
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
			V = M_2PI - V;

		radius = get_r(V, sat);

		pos_out = get_pos_ang(V, sat);

		double vel_mag = sqrt(2 * (sat.En + parent.u / radius));
		double vel_ang = acos(sat.Ar / (radius * vel_mag));
		if (V < M_PI)
			vel_ang = -vel_ang;

		if (sat.inverse)
			vel_mag = -vel_mag;

		vel_ang = ang_wrap(-V - vel_ang + sat.normal + M_PI2);

		vel_out.x = -vel_mag * cos(vel_ang);
		vel_out.y = -vel_mag * sin(vel_ang);
	}

	void get_expiry(body &sat, vector<body*> list, double w_time, vector<body*> &list_out, vector<vector<double>> &mins_out, bool &will_expire, int &new_parent, double &expire_time)
	{
		vector<body*> list_in = list;
		vector<double> starts;			//Current V-values for bodies
		vector<double> ends;			//V-values for bodies at the end of the relevant time span
		vector<int>	its;				//The iterators for the bodies
		list = vector<body*>(0);
		int precision = 1000;

		vector<vector<double>> pairs_dist;
		vector<vector<double>> pairs_time;

		vector<double> min_dists;
		vector<double> min_times;

		for (int i = 0; i < list_in.size(); i++)
		{
			body &pln = *list_in[i];
			if (pln.parent == sat.parent && pln.self != sat.self && !pln.isSun)
			{
				vec_n pln_pos = pln.pos - (*pln.parent).pos;
				vec_n end_pos;
				list.push_back(pln.self);
				starts.push_back(pln.normal - atan2(pln_pos));		//Get the current true anomaly of the body. atan2(pos) = normal - V --> V = normal - atan2(pos)

				do_orbit(pln, sat.expiry, end_pos);
				ends.push_back(pln.normal - atan2(end_pos));
				its.push_back(0);

				pairs_dist.push_back(vector<double>());
				pairs_time.push_back(vector<double>());
				min_times.push_back(DBL_MAX);
				min_dists.push_back(DBL_MAX);
			}
		}
		double start = sat.normal - atan2(sat.pos - (*sat.parent).pos);
		vec_n end_pos_plyr;
		do_orbit(sat, sat.expiry, end_pos_plyr);
		double end = sat.normal - atan2(end_pos_plyr);
		
		vec_n sat_pos_old;
		double sat_tim_old;

		for (int i2 = 0; i2 < precision; i2++)
		{
			int i_sat = 0;
			int i_pln = 0;
			double Vs = to_rad(i2, precision, end-start);

			vec_n sat_pos = get_pos_ang(start + Vs, sat);
			double sat_tim = M_to_time(sat, get_M(Vs, sat.ecc), sat.t_l);

			vec_n sat_pos_diff = sat_pos - sat_pos_old;
			double sat_tim_diff = sat_tim - sat_tim_old;

			for (int i = 0; i < list.size(); i++)
			{
				body &pln = *list[i];
				double Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
				double pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);
				vec_n pln_pos = get_pos_ang(Vp, pln);

				while (pln_tim < sat_tim_old)
				{
					its[i]++;
					Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
					pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);
				}
				while (pln_tim < sat_tim)
				{
					double diff_part = (pln_tim - sat_tim_old) / sat_tim_diff;
					vec_n sat_pos_est = sat_pos_old + sat_pos_diff * diff_part;

					double dist_est = vmag(sat_pos_est - pln_pos);
					pairs_dist[i].push_back(dist_est);
					pairs_time[i].push_back(pln_tim);

					its[i]++;
					Vp = starts[i] + to_rad(its[i], precision, ends[i] - starts[i]);
					pln_tim = M_to_time(pln, get_M(Vp, pln.ecc), pln.t_l);
				}
			}
			sat_pos_old = sat_pos;
			sat_tim_old = sat_tim;
		}


		body *parent_next_p = (*sat.parent).parent;
		double end_time = sat.expiry;
		double end_body = -1;
		bool ending = false;

		for (int i = 0; i < list.size(); i++)
		{
			for (int i2 = 0; i2 < pairs_time.size(); i2++)
			{
				if (pairs_dist[i][i2] < min_dists[i])
				{
					min_dists[i] = pairs_dist[i][i2];
					min_times[i] = pairs_time[i][i2];
				}
				if (pairs_dist[i][i2] < (*list[i]).SOI)
					if (pairs_time[i][i2] < end_time)
					{
						ending = true;
						end_time = pairs_time[i][i2];
						end_body = i;
						break;
					}
			}
		}

		list_out = list;
		will_expire = ending;
		expire_time = end_time;
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

		if (sat.shape || r_vec.two > parent.SOI)
		{
			double V_max = get_V_r(parent.SOI, sat);
			double M_max = get_M(V_max, sat.ecc);
			sat.expiry = M_to_time(sat, M_max, sat.epoch);
		}
	}

	//End orbital algorithms

	vector<vec_n> do_phys_tick(vector<body*> bodies, double w_time, bool phys_mode = 0)
	{
		vector<vec_n> out;
		bool emodelast = false;


		for (int i = 0; i < bodies.size(); i++)
		{
			body &sat = *bodies[i];
			vec_n pos;
			vec_n vel;
			do_orbit(sat, w_time, pos, vel);

			sat.pos = pos + (*sat.parent).pos;
			sat.vel = vel + (*sat.parent).vel;
			sat.t_l = w_time;
			if (!sat.expire)
				sat.expiry = sat.t_l + sat.t_p;

			out.push_back(sat.pos);
		}

		//get_expiry(*bodies[bodies.size() - 1], bodies, w_time);
		return out;
	}

	void run_engine()
	{
		using shared::screen_state;

		zoom_mem *= pow(1.1, input::keyboard.scroll);

		screen_state.zoom = zoom_mem;
		screen_state.focus = 4;
		screen_state.bodies = do_phys_tick(bodies, shared::r_time * 0.0001, 0);

		vector<vector<vec_n>> tails_out;

		for (int i = 0; i < tails.size(); i++)
		{
			vector<vec_n> temp;
			vec_n par_pos = (*(*bodies[i]).parent).pos;
			tails_out.push_back(temp);
			for (int i2 = 0; i2 < tails.size(); i2++)
			{
				tails_out[i].push_back(tails[i][i2] + par_pos);
			}
		}

		screen_state.paths = tails;
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
			vector<vec_n> planet;
			paths.push_back(planet);

			for (int i2 = 0; i2 <= subdiv; i2++)
				paths[i].push_back( get_pos_ang( to_rad(i2, subdiv), *list[i] ) );
		}
		return paths;
	}

	void phys_init()
	{
		sun.pos.x = 0;
		sun.pos.y = 0;
		sun.vel.x = 0;
		sun.vel.y = 0;
		sun.u = pow(10, 14);
		sun.name = "Sun";

		moon.pos.x = -25300;
		moon.pos.y = 0;
		moon.vel.x = 0;
		moon.vel.y = 40000;
		moon.u = pow(10, 8);
		moon.name = "moon";

		planet.pos.x = -25000;
		planet.pos.y = 0;
		planet.vel.x = 0;
		planet.vel.y = 39900;
		planet.u = pow(10, 10);
		planet.name = "planet";

		pluto.pos.x = -30000;
		pluto.pos.y = -0;
		pluto.vel.x = 0;
		pluto.vel.y = 40000;
		pluto.u = pow(10, 10);
		pluto.name = "pluto";

		dune.pos.x = -35000;
		dune.pos.y = 0;
		dune.vel.x = 0;
		dune.vel.y = 40000;
		dune.u = pow(10, 10);
		dune.name = "dune";

		yavin.pos.x = -40000;
		yavin.pos.y = 0;
		yavin.vel.x = 0;
		yavin.vel.y = 40000;
		yavin.u = pow(10, 10);
		yavin.name = "yavin";

		plyr.pos.x = -40200;
		plyr.pos.y = 0;
		plyr.vel.x = 0;
		plyr.vel.y = 43000;
		plyr.u = 1;
		plyr.name = "plyr";

		vector<body*> unsorted;

		unsorted.push_back(&sun);		//White
		unsorted.push_back(&moon);		//Blue
		unsorted.push_back(&planet);	//Yellow
		unsorted.push_back(&pluto);		//Green
		unsorted.push_back(&dune);		//Cyan
		unsorted.push_back(&yavin);		//Magenta
		unsorted.push_back(&plyr);		//Red

		bodies = sort_bodies(unsorted);

		for (int i2 = 0; i2 < 5; i2++)
			for (int i = 0; i < bodies.size(); i++)
			{
				make_orbit(*bodies[i], 0);								//Get everything's orbit around the sun
				(*bodies[i]).parent = get_parent(*bodies[i], bodies);	//Check if the newly-determined SOIs make them into moons
			}

		tails = get_tails_basic(bodies);
	}

	void engine_init()
	{
		phys_init();
	}
}