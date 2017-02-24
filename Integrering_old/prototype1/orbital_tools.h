#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <windows.h>
#define skip if (false)

using namespace std;


clock_t cur_time = clock();
vector <sf::CircleShape> tail;
int cps = CLOCKS_PER_SEC;
const long double M_PI = acos(-1.0L);
const long double M_2PI = M_PI * 2L;
const long double M_PI2 = M_PI / 2L;
const long double M_3PI2 = M_PI2 * 3L;

const long double M_E = exp(1.0);


sf::RenderWindow window2(sf::VideoMode(1080, 860), "Orbitals");

struct vec_n
{
	long double x = 0;
	long double y = 0;
	long double z = 0;
};

void draw_circle(double x, double y, bool beblue) //Render either of the planets
{

	sf::CircleShape planet;

	planet.setPosition(x - 5, y - 5); // Subtract five to compensate for circle size
	planet.setFillColor(beblue ? sf::Color::Blue : sf::Color::Red); // Blue for orbiting body, red for central
	planet.setRadius(10);


	window2.draw(planet);
}

void draw_tail(double x, double y, bool beblue, bool draw = false, bool add = true) //Render either of the tails
{
	sf::CircleShape holder = tail[0]; //Declare tail particle
	holder.setPosition(x, y); //Apply properties
	holder.setRadius(1.5); //Radius
	holder.setFillColor((beblue ? sf::Color::Red : sf::Color::Cyan)); //Green for orbiting body, cyan for central

	if (add)
	{
		tail.push_back(holder); //Store in main vector
		tail.erase(tail.begin()); //Erase oldest tail particle.
	}

	//cout << "ONE \n";

	if (draw) for (unsigned short i = 0; i < tail.size(); i++)
	{
		holder = tail[i];
		//holder.setPosition();
		window2.draw(holder); //Draw all tail particles
	}

}



struct body
{
	double ax_a;	//Semimajor axis
	double ax_b;	//Semiminor axis
	double ecc;		//Orbital eccentrictity
	double t_0;		//Starting timestamp
	double t_l;		//Timestamp at a given time
	double u;		//Gravitational parameter
	double normal;	//angle of the major axis
	double area;	//Orbital area
	double SOI;
	bool closed;
	bool inverse;
	vec_n vel;		//Velocity at a given time
	vec_n pos;		//Position at a given time

	double En;		//Orbital energy
	double Ar;		//Area swept per unit of time
	double Mn;		//Mean anomaly swept per unit of time

	body* parent;
};

double sqr(double in)
{
	return pow(in, 2);
}

double cot(double x)
{
	return tan(M_PI2 - x);
}

double acot(double x)
{
	return M_PI2 - atan(x);
}

vec_n t_to_screen(vec_n in)
{
	vec_n out;
	out.y = 500 - in.y;
	out.x = 500 + in.x;

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

double dist(double x, double y)
{
	return sqrt(sqr(x) + sqr(y));
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
	return acos(sat.ax_a * (1 - sqr(sat.ecc)) / (sat.ecc * r) - 1 / sat.ecc);
}

void fill_tail() //Generate one thousand tail particles, fill main vector.
{
	for (int i = 0; i < 10000; i++)
	{
		sf::CircleShape holder;
		holder.setPosition(0, 0);
		holder.setRadius(1.5);
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

double get_M(double tV, double ecc)
{
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

double get_area(double tV, body sat)
{
	return get_M(tV, sat.ecc) / (M_2PI)* sat.area;
}

double get_area_E(double E, body sat)
{
	return get_M_E(E, sat.ecc) / M_2PI * sat.area;
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

double get_hM(double V, double ecc)
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

vec_n rot_vec(vec_n in, double rot)
{
	double l = dist(in.x, in.y);
	double angle = atan2(in.y, in.x);

	angle += rot;

	in.y = l * sin(angle);
	in.x = l * cos(angle);

	return in;
}

double do_orbit_precise(double part, body sat, int precision)
{
	double ceil = M_2PI;
	double floor = 0;
	double mid = 0;

	for (int i = 0; i < precision; i++)
	{
		mid = (ceil + floor) / 2;

		if (get_M(mid, sat.ecc) > part)
			ceil = mid;
		else
			floor = mid;
	}

	return (ceil + floor) / 2;
	//return get_V((ceil + floor) / 2, sat.ecc);
}

double do_orbit_precise_H(double part, body sat, int precision)
{
	double ceil = M_PI;
	double floor = -M_PI;
	double mid = 0;

	for (int i = 0; i < precision; i++)
	{
		mid = (ceil + floor) / 2;

		if (get_hM(mid, sat.ecc) > part)
			ceil = mid;
		else
			floor = mid;
	}

	return (ceil + floor) / 2;
	//return get_V((ceil + floor) / 2, sat.ecc);
}

void do_orbit_phys(body &sat, double d_time)
{
	double m_time = d_time - sat.t_l;
	sat.t_l = d_time;

	body &parent = *sat.parent;
	vec_n pos_rel;
	vec_n force_vec;
	pos_rel.x = sat.pos.x - parent.pos.x;
	pos_rel.y = sat.pos.y - parent.pos.y;

	double angle = atan2(pos_rel.y, pos_rel.x);
	double pos_mag = dist(pos_rel.x, pos_rel.y);
	double F = -parent.u / sqr(pos_mag);

	force_vec.x = F * cos(angle);
	force_vec.y = F * sin(angle);

	sat.vel.x += force_vec.x * m_time;
	sat.vel.y += force_vec.y * m_time;

	sat.pos.x += sat.vel.x * m_time;
	sat.pos.y += sat.vel.y * m_time;
}

void do_orbit(body &sat, double d_time)
{
	sat.t_l = d_time;
	d_time = d_time - sat.t_0;

	body &parent = *sat.parent;

	//double part = to_rad(d_time / sat.t_p, 1);
	double part = ang_wrap(d_time * sat.Mn);
	double angle;
	double radius;

	if (sat.inverse)
		part = M_2PI - part;

	angle = do_orbit_precise(part, sat, 32);
	radius = get_r(angle, sat);

	vec_n out;
	sat.pos.y = sin(-angle) * radius;
	sat.pos.x = cos(-angle) * radius;

	sat.pos = rot_vec(sat.pos, sat.normal);

	sat.pos.x += parent.pos.x;
	sat.pos.y += parent.pos.y;

	double vel_mag = sqrt(2 * (sat.En + parent.u / radius));
	double vel_ang = acos(sat.Ar / (radius * vel_mag));
	if (angle < M_PI)
		vel_ang = -vel_ang;

	if (sat.inverse)
		vel_mag = -vel_mag;

	vel_ang = ang_wrap(-angle - vel_ang + sat.normal + M_PI2);

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
	ax_a = abs(ax_a);
	return sqrt(u / pow(ax_a, 3));
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


body make_orbit(body sat, double w_time)
{
	body &parent = *(sat.parent);
	double u = parent.u;

	vec_n pos, vel;
	double pos_mag, pos_ang, vel_mag, vel_ang, rel_ang;

	pos.x = sat.pos.x - parent.pos.x;
	pos.y = sat.pos.y - parent.pos.y;

	vel.x = sat.vel.x - parent.vel.x;
	vel.y = sat.vel.y - parent.vel.y;

	sat.inverse = false;
	do
	{
		pos_mag = dist(pos.x, pos.y);
		pos_ang = atan2(pos.y, pos.x);

		vel_mag = dist(vel.x, vel.y);
		vel_ang = atan2(vel.y, vel.x);

		rel_ang = M_PI2 + vel_ang - pos_ang;
		rel_ang = ang_wrap(rel_ang, 2);

		if (abs(rel_ang) > M_PI2)
		{
			vel.x = -vel.x;
			vel.y = -vel.y;
			sat.inverse = true;
		}
	} while (abs(rel_ang) > M_PI2);

	vector<dual_val> r_vec_vector = get_extremes(vel_mag, rel_ang, pos_mag, u);

	dual_val r_vec = r_vec_vector[0];
	dual_val r_vec_vel = r_vec_vector[1];

	double major = r_vec.two + r_vec.one;

	sat.ax_a = major / 2;
	sat.ecc = (r_vec.two - r_vec.one) / (r_vec.two + r_vec.one);
	sat.ax_b = sqrt(1 - pow(sat.ecc, 2)) * major / 2;
	sat.area = M_PI * sat.ax_a * sat.ax_b;
	sat.SOI = sat.ax_a * pow(parent.u / u, 0.4);

	double tV = get_V_r(pos_mag, sat);
	if (rel_ang < 0)
		tV = ang_wrap(-tV);

	sat.normal = pos_ang + tV;

	if (sat.inverse)
		tV = ang_wrap(-tV);
	sat.En = sqr(vel_mag) / 2 - u / pos_mag;			//Store orbital Energy, velocital energy plus gravitational energy
	sat.Ar = vel_mag * pos_mag * cos(rel_ang);			//Store area swept per unit of time, altitude times velocity times cosine-of-the-angle
	sat.Mn = Mn_from_r(sat.ax_a, u);					//Store mean angular motion of the body, as per yet another equation dug up from obscure wikipedia stubs
	sat.t_0 = w_time - get_M(tV, sat.ecc) / sat.Mn;

	return sat;
}

vector<vec_n> update_all_and_convert(vector<body*> bodies, double w_time)
{
	vector<vec_n> out(1);
	out[0] = (*bodies[0]).pos;

	for (int i = 1; i < bodies.size(); i++)
	{
		body &sat = *bodies[i];
		do_orbit(sat, w_time);

		out.push_back(sat.pos);
	}

	return out;
}

vector<body*> orbital_init(double w_time)
{
	body moon, planet, sun, pluto, dune, yavin;

	moon.pos.x = -240;
	moon.pos.y = -100;
	moon.vel.x = 1000;
	moon.vel.y = -4000;
	moon.u = pow(10, 10);
	moon.parent = &sun;

	planet.pos.x = -200;
	planet.pos.y = -100;
	planet.vel.x = -4000;
	planet.vel.y = 4000;
	planet.u = pow(10, 10);
	planet.parent = &moon;

	srand(time(0));

	sun.pos.x = 0;
	sun.pos.y = 0;
	sun.vel.x = 0;
	sun.vel.y = 0;
	sun.u = pow(10, 10);

	pluto.pos.x = -300;
	pluto.pos.y = -300;
	pluto.vel.x = 0;
	pluto.vel.y = 3000;
	pluto.u = pow(10, 10);
	pluto.parent = &sun;

	dune.pos.x = -300;
	dune.pos.y = 0;
	dune.vel.x = 1000;
	dune.vel.y = 3000;
	dune.u = pow(10, 10);
	dune.parent = &sun;

	yavin.pos.x = -400;
	yavin.pos.y = 0;
	yavin.vel.x = 1000;
	yavin.vel.y = -5000;
	yavin.u = pow(10, 10);
	yavin.parent = &dune;

	moon = make_orbit(moon, w_time);
	planet = make_orbit(planet, w_time);
	pluto = make_orbit(pluto, w_time);
	dune = make_orbit(dune, w_time);
	yavin = make_orbit(yavin, w_time);

	vector<body*> bodies;
	bodies.push_back(&sun);
	bodies.push_back(&moon);
	bodies.push_back(&planet);
	bodies.push_back(&pluto);
	bodies.push_back(&dune);
	bodies.push_back(&yavin);

	return bodies;
}