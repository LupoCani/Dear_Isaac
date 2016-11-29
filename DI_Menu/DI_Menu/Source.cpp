#include <iostream>
#include <time.h>
#include <SFML\Graphics.hpp>
#include <math.h>
#include <vector>
#include <string>
#include <windows.h>

using namespace std;

clock_t cur_time = clock();
vector <sf::CircleShape> tail;
int cps = CLOCKS_PER_SEC;


sf::RenderWindow window(sf::VideoMode(1080, 860), "Orbitals");

void draw_circle(double x, double y, bool beblue) //Render either of the planets
{

	sf::CircleShape planet;

	planet.setPosition(x - 5, y - 5); // Subtract five to compensate for circle size
	planet.setFillColor(beblue ? sf::Color::Blue : sf::Color::Red); // Blue for orbiting body, red for central
	planet.setRadius(10);


	window.draw(planet);
}

void draw_tail(double x, double y, bool beblue) //Render either of the tails
{
	sf::CircleShape holder; //Declare tail particle
	holder.setPosition(x, y); //Apply properties
	holder.setRadius(1.5); //Radius
	holder.setFillColor((beblue ? sf::Color::Green : sf::Color::Cyan)); //Green for orbiting body, cyan for central
																		//holder.time = time(0); //Start timeout timer
	tail.push_back(holder); //Store in main vector

	tail.erase(tail.begin()); //Erase oldest tail particle.

	for (unsigned short i = 0; i < tail.size(); i++)
		window.draw(tail[i]); //Draw all tail particles

}

/*
void draw_line(double x, double y, bool beblue)
{
sf::CircleShape holder;
holder.setPosition(x, y);
holder.setRadius(1.5);
holder.setFillColor((beblue ? sf::Color::Green : sf::Color::Cyan));
holder.time = time(0);
tail.push_back(holder);

tail.erase(tail.begin());

for (unsigned short i = 0; i < tail.size(); i++)
window.draw(tail[i]);

}*/

void fill_tail() //Generate one thousand tail particles, fill main vector.
{
	for (int i = 0; i < 10000; i++)
	{
		sf::CircleShape holder;
		holder.setPosition(0, 0);
		holder.setRadius(1.5);
		holder.setFillColor(sf::Color::Green);
		//		holder.time = time(0);
		tail.push_back(holder);
	}
}

struct vec_2 { //Declare 2D vec structure. 
	double x = 0;
	double y = 0;
};

struct body { //Declare standard planet holder, with one vector for position, one of velocity, and a double for mass.
	vec_2 pos;
	vec_2 vel;
	double mass = 1;
};


void update_pos(body &bd1, body bd2) //Take in bodies
{
	vec_2 pos = bd1.pos; //Extract relevant position, velocity data
	vec_2 vel = bd1.vel;
	vec_2 C = bd2.pos;

	pos.x -= C.x; //Calculate body position relative other body
	pos.y -= C.y;

	double k_true = pos.y / pos.x;		//Calculate "true" K-value
	double precision = 100000;			//Define physics detail
	bool flip_r = false;				//Prepare variables for X/Y switching
	bool flip_x = false;

	if (abs(k_true) > 1.5)			//Determine if K-value is too large for effective handling
	{
		flip_r = true;				//Note down that coordinates have been rotated
		double x_l = pos.x;			//Swap X/Y Coordinates
		pos.x = pos.y;
		pos.y = x_l;
	}

	if (pos.x < 0)					//Determine if coordinates are negative
	{
		flip_x = true;				//Note that coordinates have been flipped
		pos.x = -pos.x;				//Flip X coordinate
	}

	double k = pos.y / pos.x;		//Calculate new K value, post coordinate flipping/rotating
	double v = atan(k);				//Calculate new angle

	double x_rel = cos(v);			//Split the angle into two constituent coordinates
	double y_rel = sin(v);

	double r = sqrt(pow((pos.x), 2) + pow((pos.y), 2));		//Calculate the distance between the bodies
	double F = -bd2.mass / pow(r, 2);						//Determinte the magnitude of gravity, ie. the force

	if (flip_x) //If flipped, un-flip X
	{
		x_rel = -x_rel;
	}

	if (flip_r) //If rotated, un-rotate coordinates
	{
		double x_rel_l = x_rel;
		x_rel = y_rel;
		y_rel = x_rel_l;
	}


	bd1.vel.x += x_rel * F / precision; //Apply the force to the velocity
	bd1.vel.y += y_rel * F / precision;

	bd1.pos.x += vel.x / precision; //Apply the new velocity to the position
	bd1.pos.y += vel.y / precision;
}



void main()
{
	int fcount = 0;
	int scount = 0;
	long tcount = 0;

	body mb;
	body center;

	clock_t r_time = cur_time;
	clock_t fc_time = cur_time;

	mb.vel.x = 120;
	mb.pos.x = 500;
	mb.pos.y = 600;
	mb.mass = 250000;

	center.vel.x = 0;
	center.pos.x = 500;
	center.pos.y = 500;
	center.mass = 700000;

	fill_tail();

	Sleep(5000);
	while (true)
	{
		cur_time = clock();
		update_pos(mb, center);
		update_pos(center, mb);

		fcount++;

		if (cur_time - r_time > 0.02 * cps)
		{

			vec_2 pos_avg;
			pos_avg.x = (center.pos.x + mb.pos.x) / 2;
			pos_avg.y = (center.pos.y + mb.pos.y) / 2;

			vec_2 cpos;
			cpos.x = center.pos.x - pos_avg.x + 500;
			cpos.y = center.pos.y - pos_avg.y + 500;

			vec_2 mpos;
			mpos.x = mb.pos.x - pos_avg.x + 500;
			mpos.y = mb.pos.y - pos_avg.y + 500;

			window.clear();

			draw_circle(cpos.x, cpos.y, 0);//draw_circle(center.pos.x, center.pos.y, 0);
			draw_circle(mpos.x, mpos.y, 1);//draw_circle(mb.pos.x, mb.pos.y, 1);

			draw_tail(mpos.x, mpos.y, 1);//draw_tail(mb.pos.x, mb.pos.y, 1);
			draw_tail(cpos.x, cpos.y, 0);//draw_tail(center.pos.x, center.pos.y, 0);

			window.display();
			r_time = cur_time;
		}
		if (cur_time - fc_time > 10 * cps)
		{
			cout << fcount / 10 << endl;
			fcount = 0;
			fc_time = cur_time;
			scount++;
		}
		if (scount > 1)
		{
			scount = 0;
			Sleep(2);
		}
	}
}