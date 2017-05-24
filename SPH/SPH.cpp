// SPH.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "sph_fluid_system.h"
#include <iostream>


using namespace std;

#define PARTICLE_COUNTS			40960*2


int _tmain(int argc, _TCHAR* argv[])
{

	SPH::float_3 wall_min = { -20, 4, -20 };
	SPH::float_3 wall_max = { 20, 50, 20 };
	SPH::float_3 fluid_min = { -15, 5, -15 };
	SPH::float_3 fluid_max = { 15, 35, 15 };
	SPH::float_3 ghost_min = { -20, 4, -20 };
	SPH::float_3 ghost_max = { 20, 6, 20 };
	SPH::float_3 gravity = { 0.0f, 0.0f, 0.0f };

	SPH::FluidSystem System;


	cout << &wall_min << endl;

	System.init(PARTICLE_COUNTS, &wall_min, &wall_max, &fluid_min, &fluid_max, &ghost_min, &ghost_max, &gravity);
	System.output(000);


	int end_step;
	cout << "How many steps would you like to run?" << endl;
	cin >> end_step;
	int freq;
	cout << "At what frequency would you like to output to file" << endl;
	cin >> freq;

	long double name = 9999;

	for (int i = 0; i < end_step; ++i) {
		System.tick();
		if ((i%freq) == 0) {
			System.output(name);
			--name;
		}
	}