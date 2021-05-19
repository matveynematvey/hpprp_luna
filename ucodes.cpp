#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include "ucenv.h"
#include "collide.inc"

#define RIGHT 2
#define DOWN 4
#define LEFT 8
#define UP 1
#define RP1 16
#define RP2 32
#define RP3 64
#define RP4 128
#define WIDTH 100
#define HEIGHT 100
#define ITERS 20
#define FIELD(i, j) field[(i)*(HEIGHT) + j]
#define NEW_FIELD(i, j) new_field[(i)*(HEIGHT) + j]

void fill(char *field, int i0, int j0, int i1, int j1,
		  double r, double d, double l, double u,
		  double r1, double r2, double r3, double r4, int n)
{
	int i, j;
	for (i = i0; i < i1; i++)
	{
		for (j = j0; j < j1; j++)
		{
			char cell = 0;
			if (drand48() < r)
				cell += RIGHT;
			if (drand48() < d)
				cell += DOWN;
			if (drand48() < l)
				cell += LEFT;
			if (drand48() < u)
				cell += UP;
			if (drand48() < r1)
				cell += RP1;
			if (drand48() < r2)
				cell += RP2;
			if (drand48() < r3)
				cell += RP3;
			if (drand48() < r4)
				cell += RP4;

			FIELD(i + 1, j)=cell;
		}
	}
}

extern "C"
{

	void save(int iter, InputDF &fi)
	{
		int i, j;
		char density_path[100], move_path[100], rest_path[100];
		FILE *fl_density, *fl_move, *fl_rest;

		sprintf(density_path, "density/%06d.xls", iter);
		sprintf(move_path, "move/%06d.xls", iter);
		sprintf(rest_path, "rest/%06d.xls", iter);
		fl_density=fopen(density_path, "a");
		fl_move=fopen(move_path, "a");
		fl_rest=fopen(rest_path, "a");

		char *field = fi.getData<char>();

		for (j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS-1; j++) {
		int rest=0, move=0;
		double new_density, new_move, new_rest;
		for (i=AVERAGING_RADIUS; i<HEIGHT-AVERAGING_RADIUS-1; i++) {
			sum_mass(i, j, &rest, &move);
		}

		
		new_density=1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
		new_move=1.0*move/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
		new_rest=1.0*rest/(HEIGHT-AVERAGING_RADIUS*2-1)/square;

		write(fl_density, j, ensemble, old_density, new_density);
		write(fl_move, j, ensemble, old_move, new_move);
		write(fl_rest, j, ensemble, old_rest, new_rest);

		
	}
		fclose(fl_density);
		fclose(fl_rest);
		fclose(fl_move);

		free(old_density);
		free(old_move);
		free(old_rest);
	}

	void c_iprint(int n, InputDF &fi)
	{
		printf("POLE %d\n", n);
		char *field = fi.getData<char>();
		for (int i = 0; i < 100; i++)
		{
			printf("%d ", field[i]);
		}
		printf("\n");
	}

	void init(int n, OutputDF &df)
	{
		char *field = df.create<char>(2700, 0); //полоса
		fill(field, 0, 0, 25, WIDTH, 0.7, 0.7, 0.7, 0.7, 0.25, 0, 0, 0, n);
	}

	void collide(int n, InputDF &dfi, OutputDF &dfo, OutputDF &dfo1, OutputDF &dfo2)
	{
		dfo = dfi;
		char *gr1 = dfo1.create<char>(WIDTH, 0);
		char *gr2 = dfo2.create<char>(WIDTH, 0);
		char *field = dfi.getData<char>();
		for (int i = 100; i < 2600; i++)
				field[i] = collideL(field[i]);
		for (int i = 100; i < 200; i++)
			gr1[i - 100] = field[i];
		for (int j = 2500, i = 0; j < 2600; j++, i++)
			gr2[i] = field[j];
		
	}

	void fillShadowEdges(int n, InputDF &dfi, InputDF &dfi1, InputDF &dfi2, OutputDF &dfo)
	{
		//char *fieldg = dfo.create<char>(2700, 0);
		dfo = dfi;
		char *field = dfi.getData<char>();
		char *gr1 = dfi1.getData<char>();
		char *gr2 = dfi2.getData<char>();
		
		for (int i = 0; i < 100; i++)
				field[i] = gr1[i];
		for (int j = 2600, i = 0; j < 2700; j++, i++)
				field[j] = gr2[i];
		
	}

	void nextLayer(int n, InputDF &dfi, OutputDF &dfo)
{
		char *field = dfi.getData<char>();
		char *new_field = dfo.create<char>(2700, 0);	
		dfo = dfi;
		int HEIGHT_N = HEIGHT / 4;
		for (int i=1; i < HEIGHT_N+1; i++) 
		{
			int j;
			for (j=0; j<WIDTH; j++) 
			{
			char cell=FIELD(i, j)&0xf0; // keep rest mass
			// propogate particles
			if (FIELD(i, (j+WIDTH-1)%WIDTH)&RIGHT) cell+=RIGHT; //field[(i)*(WIDTH) + j]
			if (FIELD(i, (j+1)%WIDTH)&LEFT) cell+=LEFT;
			if (FIELD((i+HEIGHT_N-1)%HEIGHT_N, j)&UP) cell+=UP;
			if (FIELD((i+1)%HEIGHT_N, j)&DOWN) cell+=DOWN;

			NEW_FIELD(i, j)=cell;
			}
		}
} 
}
