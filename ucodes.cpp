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
#define AVERAGING_RADIUS 2
#define FIELD(i, j) field[(i)*(HEIGHT) + j]
#define NEW_FIELD(i, j) new_field[(i)*(HEIGHT) + j]

size_t splits, iters;

void fill(char *field, int i0, int j0, int i1, int j1,
		  double r, double d, double l, double u,
		  double r1, double r2, double r3, double r4)
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

	inline int GetMove(char cell)
	{
		int mass=0;
		if (cell&RIGHT) mass+=1;
		if (cell&DOWN) mass+=1;
		if (cell&LEFT) mass+=1;
		if (cell&UP) mass+=1;
		return mass;
	}

	inline int GetRest(char cell)
	{
		int mass=0;
		if (cell&RP1) mass+=2;
		if (cell&RP2) mass+=4;
		if (cell&RP3) mass+=8;
		if (cell&RP4) mass+=16;
		return mass;
	}

	void SumMass(int row, int col, int* rest, int* move, char *field)
	{
		int i, j;
		for (i=-AVERAGING_RADIUS; i<=AVERAGING_RADIUS; i++) {
			for (j=-AVERAGING_RADIUS; j<=AVERAGING_RADIUS; j++) {
				*rest+=GetRest(FIELD(i+row, j+col));
				*move+=GetMove(FIELD(i+row, j+col));
			}
		}
	}

	double* load(char const* path)
	{
		int j;
		FILE *fl;
		double *data=(double*)malloc(
			(WIDTH-2*AVERAGING_RADIUS-1)*sizeof(double)
		);

		fl=fopen(path, "r");

		for(j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS-1; j++) {
			int row;
			double val;
			if (fscanf(fl, "%d\t%lf\n", &row, &val)!=2 || row!=j) {
				fprintf(stderr, "ERROR: File '%s' has illegal format",
					path);
				abort();
			} else {
				data[j-AVERAGING_RADIUS]=val;
			}
		}

		fclose(fl);
		return data;
	}

	void write(FILE *fl, int j, int ensemble, double *old, double val)
	{
		assert(ensemble? !!old: !old);

		fprintf(fl, "%d\t%lf\n", j, 
			ensemble?
				(val+ensemble*old[j-AVERAGING_RADIUS])/(ensemble+1) 
				: val
		);
	}

	void SaveToFiles(int iter, int layer, int ensemble, InputDF &fi)
	{
		char *field = fi.getData<char>();
		int i, j, rad1 = 0, rad2 = 0;
		char density_path[100], move_path[100], rest_path[100];
		FILE *fl_density, *fl_move, *fl_rest;
		double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1);
		double *old_density=0, *old_move=0, *old_rest=0;

		sprintf(density_path, "density/%06d.xls", iter);
		sprintf(move_path, "move/%06d.xls", iter);
		sprintf(rest_path, "rest/%06d.xls", iter);
		
		if (ensemble>0) {
			old_density=load(density_path);
			old_move=load(move_path);
			old_rest=load(rest_path);
		}

		if (layer == 1) 
		{
			rad1 = AVERAGING_RADIUS;
			rad2 = 0;
			fl_density=fopen(density_path, "w");
			fl_move=fopen(move_path, "w");
			fl_rest=fopen(rest_path, "w");
		}
		else
		{
		fl_density=fopen(density_path, "a");
		fl_move=fopen(move_path, "a");
		fl_rest=fopen(rest_path, "a");
		}

		if (layer == splits)
		{
			rad2 = AVERAGING_RADIUS;
			rad1 = 0;
		}

	
		for (i=rad1; i < (HEIGHT / splits) - rad2; i++) {
			int rest=0, move=0;
			double new_density, new_move, new_rest;
			for (j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS; j++) {
				SumMass(i, j, &rest, &move, field);
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
	}


	void c_iprint(int n)
	{
		printf("%d\n", n);
		//for (int i = 0; i < 10; i++)
		//{
		//	printf("%d ", test[i]);
		//}
	}

	void sources(int iter)
	{
			fill(0, WIDTH/2-SOURCE_WIDTH/2, HEIGHT, WIDTH/2+SOURCE_WIDTH/2, 1, 1, 1, 1, 0.75, 0, 0, 0);
	}

	void SetField(OutputDF &df)
	{
		char *field = df.create<char>(2700, 0); //полоса
		fill(field, 0, 0, HEIGHT / splits, WIDTH, 0.7, 0.7, 0.7, 0.7, 0.25, 0, 0, 0);
	}

	void InitIndex(OutputDF &df, int val)
	{
		df.setValue<int>(iters * 3 * val + val);
	}

	void InitSplitsAndIters(OutputDF &df1, int val1, OutputDF &df2, int val2)
	{
		splits = val2;
		iters = val1;
		df1.setValue<int>(val1);
		df2.setValue<int>(val2);
	}

	void CollideCells(InputDF &dfi, OutputDF &dfo, OutputDF &dfo1, OutputDF &dfo2)
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

	void FillEdges(InputDF &dfi, InputDF &dfi1, InputDF &dfi2, OutputDF &dfo)
	{
		dfo = dfi;
		char *field = dfi.getData<char>();
		char *gr1 = dfi1.getData<char>();
		char *gr2 = dfi2.getData<char>();
		
		for (int i = 0; i < 100; i++)
				field[i] = gr1[i];
		for (int j = 2600, i = 0; j < 2700; j++, i++)
				field[j] = gr2[i];
		
	}

	void AssemblyNewLayer(InputDF &dfi, OutputDF &dfo)
{
		char *field = dfi.getData<char>();
		char *new_field = dfo.create<char>(2700, 0);	
		dfo = dfi;
		int HEIGHT_N = HEIGHT / splits;
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
