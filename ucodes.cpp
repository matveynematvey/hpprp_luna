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
#define WIDTH 500
#define HEIGHT 500
#define AVERAGING_RADIUS 1
#define SOURCE_WIDTH 10
#define FIELD(i, j) field[(i) * (WIDTH) + j] //i-строка j-столбец
#define NEW_FIELD(i, j) new_field[(i) * (WIDTH) + j]

size_t splits, iters, max_ensemble;

void SegmentIntersection(int a1, int b1, int a2, int b2, int &a3, int &b3)
{
	a3 = std::max(a1, a2);
	b3 = std::min(b1, b2);
	if (a3 > b3)
		a3 = b3 = -1;
}

void Fill(char *field, int i0, int j0, int i1, int j1,
		  double r, double d, double l, double u,
		  double r1, double r2, double r3, double r4, int split)
{

	int i00 = (split - 1) * HEIGHT / splits, i11 = i00 + HEIGHT / splits; //start and end of the strip in global area
	int j00 = 0, j11 = WIDTH;											  //start and end of the strip strip in global area
	int i_s, i_e, j_s, j_e;
	SegmentIntersection(i00, i11, i0, i1, i_s, i_e);
	SegmentIntersection(j00, j11, j0, j1, j_s, j_e);
	i_s -= i00; //start of the strip in local area
	i_e -= i00; //end of the strip in local area

	int i, j;
	for (i = i_s + AVERAGING_RADIUS; i < i_e + AVERAGING_RADIUS; i++)
	{
		for (j = j_s; j < j_e; j++)
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

			/*if(split == 3){
				if(i ){
					printf("%d\t%d\t%d\t%d\n", split, i + 1 - (HEIGHT / splits) * (split - 1), j, cell);
					fflush(stdout);
				}
			}*/

			FIELD(i + AVERAGING_RADIUS, j) = cell;
		}
	}
}

extern "C"
{

	inline int GetMove(char cell)
	{
		int mass = 0;
		if (cell & RIGHT)
			mass += 1;
		if (cell & DOWN)
			mass += 1;
		if (cell & LEFT)
			mass += 1;
		if (cell & UP)
			mass += 1;
		return mass;
	}

	inline int GetRest(char cell)
	{
		int mass = 0;
		if (cell & RP1)
			mass += 2;
		if (cell & RP2)
			mass += 4;
		if (cell & RP3)
			mass += 8;
		if (cell & RP4)
			mass += 16;
		return mass;
	}

	void SumMass(int row, int col, int *rest, int *move, char *field)
	{
		int i, j;
		for (i = -AVERAGING_RADIUS; i <= AVERAGING_RADIUS; i++)
		{
			for (j = -AVERAGING_RADIUS; j <= AVERAGING_RADIUS; j++)
			{
				*rest += GetRest(FIELD(i + row, j + col));
				*move += GetMove(FIELD(i + row, j + col));
			}
		}
	}


	void Calculation(int ensemble, int split, int iter, InputDF &fi1, InputDF &fi2, OutputDF &dfo)
	{
		double *DenMoveRest = dfo.create<double>(WIDTH * 6, 0);
		double *OldDenMoveRest = fi2.getData<double>();
		char *field = fi1.getData<char>();
		int i, j, rad1 = 0, rad2 = 0, index = WIDTH;
		double square = (AVERAGING_RADIUS * 2 + 1) * (AVERAGING_RADIUS * 2 + 1);
		char density_path[100], move_path[100], rest_path[100];
		FILE *fl_density, *fl_move, *fl_rest;

		char field_path[100], old_path[100];

		if (split == 1)
		{
			rad1 = AVERAGING_RADIUS;
			rad2 = 0;
		}
		if (split == splits)
		{
			rad2 = AVERAGING_RADIUS;
			rad1 = 0;
		}
		if (ensemble == max_ensemble)
		{
			sprintf(density_path, "density/%06d.xls", iter);
			sprintf(move_path, "move/%06d.xls", iter);
			sprintf(rest_path, "rest/%06d.xls", iter);
			fl_density = fopen(density_path, "w");
			fl_move = fopen(move_path, "w");
			fl_rest = fopen(rest_path, "w");
		}

		for (j = AVERAGING_RADIUS; j < WIDTH - AVERAGING_RADIUS - 1; j++)
		{
			int rest = 0, move = 0;
			for (i = AVERAGING_RADIUS + rad1; i < (HEIGHT / splits) - rad2 - AVERAGING_RADIUS - 1; i++)
			{
				SumMass(i, j, &rest, &move, field);
			}

			if (split == splits)
			{
				double new_density = 1.0 * (rest + move) / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square;
				double new_move = 1.0 * move / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square;
				double new_rest = 1.0 * rest / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square;
				DenMoveRest[j - AVERAGING_RADIUS + index * 3] = (new_density + OldDenMoveRest[j - AVERAGING_RADIUS] + (ensemble - 1) * OldDenMoveRest[j - AVERAGING_RADIUS + index * 3]) / (ensemble);
				DenMoveRest[j - AVERAGING_RADIUS + index * 4] = (new_move + OldDenMoveRest[j - AVERAGING_RADIUS + index] + (ensemble - 1) * OldDenMoveRest[j - AVERAGING_RADIUS + index * 4]) / (ensemble);
				DenMoveRest[j - AVERAGING_RADIUS + index * 5] = (new_rest + OldDenMoveRest[j - AVERAGING_RADIUS + index * 2] + (ensemble - 1) * OldDenMoveRest[j - AVERAGING_RADIUS + index * 5]) / (ensemble);
			}
			else
			{
				DenMoveRest[j - AVERAGING_RADIUS] = 1.0 * (rest + move) / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square + OldDenMoveRest[j - AVERAGING_RADIUS];
				DenMoveRest[j - AVERAGING_RADIUS + index] = 1.0 * move / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square + OldDenMoveRest[j - AVERAGING_RADIUS + index];
				DenMoveRest[j - AVERAGING_RADIUS + index * 2] = 1.0 * rest / (HEIGHT - AVERAGING_RADIUS * 2 - 1) / square + OldDenMoveRest[j - AVERAGING_RADIUS + 2 * index];
				DenMoveRest[j - AVERAGING_RADIUS + index * 3] = OldDenMoveRest[j - AVERAGING_RADIUS + index * 3];
				DenMoveRest[j - AVERAGING_RADIUS + index * 4] = OldDenMoveRest[j - AVERAGING_RADIUS + index * 4];
				DenMoveRest[j - AVERAGING_RADIUS + index * 5] = OldDenMoveRest[j - AVERAGING_RADIUS + index * 5];
			}

			if (ensemble == max_ensemble)
			{
				fprintf(fl_density, "%d\t%lf\n", j, DenMoveRest[j - AVERAGING_RADIUS + index * 3]);
				fprintf(fl_move, "%d\t%lf\n", split, DenMoveRest[j - AVERAGING_RADIUS + index * 4]);
				fprintf(fl_rest, "%d\t%lf\n", split, DenMoveRest[j - AVERAGING_RADIUS + index * 5]);
			}
		}
		if (ensemble == max_ensemble)
		{
			fclose(fl_density);
			fclose(fl_rest);
			fclose(fl_move);
		}
		if (!(iter%10) && split == splits && ensemble == max_ensemble)
			printf("Number of iteration: %d\r",iter);
	}

	void SetField(int split, OutputDF &df)
	{
		char *field = df.create<char>(((HEIGHT / splits) + 2 * AVERAGING_RADIUS) * WIDTH, 0); //полоса
		Fill(field, 0, 0, HEIGHT, WIDTH, 0.7, 0.7, 0.7, 0.7, 0.25, 0, 0, 0, split);
	}

	void SourceField(int split, int iter, InputDF &fi, OutputDF &df)
	{
		char *field = fi.getData<char>();
		df = fi;
		if (iter == 10)
			Fill(field, 0, WIDTH / 2 - SOURCE_WIDTH / 2, HEIGHT, WIDTH / 2 + SOURCE_WIDTH / 2, 1, 1, 1, 1, 0.75, 0, 0, 0, split);
	}

	void InitFiction(OutputDF &df)
	{
		double *fict = df.create<double>(WIDTH * 6, 0);
	}

	void InitSplitsAndIters(OutputDF &df1, int val1, OutputDF &df2, int val2, OutputDF &df3, int val3)
	{
		splits = val2;
		while (HEIGHT%splits)
		{
			printf("Incorrect value of splits\nEnter the correct value of splits: ");
			scanf("%lu", &splits);
			val2 = splits;
		}	
		iters = val1;
		max_ensemble = val3;
		df1.setValue<int>(val1);
		df2.setValue<int>(val2);
		df3.setValue<int>(val3);
		printf("MaxEnseble : %lu, Iterations :  %lu, Splits :  %lu, WIDTH :  %lu, HEIGHT :  %lu\n", max_ensemble, iters, splits, WIDTH, HEIGHT);
	}

	void CollideCells(int split, InputDF &dfi, OutputDF &dfo, OutputDF &dfo1, OutputDF &dfo2)
	{
		char *gr1 = dfo1.create<char>(WIDTH * AVERAGING_RADIUS, 0);
		char *gr2 = dfo2.create<char>(WIDTH * AVERAGING_RADIUS, 0);
		char *field = dfi.getData<char>();
		dfo = dfi;
		char *new_field = dfi.getData<char>();

		for (int i = WIDTH * AVERAGING_RADIUS; i < ((HEIGHT / splits) + AVERAGING_RADIUS) * WIDTH; i++)
			new_field[i] = collideL(field[i]);
		for (int i = WIDTH * AVERAGING_RADIUS; i < WIDTH * (2 * AVERAGING_RADIUS); i++)
			gr1[i - WIDTH * AVERAGING_RADIUS] = new_field[i];
		for (int j = ((HEIGHT / splits)) * WIDTH, i = 0; j < ((HEIGHT / splits) + AVERAGING_RADIUS) * WIDTH; j++, i++)
			gr2[i] = new_field[j];
	}

	void FillEdges(int split, InputDF &dfi, InputDF &dfi1, InputDF &dfi2, OutputDF &dfo)
	{
		char *field = dfi.getData<char>();
		char *gr1 = dfi1.getData<char>();
		char *gr2 = dfi2.getData<char>();
		dfo = dfi;
		char *new_field = dfi.getData<char>();
		//for(int i=0; i<50; i++)
		//	printf("%d\t%d\n", i, new_field[i]);
		for (int i = 0; i < WIDTH * AVERAGING_RADIUS; i++)
			new_field[i] = gr1[i];
		for (int j = ((HEIGHT / splits) + AVERAGING_RADIUS) * WIDTH, i = 0; j < ((HEIGHT / splits) + 2 * AVERAGING_RADIUS) * WIDTH; j++, i++)
			new_field[j] = gr2[i];
	}

	void AssemblyNewLayer(int split, InputDF &dfi, OutputDF &dfo)
	{
		char *field = dfi.getData<char>();
		char *new_field = dfo.create<char>(((HEIGHT / splits) + 2 * AVERAGING_RADIUS) * WIDTH, 0);

		int HEIGHT_N = HEIGHT / splits;

		char field_path1[100];
		FILE *fl_field1;

		for (int i = AVERAGING_RADIUS; i < HEIGHT_N + AVERAGING_RADIUS; i++)
		{
			int j;
			for (j = 0; j < WIDTH; j++)
			{
				char cell = FIELD(i, j) & 0xf0; // keep rest mass
				if (FIELD(i, (j + WIDTH - 1) % WIDTH) & RIGHT)
					cell += RIGHT; //field[(i)*(WIDTH) + j]
				if (FIELD(i, (j + 1) % WIDTH) & LEFT)
					cell += LEFT;
				if (FIELD((i + HEIGHT_N - 1) % HEIGHT_N, j) & UP)
					cell += UP;
				if (FIELD((i + 1) % HEIGHT_N, j) & DOWN)
					cell += DOWN;
				NEW_FIELD(i, j) = cell;
			}
		}
	}
}
