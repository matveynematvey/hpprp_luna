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
#define ITERS 4
#define AVERAGING_RADIUS 1
#define SOURCE_WIDTH 20
#define FIELD(i, j) field[(i)*(WIDTH) + j] //i-строка j-столбец
#define NEW_FIELD(i, j) new_field[(i)*(WIDTH) + j]

size_t splits, iters, max_ensemble;

void segment_intersection(int a1, int b1, int a2, int b2, int &a3, int &b3){
	a3 = std::max(a1, a2);
	b3 = std::min(b1, b2);
	if(a3>b3)
		a3=b3=-1;
}

void fill(char *field, int i0, int j0, int i1, int j1,
		  double r, double d, double l, double u,
		  double r1, double r2, double r3, double r4, int split)
{

	int i00 = (split - 1)*HEIGHT/splits, i11=i00 + HEIGHT/splits; //start and end of the strip in global area
	int j00 = 0, j11 = WIDTH; //start and end of the strip strip in global area
	int i_s, i_e, j_s, j_e;
	segment_intersection(i00, i11, i0, i1, i_s, i_e);
	segment_intersection(j00, j11, j0, j1, j_s, j_e);
	i_s -= i00; //start of the strip in local area
	i_e -= i00; //end of the strip in local area
	
	int i, j;
	for (i = i_s; i < i_e; i++)
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
			
			FIELD(i, j)=cell;
		}
	}
}

void fillS(char *field, int i0, int j0, int i1, int j1,
		  double r, double d, double l, double u,
		  double r1, double r2, double r3, double r4, int split, int iter, int ens)
{
	int i00 = (split - 1)*HEIGHT/splits, i11=i00 + HEIGHT/splits; //start and end of the strip in global area
	int j00 = 0, j11 = WIDTH; //start and end of the strip strip in global area
	int i_s, i_e, j_s, j_e;
	segment_intersection(i00, i11, i0, i1, i_s, i_e);
	segment_intersection(j00, j11, j0, j1, j_s, j_e);
	i_s -= i00; //start of the strip in local area
	i_e -= i00; //end of the strip in local area
	
	int i, j;
	for (i = i_s; i < i_e; i++)
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
			
			//printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", FIELD(i + 1 - (HEIGHT / splits) * (split - 1), j), split, iter, i, j, ens, 0);
			//fflush(stdout);
			
			FIELD(i, j)=cell;
			
			//printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", FIELD(i + 1 - (HEIGHT / splits) * (split - 1), j), split, iter, i, j, ens, 1);
			//fflush(stdout);
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
				//printf("%d\t\n", (i + row)*HEIGHT + j + col);
				//fflush(stdout);
				*rest+=GetRest(FIELD(i+row, j+col));
				*move+=GetMove(FIELD(i+row, j+col));
			}
		}
	}

/*
	void SaveToFiles(int iter, int split, InputDF &fi)
	{
		int i, j;
		char density_path[100], move_path[100], rest_path[100];
		FILE *fl_density, *fl_move, *fl_rest;
		double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1), rad = 0;

		sprintf(density_path, "density/%06d.xls", iter);
		sprintf(move_path, "move/%06d.xls", iter);
		sprintf(rest_path, "rest/%06d.xls", iter);


		if (split == 1) 
		{
			rad = AVERAGING_RADIUS;
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
		if (split == splits) rad = AVERAGING_RADIUS;


		for (i = 0; i < (HEIGHT / splits) - rad; i++)
		{
		
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

*/

	void Calculation(int ensemble, int split, int iter, InputDF &fi1, InputDF &fi2, OutputDF &dfo)
	{
		//printf("%d\t%d\t%d\n", ensemble, split, iter);
		//fflush(stdout);
		double *DenMoveRest = dfo.create<double>(WIDTH*6, 0);
		double *OldDenMoveRest = fi2.getData<double>();
		char *field = fi1.getData<char>();
		int i, j, rad1 = 0, rad2 = 0, index = WIDTH;
		double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1);
		char density_path[100], move_path[100], rest_path[100];
		FILE *fl_density, *fl_move, *fl_rest;

		char field_path[100], old_path[100];
		
		FILE *fl_field, *fl_old;
		/*
		sprintf(field_path, "Testing/Calc/%d_%d_%d.txt", ensemble, split, iter);
		fl_field=fopen(field_path, "w");
		for(int i=0; i<((HEIGHT / splits) + 2) * WIDTH; i++){
			fprintf(fl_field, "%d\t%d\n", i, field[i]);
		}
		fclose(fl_field);
*/
		//НЕ ЗАПИСЫВАТЬ, ЕСЛИ СЛОЙ НЕ ПОСЛЕДНИЙ!!!!
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
		if(ensemble == max_ensemble){
				sprintf(density_path, "density/%06d.xls", iter);
				sprintf(move_path, "move/%06d.xls", iter);
				sprintf(rest_path, "rest/%06d.xls", iter);
				fl_density=fopen(density_path, "w");
				fl_move=fopen(move_path, "w");
				fl_rest=fopen(rest_path, "w");
		}
		
		for (j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS; j++) {
			int rest=0, move=0;
			for (i=rad1; i < (HEIGHT / splits) - rad2; i++ ) {
				SumMass(i, j, &rest, &move, field);
			}

			if(split == splits){
				double new_density=1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
				double new_move=1.0*move/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
				double new_rest=1.0*rest/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
				DenMoveRest[j + index * 3] = (new_density + OldDenMoveRest[j] +(ensemble - 1)*OldDenMoveRest[j + index * 3])/(ensemble);
				DenMoveRest[j + index * 4] = (new_move + OldDenMoveRest[j + index] +(ensemble - 1)*OldDenMoveRest[j + index * 4])/(ensemble);
				DenMoveRest[j + index * 5] = (new_rest + OldDenMoveRest[j + index * 2] + (ensemble - 1)*OldDenMoveRest[j + index * 5])/(ensemble);
			}
			else{
				DenMoveRest[j] = 1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square + OldDenMoveRest[j];
				DenMoveRest[j + index] = 1.0*move/(HEIGHT-AVERAGING_RADIUS*2-1)/square + OldDenMoveRest[j + index];
				DenMoveRest[j + index * 2] = 1.0*rest/(HEIGHT-AVERAGING_RADIUS*2-1)/square + OldDenMoveRest[j + 2*index];
				DenMoveRest[j + index * 3] = OldDenMoveRest[j + index * 3];
				DenMoveRest[j + index * 4] = OldDenMoveRest[j + index * 4];
				DenMoveRest[j + index * 5] = OldDenMoveRest[j + index * 5];
			}
			
			//printf("%d\t%d\t%d\t%d\n", split, iter, i, DenMoveRest[j]);
			//fflush(stdout);
			if(ensemble == max_ensemble){
				fprintf(fl_density, "%d\t%lf\n", j , DenMoveRest[j + index * 3]);
				fprintf(fl_move, "%d\t%lf\n", split, DenMoveRest[j + index * 4]);
				fprintf(fl_rest, "%d\t%lf\n", split, DenMoveRest[j + index * 5]);
			}
		}
		if(ensemble == max_ensemble){
			fclose(fl_density);
			fclose(fl_rest);
			fclose(fl_move);
		}
		
		/*sprintf(old_path, "Testing/Calc/old_%d_%d_%d.txt", ensemble, split, iter);
		fl_old=fopen(old_path, "w");
		for(int i =0; i< WIDTH * 3; i++)
			fprintf(fl_old, "%d\t%lf\t%lf\t%lf\t%lf\n", 
			i, OldDenMoveRest[i], OldDenMoveRest[i + index * 3], DenMoveRest[i], DenMoveRest[i  + index * 3]);
		fclose(fl_old);*/
	}

	void c_iprint(int n, int k, InputDF &fi)
	{
		//int data = fi.getValue<int>();
		char* data = fi.getData<char>();	
		for (int i=0; i < n; i++)
		{
			printf("%d\t%d\t%d\n",i, k, data[i]);
			fflush(stdout);
		}
		printf("\n\n");
			fflush(stdout);
	}

	void SetField(int split, int ens, OutputDF &df)
	{
		char *field = df.create<char>(((HEIGHT / splits) + 2) * WIDTH, 0); //полоса
		fill(field, (HEIGHT / splits) * (split - 1), 0, (HEIGHT / splits) * split, WIDTH, 0.7, 0.7, 0.7, 0.7, 0.25, 0, 0, 0, split);

		/*char field_path[100];
		FILE *fl_field;

		sprintf(field_path, "Testing/SetField/%d_%d.txt", ens, split);
		fl_field=fopen(field_path, "w");
		for(int i=0; i<((HEIGHT / splits) + 2) * WIDTH; i++){
			fprintf(fl_field, "%d\t%d\n", i, field[i]);
		}
		fclose(fl_field);
		*/
		
		/*
		char field_sum_path[100];
		FILE *fl_field_sum;

		sprintf(field_sum_path, "Testing/SetField/%d_%d_sum.txt", ens, split);
		fl_field_sum=fopen(field_sum_path, "w");
		int rad1 = 0, rad2 = 0;
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
		printf("%d\t%d\t%d\t%d\t%d\n", split, 1 + rad1, (HEIGHT / splits) - rad2 - 1, AVERAGING_RADIUS, WIDTH-AVERAGING_RADIUS);
		fflush(stdout);
		for (int i=1 + rad1; i < (HEIGHT / splits) - rad2 - 1; i++) {
			int rest=0, move=0;
			for (int j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS; j++) {
				SumMass(i, j, &rest, &move, field);
			}
			double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1);
			double sum = 1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
			fprintf(fl_field_sum, "%d\t%lf\n", i, sum);
		}*/
	}

	void SourceField(int split, int iter, int ens, InputDF &fi, OutputDF &df)
	{
		char *field = fi.getData<char>();
		df = fi;
		if (iter == 1){
			fillS(field, (HEIGHT / splits) * (split - 1), WIDTH/2-SOURCE_WIDTH/2, (HEIGHT / splits) * split, WIDTH/2+SOURCE_WIDTH/2, 1, 1, 1, 1, 0.75, 0, 0, 0, split, iter, ens);
			//printf("%d\t%d\t%d\t%d\t%d\t%d\n", split, iter, (HEIGHT / splits) * (split - 1), WIDTH/2-SOURCE_WIDTH/2, (HEIGHT / splits) * split, WIDTH/2+SOURCE_WIDTH/2);
			//fflush(stdout);
		}

		/*
		char field_path[100];
		FILE *fl_field;

		sprintf(field_path, "Testing/SourceField/%d_%d_%d.txt", ens, split, iter);
		fl_field=fopen(field_path, "w");
		for(int i=0; i<((HEIGHT / splits) + 2) * WIDTH; i++){
			fprintf(fl_field, "%d\t%d\n", i, field[i]);
		}
		fclose(fl_field);*/

	}

	void InitFiction(OutputDF &df)
	{
		double *fict = df.create<double>(WIDTH*6, 0);
	}

	void InitSplitsAndIters(OutputDF &df1, int val1, OutputDF &df2, int val2, OutputDF &df3, int val3)
	{
		splits = val2;
		iters = val1;
		max_ensemble = val3;
		df1.setValue<int>(val1);
		df2.setValue<int>(val2);
		df3.setValue<int>(val3);
	}

	void CollideCells(int split, int iter, int ens, InputDF &dfi, OutputDF &dfo, OutputDF &dfo1, OutputDF &dfo2)
	{
		
		char *gr1 = dfo1.create<char>(WIDTH, 0);
		char *gr2 = dfo2.create<char>(WIDTH, 0);
		char *field = dfi.getData<char>();
		dfo=dfi;
		char *new_field = dfi.getData<char>();
		
		//for(int i=0; i<50; i++)
		//	printf("%d\t%d\n", i, new_field[i]);
		for (int i = WIDTH; i < ((HEIGHT / splits) + 1) * WIDTH; i++)
			new_field[i] = collideL(field[i]);
		for (int i = WIDTH; i < WIDTH*2; i++)
			gr1[i - WIDTH] = new_field[i];
		for (int j = ((HEIGHT / splits)) * WIDTH, i = 0; j < ((HEIGHT / splits) + 1) * WIDTH; j++, i++)
			gr2[i] = new_field[j];

	/*
	if(iter > 1){
		char field_path[100];
		FILE *fl_field;

		sprintf(field_path, "Testing/CollideField/%d_%d_%d.txt", ens, split, iter);
		fl_field=fopen(field_path, "w");
		for(int i=HEIGHT; i<((HEIGHT / splits) + 1) * WIDTH; i++){
			fprintf(fl_field, "%d\t%d\n", i, field[i]);
		}
		fclose(fl_field);
	}*/
		/*
		char field_sum_path[100];
		FILE *fl_field_sum;

		sprintf(field_sum_path, "Testing/CollideField/%d_%d_sum.txt", ens, split);
		fl_field_sum=fopen(field_sum_path, "w");
		int rad1 = 0, rad2 = 0;
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
		//printf("%d\t%d\t%d\t%d\t%d\n", split, 1 + rad1, (HEIGHT / splits) - rad2 - 1, AVERAGING_RADIUS, WIDTH-AVERAGING_RADIUS);
		fflush(stdout);
		for (int j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS; j++) {
			int rest=0, move=0;
			for (int i=1 + rad1; i < (HEIGHT / splits) - rad2 - 1; i++ ) {
				SumMass(i, j, &rest, &move, field);
			}
			double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1);
			double sum = 1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
			fprintf(fl_field_sum, "%d\t%lf\n", j, sum);
		}
		fclose(fl_field_sum);
		*/
	}

	void FillEdges(int split, int iter, int ens, InputDF &dfi, InputDF &dfi1, InputDF &dfi2, OutputDF &dfo)
	{
		
		char *field = dfi.getData<char>();
		char *gr1 = dfi1.getData<char>();
		char *gr2 = dfi2.getData<char>();
		dfo = dfi;
		char *new_field = dfi.getData<char>();
		//for(int i=0; i<50; i++)
		//	printf("%d\t%d\n", i, new_field[i]);
		for (int i = 0; i < WIDTH; i++)
				new_field[i] = gr1[i];
		for (int j = ((HEIGHT / splits) + 1) * WIDTH, i = 0; j < ((HEIGHT / splits) + 2) * WIDTH; j++, i++)
				new_field[j] = gr2[i];

		/*
		char field_path[100];
		FILE *fl_field;

		sprintf(field_path, "Testing/EdgesField/%d_%d_%d.txt", ens, split, iter);
		fl_field=fopen(field_path, "w");
		for(int i=0; i<((HEIGHT / splits) + 2) * WIDTH; i++){
			fprintf(fl_field, "%d\t%d\n", i, field[i]);
		}
		fclose(fl_field);*/
	}

	void AssemblyNewLayer(int split, int iter, int ens, InputDF &dfi, OutputDF &dfo)
{
		char *field = dfi.getData<char>();
		char *new_field = dfo.create<char>(((HEIGHT / splits) + 2) * WIDTH, 0);	
		//dfo = dfi;
		//char *new_field = dfi.getData<char>();
		
		//for(int i=0; i < ((HEIGHT / splits) + 2) * WIDTH; i++)
		//	new_field[i] = field[i];
		int HEIGHT_N = HEIGHT / splits;
		
		//if(iter == 2 && ens == 2 && split == 1)
		//for(int i=0; i< 200; i++)
		//	printf("%d\t%d\t%d\n", i, new_field[i], field[i]);
		char field_path1[100];
		FILE *fl_field1;
		sprintf(field_path1, "Testing/NewField/%d_%d_%d_n.txt", ens, split, iter);
		fl_field1=fopen(field_path1, "w");
		for (int i=1; i < HEIGHT_N+1; i++) 
		{
			int j;
			for (j=0; j<WIDTH; j++) 
			{
				char cell=FIELD(i, j)&0xf0; // keep rest mass
				fprintf(fl_field1, "%d\t%d\t", FIELD(i, j), cell);
				fflush(stdout);
				// propogate particles
				if (FIELD(i, (j+WIDTH-1)%WIDTH)&RIGHT) cell+=RIGHT; //field[(i)*(WIDTH) + j]
				if (FIELD(i, (j+1)%WIDTH)&LEFT) cell+=LEFT;
				if (FIELD((i+HEIGHT_N-1)%HEIGHT_N, j)&UP) cell+=UP;
				if (FIELD((i+1)%HEIGHT_N, j)&DOWN) cell+=DOWN;
				NEW_FIELD(i, j)=cell;
				fprintf(fl_field1, "%d\n", cell);
				fflush(stdout);
			}
		}
		fclose(fl_field1);
		/*if(iter == 1 && ens == 1 && split == 1)
		for(int i=0; i<100; i++){
			printf("%d\t%d\n", i, new_field[i]);
			
		}*/
		/*
		char field_path[100];
		FILE *fl_field;
		//if(iter == 1 && ens == 1 && split == 1)
		//for(int i=0; i< 100; i++)
			//printf("%d\t%d\t%d\n", i, new_field[i], 2);
		
		sprintf(field_path, "Testing/NewField/%d_%d_%d.txt", ens, split, iter);
		fl_field=fopen(field_path, "w");
		for(int i=0; i<((HEIGHT / splits) + 2) * WIDTH; i++){
				//printf("%d\t%d\n", i, new_field[i]);
			fprintf(fl_field, "%d\t%d\t%d\n", i, field[i], new_field[i]);
		}
		fclose(fl_field);*/
} 
}
