#include "ucenv.h"

extern "C" {

char collideL(char cell)
{
	double r=drand48();
	switch(cell) {
		case 5:
			if (r<0.3333) return 5;
			if (r<0.66665) return 10;
			return 16;
		case 7:
			if (r<0.5) return 7;
			return 18;
		case 10:
			if (r<0.3333) return 10;
			if (r<0.66665) return 5;
			return 16;
		case 11:
			if (r<0.5) return 11;
			return 17;
		case 13:
			if (r<0.5) return 13;
			return 24;
		case 14:
			if (r<0.5) return 14;
			return 20;
		case 15:
			if (r<0.3333) return 15;
			if (r<0.66665) return 21;
			return 26;
		case 16:
			if (r<0.3333) return 16;
			if (r<0.66665) return 5;
			return 10;
		case 17:
			if (r<0.5) return 17;
			return 11;
		case 18:
			if (r<0.5) return 18;
			return 7;
		case 20:
			if (r<0.5) return 20;
			return 14;
		case 21:
			if (r<0.3333) return 21;
			if (r<0.66665) return 15;
			return 26;
		case 24:
			if (r<0.5) return 24;
			return 13;
		case 26:
			if (r<0.3333) return 26;
			if (r<0.66665) return 15;
			return 21;
		default:
			return cell;
	}
}

}
