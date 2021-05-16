char collide1(char cell)
{
	double r=drand48();
	switch(cell) {
		case 5:
			if (r<0.33335) return 16;
			if (r<0.6667) return 10;
			return 5;
		case 7:
			if (r<0.5) return 18;
			return 7;
		case 10:
			if (r<0.33335) return 16;
			if (r<0.66665) return 10;
			return 5;
		case 11:
			if (r<0.5) return 17;
			return 11;
		case 13:
			if (r<0.5) return 24;
			return 13;
		case 14:
			if (r<0.5) return 20;
			return 14;
		case 15:
			if (r<0.33335) return 26;
			if (r<0.6667) return 21;
			return 15;
		case 16:
			if (r<0.3333) return 16;
			if (r<0.66665) return 10;
			return 5;
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
			if (r<0.33335) return 26;
			if (r<0.66665) return 21;
			return 15;
		case 24:
			if (r<0.5) return 24;
			return 13;
		case 26:
			if (r<0.3333) return 26;
			if (r<0.66665) return 21;
			return 15;
		default:
			return cell;
	}
}
