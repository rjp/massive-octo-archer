#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

/* C logging from CCAN/Junkcode from:
 * https://github.com/rustyrussell/ccan/tree/master/junkcode/henryeshbaugh%40gmail.com-log
 */
#include "log/log.h"

#define print_log(M, ...)

#define MAXPOINTS 8192
#define OUTSIZE 512
#define OUTHEIGHT 255

double l_ocean = 0.10;
double l_water = 0.20;
double l_sand  = 0.25;
double l_dirt  = 0.65;
double l_rocks = 0.70;
double l_snow  = 0.80;

typedef struct { int r; int g; int b; int a; } colour;
colour c_blank = { 0, 0, 0, 0 };
colour c_ocean = { 0, 0, 128, 1 };
colour c_water = { 0, 0, 255, 1 };
colour c_sand  = { 255, 192, 64, 1 };
colour c_dirt  = { 128, 64, 0, 1 };
colour c_rocks = { 128, 128, 128, 1 };
colour c_grass = { 0, 255, 0, 1 };
colour c_snow  = { 255, 255, 255, 1 };
colour c_bork  = { 255, 0, 0, 1};

/* Strictly a ppm now but meh */
colour pgm[OUTSIZE][OUTSIZE];

typedef struct { float x; float y; float h; int iteration; int parent; colour c; } point;
typedef struct { point points[MAXPOINTS]; int howmany; } landscape;

/* Working array for temporary points */
point temp[MAXPOINTS];
int temp_howmany;

int greyscale = 0;

/* Lifted from
 * http://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
 */
double
randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}

double
m1_p1(void)
{
	double n = randn(0.0, 1.0);
	return n;	
}

/* Return a uniform random range between [0,maxi) */
int
randi(int maxi)
{
	double n = rand();
	double n_scaled = n / (RAND_MAX+1.0);
	int to_scale = n_scaled * maxi;
	
	assert(to_scale > -1);
	assert(to_scale < maxi);

	return to_scale;
}

void
dump_points(landscape w)
{
	int i;
	for(i=0; i<w.howmany; i++) {
		print_log(LOG_INFO, "%d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, w.points[i].x, w.points[i].y, w.points[i].h, w.points[i].iteration, w.points[i].parent);
	}
}

void
dump_temp()
{
	int i;
	for(i=0; i<temp_howmany; i++) {
		print_log(LOG_INFO, "%d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, temp[i].x, temp[i].y, temp[i].h, temp[i].iteration, temp[i].parent);
	}
}

typedef struct { double d; int i; } di;
di sorted[MAXPOINTS];

int
compare_dist(const void *a, const void *b)
{
	const di *va = a;
	const di *vb = b;

	if (vb->d > va->d) { return -1; }
	if (vb->d < va->d) { return +1; }
	return 0;
}

colour
shade_of_x(int r, int g, int b, int ratio)
{
	int m = randi(2*ratio) - ratio;
	r = fmax(0, fmin(255, r + m));
	g = fmax(0, fmin(255, g + m));
	b = fmax(0, fmin(255, b + m));
	return (colour){ r, g, b, 1 };
}

colour
colour_of_x(int r, int g, int b, int rr, int gr, int br)
{
	r = fmax(0, fmin(255, r + randi(2*rr) - rr));
	g = fmax(0, fmin(255, g + randi(2*gr) - gr));
	b = fmax(0, fmin(255, b + randi(2*br) - br));
	return (colour){ r, g, b, 1 };
}

colour
colour_of_dirt(void)
{
	return shade_of_x(c_dirt.r, c_dirt.g, c_dirt.b, 32);
}

colour
colour_of_grass(void)
{
	/* 1d20 > 18 for grass -> dirt transformation */
	if (randi(20) > 19) {
		return colour_of_dirt();
	}
		
	return colour_of_x(32, 192, 32, 32, 32, 32);
}

colour
colour_of_snow(void)
{
	/* We want all the sliders to move up and down in lockstep */
	return shade_of_x(255, 255, 255, 32);
}

colour
colour_of_rocks(void)
{
	return shade_of_x(c_rocks.r, c_rocks.g, c_rocks.b, 32);
}

colour
colour_of_sand(void)
{
	return shade_of_x(c_sand.r, c_sand.g, c_sand.b, 16);
}

colour
colour_of_water(void)
{
	return colour_of_x(c_water.r, c_water.g, c_water.b, 32, 32, 32);
}

/* We need a clever multi-stage MACRO here to auto-generate these for us */
int
clamp(int in, int low, int high)
{
	if (in < low) { return low; }
	if (in > high) { return high; }
	return in;
}

colour
colour_by_height(double unscaled_height)
{
	double scaled_height = fmax( fmin( (unscaled_height + 0.2)/0.4, 1.0 ), 0.0 );
	colour blocks = colour_of_grass();

	if (greyscale) {
		int x = clamp((int)(OUTHEIGHT * scaled_height), 0, OUTHEIGHT);
		return (colour){x, x, x, 1};
	}

	if (scaled_height < l_sand) {
		blocks = colour_of_sand();
	}

	if (scaled_height < l_water) {
		blocks = colour_of_water();
	}

	if (scaled_height < l_ocean) {
		blocks = c_ocean;
	}

	if (scaled_height > l_dirt) { 
		blocks = colour_of_dirt();
	}

	if (scaled_height > l_rocks) { 
		blocks = colour_of_rocks();
	}

	if (scaled_height > l_snow) { 
		 blocks = colour_of_snow();
	}

	return blocks;
}


void
interpolate(point t[], int howmany, landscape w, int generation)
{
	int i,j;
	double reduction = pow(0.5, generation);
	
print_log(LOG_INFO, "generation %d, %d world points, %d temp points\n", generation, w.howmany, howmany);

	for(i=0; i<howmany; i++) {
		double dsq[MAXPOINTS]; /* enough headroom */
		double total_d = 0.0;
		point ip = t[i];
		double i_height = 0.0;
print_log(LOG_INFO, "interpolating new point %d: ", i);

		/* First we work out the individual distances and store them */
		for(j=0; j<w.howmany; j++) {	
			point wp = w.points[j];
			double d = pow(wp.x-ip.x,2) + pow(wp.y-ip.y,2);
			sorted[j] = (di){ d, j };
print_log(LOG_INFO, "<%d, %.4f, %.4f> ", j, d, wp.h);
		}

		qsort(sorted, w.howmany, sizeof(di), compare_dist);
print_log(LOG_INFO, "<%d, %.4f>\n", sorted[0].i, sorted[0].d);
			
		/* We only get here after one generation of spawning which
		   means we have at least 12 points to consider
 		*/
		int lim = w.howmany > 10 ? 10 : w.howmany;

		for(j=0; j<lim; j++) {
			double d = sorted[j].d;
			if (d > 0.0) { /* && d < 1.0) { */
				double dr2 = 1.0 / pow(d, 2);
				total_d = total_d + dr2;
				dsq[j] = dr2;
			}
			else {
				dsq[j] = 0.0;
			}
		}

print_log(LOG_INFO, "%d points, total=%.4f, first=%.4f ratio=%.4f\n", lim, total_d, sorted[0].d, 1.0/(pow(sorted[0].d,2)));

		for(j=0; j<lim; j++) {
			double r = dsq[j];
/*			double n_height = w.points[j].h * r; */
			double n_height = w.points[sorted[j].i].h * r;
			i_height = i_height + n_height;
print_log(LOG_INFO, "point %d, h=%.4f, invsqlaw=%.4f, adding=%.4f, total=%.4f\n", j, w.points[j].h, r, n_height, i_height);
		}

		t[i].h = i_height / total_d;

print_log(LOG_INFO, "final height = %.4f / %.4f = %.4f\n", i_height, total_d, t[i].h);

 		t[i].h = t[i].h + reduction * 0.1 * m1_p1(); 
		t[i].c = colour_by_height(t[i].h);
	}
}

int main(int argc, char **argv) {
	int i, j, q, gen;
	landscape world;
	int seed = time(NULL);

	/* In theory, this prevents the debugging from being printed */
	// set_log_mode(LOG_WARNING);

	if (argc > 1) {
		seed = atoi(argv[1]);
		print_log(LOG_INFO, "SEED %d\n", seed);
	}

	srand(seed);

	/* Generate the initial set of 3 points around the origin */
	for(i=0; i<5; i++) {
		world.points[i] = (point){
			0.25*m1_p1(), 0.25*m1_p1(), 0.1*m1_p1(), 0, -1
		};
		world.points[i].c = colour_by_height(world.points[i].h);
		world.howmany = i+1;
	}
	puts("Initial world");
	dump_points(world);

	for(gen=1; gen<6; gen++) {
		double reduction = pow(0.8, gen);

			/* Iteration step */
			print_log(LOG_INFO, "Iteration step %d\n", gen);

			/* 0. Reset the count of how many temporary points we have */
			temp_howmany = 0;

			/* 1. Generate parent points for 3x new points */
			for(i=0; i<3*world.howmany; i++) {
				int parent_index = randi(world.howmany);
				point parent = world.points[parent_index];
				temp[i] = (point){ 
					fmax( fmin( parent.x + reduction*0.5*m1_p1(), 1 ), -1 ),
					fmax( fmin( parent.y + reduction*0.5*m1_p1(), 1 ), -1 ),
					parent.h,
					1,
					parent_index
				};
				temp_howmany++;
			}

			/* Interpolate the heights on the new points based on the old points */
			interpolate(temp, temp_howmany, world, gen);

			/* Push the new points onto the end of the old points */
			for(i=0; i<temp_howmany; i++) {
				world.points[world.howmany] = temp[i];
				world.howmany++;
			}
	}

	puts("Final world");
	dump_points(world);
print_log(LOG_INFO, "World has %d points\n", world.howmany);
	/* pgm_voronoi(world); */

	for(i=0; i<OUTSIZE; i++) {
		for(j=0; j<OUTSIZE; j++) {
			pgm[i][j] = c_blank;
		}
	}

	double min_height = 10.0, max_height = -10.0;

	for(i=0; i<world.howmany; i++) {
		int px = (int)((world.points[i].x + 1.0) * (OUTSIZE/2.0));
		int py = (int)((world.points[i].y + 1.0) * (OUTSIZE/2.0));

		min_height = fmin(min_height, world.points[i].h);
		max_height = fmax(max_height, world.points[i].h);

		pgm[px][py] = colour_by_height(world.points[i].h);
	}

	for(i=0; i<OUTSIZE; i++) {
		for(j=0; j<OUTSIZE; j++) {
			colour t = pgm[i][j];
/*			if (t.a == 0) { */
				/* Uncoloured pixel needs voronoising */
				double min_dist = 999.0;
				point min_point;
				double tx = i/(OUTSIZE/2.0) - 1.0;
				double ty = j/(OUTSIZE/2.0) - 1.0;

				for(q=0; q<world.howmany; q++) {
					point p = world.points[q];
					double d = pow(p.x-tx,2) + pow(p.y-ty,2);
					if (d < min_dist) { min_dist = d; min_point = p; }
				}
				pgm[i][j] = min_point.c;
		 /* } */	
			/* Highlight the original points */
/* 			if (t.a == 1) { pgm[i][j] = c_blank; } */
		}
	}

	fprintf(stderr, "P3 %d %d 255\n", OUTSIZE, OUTSIZE);
	for(i=0; i<OUTSIZE; i++) {
		for(j=0; j<OUTSIZE; j++) {
			colour t = pgm[i][j];
			fprintf(stderr, "%d %d %d ", t.r, t.g, t.b);
			if (j % 15 == 14) { fprintf(stderr, "\n"); }
		}
		fprintf(stderr, "\n");
	}

	print_log(LOG_INFO, "minH = %.5f, maxH = %.5f\n", min_height, max_height);
}
