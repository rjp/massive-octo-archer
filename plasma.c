#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define MAXPOINTS 131072
#define OUTSIZE 512
#define OUTHEIGHT 255

double l_ocean = 0.15 * OUTHEIGHT;
double l_water = 0.20 * OUTHEIGHT;
double l_sand  = 0.25 * OUTHEIGHT;
double l_dirt  = 0.55 * OUTHEIGHT;
double l_rocks = 0.60 * OUTHEIGHT;
double l_snow  = 0.80 * OUTHEIGHT;

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

typedef struct { float x; float y; float h; int iteration; int parent; } point;
typedef struct { point points[MAXPOINTS]; int howmany; } landscape;

/* Working array for temporary points */
point temp[MAXPOINTS];
int temp_howmany;

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
		printf("%d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, w.points[i].x, w.points[i].y, w.points[i].h, w.points[i].iteration, w.points[i].parent);
	}
}

void
dump_temp()
{
	int i;
	for(i=0; i<temp_howmany; i++) {
		printf("%d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, temp[i].x, temp[i].y, temp[i].h, temp[i].iteration, temp[i].parent);
	}
}

void
interpolate(point t[], int howmany, landscape w, int generation)
{
	int i,j;
	double reduction = pow(0.5, generation);
	
	for(i=0; i<howmany; i++) {
		double dsq[MAXPOINTS]; /* enough headroom */
		double total_d = 0;
		point ip = t[i];
		double i_height = 0.0;
		double d_threshold = 9999999.0;
		
		if (w.howmany > 20) {
			d_threshold = 4.0;
		}

		/* First we work out the individual distances and total distance */
		for(j=0; j<w.howmany; j++) {	
			point wp = w.points[j];
			double d = pow(wp.x-ip.x,2) + pow(wp.y-ip.y,2);
			double dr2;
			if (d > 0.0) {
				if (d <= d_threshold) {
					dr2 = 1.0 / pow(d, 2);
					total_d = total_d + dr2;
					dsq[j] = dr2;
				}
				else {
					total_d = total_d + 0.0;
					dsq[j] = 0.0;
				}
			}
			else {
				total_d = total_d + 0.0;
				dsq[j] = 0.0;
			}
		}

		for(j=0; j<w.howmany; j++) {
			double d = dsq[j];
			double r = d;
			i_height = i_height + w.points[j].h * r;
		}

		assert(total_d != 0.0);

		t[i].h = i_height / total_d;
if(isnan(t[i].h)) {
	for(j=0; j<w.howmany; j++) {
		printf("%d dr2=%.6f\n", j, dsq[j]);
	}
	printf("%.6f / %.6f isnan\n", i_height, total_d);
}
assert(! isnan(t[i].h));
		t[i].h = t[i].h + reduction * 0.1 * m1_p1();
	}
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
colour_of_grass(void)
{
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
colour_of_dirt(void)
{
	return shade_of_x(c_dirt.r, c_dirt.g, c_dirt.b, 32);
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

colour
colour_by_height(double unscaled_height)
{
	double scaled_height = (int)(OUTHEIGHT * (unscaled_height + 0.2)/0.4);
	colour blocks = colour_of_grass();

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

int main(void) {
	int i, j, q, gen;
	landscape world;

	srand(time(NULL));

	/* Generate the initial set of 3 points around the origin */
	for(i=0; i<3; i++) {
		world.points[i] = (point){
			0.25*m1_p1(), 0.25*m1_p1(), 0.1*m1_p1(), 0, -1
		};
		world.howmany = i+1;
	}
	puts("Initial world");
	dump_points(world);

	for(gen=1; gen<6; gen++) {
		double reduction = pow(0.75, gen);

			/* Iteration step */
			printf("Iteration step %d\n", gen);

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
printf("World has %d points\n", world.howmany);
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
if (px == 0) {
	printf("Left EDGE: <%.5f,%.5f>\n", world.points[i].x, world.points[i].y);
}
		min_height = fmin(min_height, world.points[i].h);
		max_height = fmax(max_height, world.points[i].h);

		pgm[px][py] = colour_by_height(world.points[i].h);
	}

	for(i=0; i<OUTSIZE; i++) {
		for(j=0; j<OUTSIZE; j++) {
			colour t = pgm[i][j];
			if (t.a == 0) {
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
				srand(1000*min_point.x + 1000*min_point.y);
				pgm[i][j] = colour_by_height(min_point.h);
			}

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

	printf("minH = %.5f, maxH = %.5f\n", min_height, max_height);
}
