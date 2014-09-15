#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/queue.h>
#include <assert.h>
#include <getopt.h>

#define OUTSIZE 512
#ifdef DEBUG
#define debug(M, ...) fprintf(stdout, "DEBUG %s:%d: " M, __FILE__, __LINE__, ##__VA_ARGS__)
#else
#define debug(M, ...)
#endif

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
colour **pgm;

typedef struct { float x; float y; float h; int iteration; int parent; colour c; } point;
typedef struct { int howmany; } landscape;

/* Define a linked list of points */
struct list_item { point p; TAILQ_ENTRY(list_item) entries; };
/* We need an iteration pointer */
struct list_item *np, *np_2;

/* List of points we've created and interpolated */
TAILQ_HEAD(mainlist_h, list_item);

/* Working array for temporary points */
TAILQ_HEAD(templist_h, list_item);

int temp_howmany;

int greyscale = 0;

double *dsq;

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
dump_points(struct mainlist_h *l)
{
	int i = 0;
    TAILQ_FOREACH(np, l, entries) {
		debug("[MAIN] %d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, np->p.x, np->p.y, np->p.h, np->p.iteration, np->p.parent);
        i++;
	}
}

void
dump_temp(struct templist_h *l)
{
	int i = 0;
    TAILQ_FOREACH(np, l, entries) {
		debug("[TEMP] %d = <%.3f,%.3f,%.3f> g=%d p=%d\n", i, np->p.x, np->p.y, np->p.h, np->p.iteration, np->p.parent);
        i++;
	}
}

typedef struct { double d; int i; double h; } di;
di *sorted;

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
interpolate(struct templist_h *t, int howmany, struct mainlist_h *m, int w_howmany, int generation)
{
	int i=0,j;
	double reduction = pow(0.5, generation);
	
debug("generation %d, %d world points, %d temp points\n", generation, w.howmany, howmany);

    /* Loop over the templist to process each new point */
    TAILQ_FOREACH(np, t, entries) {
		double total_d = 0.0;
		point ip = np->p;
		double i_height = 0.0;
debug("interpolating new point %d: ", i);

        int j = 0;

		/* First we work out the individual distances and store them */
        TAILQ_FOREACH(np_2, m, entries) {
			point wp = np_2->p;
			double d = pow(wp.x-ip.x,2) + pow(wp.y-ip.y,2);
            if (d > 0.0) {
                sorted[j] = (di){ d, j, wp.h };
    debug("[PDIST] <%d, %.4f, %.4f> <%.3f,%.3f> <=> <%.3f,%.3f> %s\n", j, d, wp.h, wp.x,wp.y, ip.x,ip.y, (wp.x==ip.x && wp.y==ip.y) ? "SAME" : "NOSAME");
                j++;
            }
		}

        assert(w_howmany > 0);
		qsort(sorted, w_howmany, sizeof(di), compare_dist);
debug("[MINDIST] <%d, %.4f>\n", sorted[0].i, sorted[0].d);
			
		/* We only get here after one generation of spawning which
		   means we have at least 12 points to consider
 		*/
		int lim = w_howmany > 10 ? 10 : w_howmany;

		for(j=0; j<lim; j++) {
			double d = sorted[j].d;
            debug("[DIST] %d = %.3f (h=%.3f)\n", j, d, sorted[j].h);
			if (d > 0.0) { /* && d < 1.0) { */
				double dr2 = 1.0 / pow(d, 2);
				total_d = total_d + dr2;
				dsq[j] = dr2;
			}
			else {
				dsq[j] = 0.0;
			}
		}

debug("%d points, total=%.4f, first=%.4f ratio=%.4f\n", lim, total_d, sorted[0].d, 1.0/(pow(sorted[0].d,2)));

		for(j=0; j<lim; j++) {
			double r = dsq[j];
/*			double n_height = w.points[j].h * r; */
			double n_height = sorted[j].h * r;
			i_height = i_height + n_height;
debug("point %d, h=%.4f, invsqlaw=%.4f, adding=%.4f, total=%.4f\n", j, w.points[j].h, r, n_height, i_height);
		}

        assert(total_d != 0);
		np->p.h = i_height / total_d;

debug("final height = %.4f / %.4f = %.4f\n", i_height, total_d, t[i].h);

        assert(!isnan(np->p.h));
 		np->p.h = np->p.h + reduction * 0.1 * m1_p1(); 
		np->p.c = colour_by_height(np->p.h);
	}
}

int main(int argc, char **argv) {
	int i, j, q, gen;
	landscape world;
	static int seed;
    int initial_points = 5;
    static int iterations = 6;
    int maxpoints = initial_points;
    struct mainlist_h mainlist = TAILQ_HEAD_INITIALIZER(mainlist);
    struct templist_h templist = TAILQ_HEAD_INITIALIZER(templist);
    static int size = OUTSIZE;
    static double span = 0.75;

    seed = time(NULL);

    int c;
    while (1) {
        static struct option long_options[] = {
            {"iterations", optional_argument, &iterations, 'i'},
            {"size", optional_argument, &size, 'w'},
            {"seed", optional_argument, &seed, 's'},
            {"span", optional_argument, &span, 'v'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "i:s:w:v:", long_options, &option_index);

        if (c == -1) { break; }

        switch(c) {
            case 'i':
                printf("ITERATIONS = %s\n", optarg);
                iterations = atoi(optarg);
                break;
            case 'w':
                printf("OUTSIZE = %s\n", optarg);
                size = atoi(optarg);
                break;
            case 's':
                printf("SEED = %s\n", optarg);
                seed = atoi(optarg);
                break;
            case 'v':
                printf("SPAN = %s\n", optarg);
                span = atof(optarg);
                break;

            default:
                abort();
        }
    }

    for(i=1; i<iterations; i++) {
        maxpoints = maxpoints + 3*maxpoints;
    }
    printf("MAXPOINTS FOR n=%d, i=%d is %d\n", initial_points, iterations, maxpoints);

    dsq = (double *)malloc((maxpoints+1)*sizeof(double));
    assert(dsq != NULL);

    sorted = (di *)malloc((maxpoints+1)*sizeof(di));
    assert(sorted != NULL);

    colour (*pgm)[size] = malloc(sizeof *pgm * size);
    assert(pgm != NULL);

    debug("SEED %d\n", seed);
	srand(seed);

	/* Generate the initial set of 3 points around the origin */
	for(i=0; i<initial_points; i++) {
        struct list_item *f = malloc(sizeof(struct list_item));
		point new_point = (point){
			0.25*m1_p1(), 0.25*m1_p1(), 0.1*m1_p1(), 0, -1, {-1,-1,-1}
		};
		new_point.c = colour_by_height(new_point.h);
        f->p = new_point;
        TAILQ_INSERT_TAIL(&mainlist, f, entries);

		world.howmany = i+1;

	}
	dump_points(&mainlist);

	for(gen=1; gen<iterations; gen++) {
		double reduction = pow(0.8, gen);

			/* Iteration step */
			debug("Iteration step %d\n", gen);

			/* 0. Reset the count of how many temporary points we have */
			temp_howmany = 0;

			/* 1. Generate parent points for 3x new points */
			for(i=0; i<3*world.howmany; i++) {
				int parent_index = randi(world.howmany);
                int q = -1;
				point parent;

                TAILQ_FOREACH(np, &mainlist, entries) {
                    q++;
                    if (q == parent_index) {
                        parent = np->p;
                        break;
                    }
                }
                assert(q == parent_index);

				point new_point = (point){ 
					fmax( fmin( parent.x + reduction*0.5*m1_p1(), 1 ), -1 ),
					fmax( fmin( parent.y + reduction*0.5*m1_p1(), 1 ), -1 ),
					parent.h,
					1,
					parent_index,
                    (colour){-1,-1,-1} 
				};
                debug("[NEWP] %.3f\n", parent.h);
                struct list_item *f = malloc(sizeof(struct list_item));
                f->p = new_point;

                TAILQ_INSERT_TAIL(&templist, f, entries);
				temp_howmany++;
			}

			/* Interpolate the heights on the new points based on the old points */
			interpolate(&templist, temp_howmany, &mainlist, world.howmany, gen);

			/* Push the new points onto the end of the old points */
            TAILQ_FOREACH(np, &templist, entries) {
                np->p.c = colour_by_height(np->p.h);
//                debug("new point %d has height %.2f and <%d,%d,%d>\n", i, temp[i].h, temp[i].c.r, temp[i].c.g, temp[i].c.b);
            }

            dump_temp(&templist);

            TAILQ_CONCAT(&mainlist, &templist, entries);
            world.howmany = world.howmany + temp_howmany;

            dump_points(&mainlist);
	}

	/* dump_points(world); */
debug("World has %d points\n", world.howmany);
	/* pgm_voronoi(world); */

	for(i=0; i<size; i++) {
		for(j=0; j<size; j++) {
			pgm[i][j] = c_blank;
		}
	}

	double min_height = 10.0, max_height = -10.0;

#if 0
	for(i=0; i<world.howmany; i++) {
		int px = (int)((world.points[i].x + 1.0) * (OUTSIZE/2.0));
		int py = (int)((world.points[i].y + 1.0) * (OUTSIZE/2.0));

		min_height = fmin(min_height, world.points[i].h);
		max_height = fmax(max_height, world.points[i].h);

		pgm[px][py] = colour_by_height(world.points[i].h);
	}
#endif

    double span2 = 2*span;

	for(i=0; i<size; i++) {
		for(j=0; j<size; j++) {
            /* Uncoloured pixel needs voronoising */
            double min_dist = 999.0;
            point *min_point = NULL;
            double tx = -span + (span2*i)/size;
            double ty = -span + (span2*j)/size;

            TAILQ_FOREACH(np, &mainlist, entries) {
                point p = np->p;
                double d = pow(p.x-tx,2) + pow(p.y-ty,2);
                debug("<%.2f,%.2f> <=> <%.2f,%.2f> = %.5f (%.5f)\n", tx,ty, p.x,p.y, d, min_dist);
                if (d < min_dist) { min_dist = d; min_point = &(np->p); }
            }
            assert(min_point != NULL);

            debug("mp <%.2f,%.2f> is %.5f from <%.2f, %.2f>, {%d,%d,%d}\n", min_point->x, min_point->y, min_dist, tx,ty, min_point->c.r, min_point->c.g, min_point->c.b);
            if (min_point->c.r == -1) {
                min_point->c = colour_by_height(min_point->h);
            }

            pgm[i][j] = min_point->c;
		}
	}

	fprintf(stderr, "P3 %d %d 255\n", size, size);
	for(i=0; i<size; i++) {
		for(j=0; j<size; j++) {
			colour t = pgm[i][j];
			fprintf(stderr, "%d %d %d ", t.r, t.g, t.b);
			if (j % 15 == 14) { fprintf(stderr, "\n"); }
		}
		fprintf(stderr, "\n");
	}

	debug("minH = %.5f, maxH = %.5f\n", min_height, max_height);
}
