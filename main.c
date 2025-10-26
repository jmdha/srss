#include <math.h>
#include <stdbool.h>
#include <float.h>

#define RGFW_IMPLEMENTATION
#include "RGFW.h"

#define WIDTH   1024
#define HEIGHT  1024
#define GWIDTH  128
#define GHEIGHT 128
#define CWIDTH  (WIDTH / GWIDTH)
#define CHEIGHT (HEIGHT / GHEIGHT)
#define max(a, b)      ((a) > (b) ? (a) : (b))
#define min(a, b)      ((a) < (b) ? (a) : (b))
#define clamp(x, v, V) (min(max((x), (v)), (V)))

const u8 WHITE[4] = { 255, 255, 255, 255 };
const u8 BLACK[4] = {   0,   0,   0, 255 };

u8*           BUF;
RGFW_window*  WIN;
RGFW_surface* SURFACE;
bool          SOLID[GWIDTH * GHEIGHT] = { 0 };
float         SMOKE[GWIDTH * GHEIGHT] = { 0 };
float         VX[GWIDTH * GHEIGHT]    = { 0 };
float         VY[GWIDTH * GHEIGHT]    = { 0 };

double ts(void) {
	struct timespec t;
	timespec_get(&t, TIME_UTC);
	return (double) t.tv_sec * 1e9 + t.tv_nsec;
}

double dt_get(void) {
	static double t0;
	const  double t1 = ts();
	const  double v  = (t1 - t0) / 1e9;
	              t0 = t1;
	return v;
}

//void render_line(i32 x0, i32 y0, i32 x1, i32 y1, const u8 color[4]) {
//	const float  h = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
//	const float  mx = 1 / h * (x1 - x0);
//	const float  my = 1 / h * (y1 - y0);
//	for (size_t i = 0; i < h; i++) {
//		i32 x = x0 + i * mx;
//		i32 y = y0 + i * my;
//		if (x < 0 || y < 0 || x >= WIDTH || y >= HEIGHT) continue;
//		i32 c = y * (4 * WIDTH) + 4 * x;
//		memcpy(&BUF[c], color, 4 * sizeof(u8));
//	}
//}

void render_rect(i32 x, i32 y, i32 w, i32 h, const u8 color[4]) {
	for(i32 _x = x; _x < w + x; _x++) {
		for(i32 _y = y; _y < h + y; _y++) {
			i32 c = _y * (4 * WIDTH) + _x * 4;
			memcpy(&BUF[c], color, 4 * sizeof(u8));
		}
	}
}

//void render_grid(const u8 color[4]) {
//	for (size_t x = 1; x < GWIDTH; x++)
//		render_line(x * CWIDTH, 0, x * CHEIGHT, HEIGHT - 1, WHITE);
//	for (size_t y = 1; y < GHEIGHT; y++)
//		render_line(0, y * CWIDTH, WIDTH - 1, y * CHEIGHT, WHITE);
//}
//
//void render_vel(const u8 color[4]) {
//	const float csize = min(CWIDTH, CHEIGHT);
//	const float chalf = csize / 2;
//	for (size_t x = 1; x < GWIDTH - 1; x++) {
//		for (size_t y = 1; y < GHEIGHT - 1; y++) {
//			float vx = clamp(VX[y * GWIDTH + x], -chalf, chalf);
//			float vy = clamp(VY[y * GWIDTH + x], -chalf, chalf);
//			i32 x0 = (x + 0.5f) * CWIDTH;
//			i32 y0 = (y + 0.5f) * CHEIGHT;
//			i32 x1 = x0 + vx;
//			i32 y1 = y0 + vy;
//			render_line(x0, y0, x1, y1, color);
//		}
//	}
//}

void render_smoke() {
	for (size_t x = 1; x < GWIDTH - 1; x++) {
		for (size_t y = 1; y < GHEIGHT - 1; y++) {
			float s  = SMOKE[y * GWIDTH + x];
			i32   px = x * CWIDTH;
			i32   py = y * CHEIGHT;
			i32   v  = clamp(s / 2.0f, 0, 255);
			u8 color[4] = { v, v, v, 255 };
			render_rect(px, py, CWIDTH, CHEIGHT, color);
		}
	}
}

void advect_cell(float* tvx, float* tvy, float* ts, double dt, size_t x, size_t y) {
	float  vx = VX[y * GWIDTH + x];
	float  vy = VY[y * GWIDTH + x];
	float  s  = SMOKE[y * GWIDTH + x];
	float  nx = x + vx * dt;
	float  ny = y + vy * dt;
	size_t fx = floor(nx);
	size_t cx = ceil(nx);
	size_t fy = floor(ny);
	size_t cy = ceil(ny);
	float  rx = nx - fx;
	float  ry = ny - fy;
	
	tvx[cy * GWIDTH + cx] += rx * ry * vx;
	tvx[cy * GWIDTH + fx] += (1 - rx) * ry * vx;
	tvx[fy * GWIDTH + cx] += rx * (1 - ry) * vx;
	tvx[fy * GWIDTH + fx] += (1 - rx) * (1 - ry) * vx;
	tvy[cy * GWIDTH + cx] += rx * ry * vy;
	tvy[cy * GWIDTH + fx] += (1 - rx) * ry * vy;
	tvy[fy * GWIDTH + cx] += rx * (1 - ry) * vy;
	tvy[fy * GWIDTH + fx] += (1 - rx) * (1 - ry) * vy;
	ts[cy * GWIDTH + cx]  += rx * ry * s;
	ts[cy * GWIDTH + fx]  += (1 - rx) * ry * s;
	ts[fy * GWIDTH + cx]  += rx * (1 - ry) * s;
	ts[fy * GWIDTH + fx]  += (1 - rx) * (1 - ry) * s;
}

void advect(double dt) {
	float tvx[GWIDTH * GHEIGHT] = { 0 };
	float tvy[GWIDTH * GHEIGHT] = { 0 };
	float ts[GWIDTH  * GHEIGHT] = { 0 };
	for (size_t x = 1; x < GWIDTH - 1; x++)
		for (size_t y = 1; y < GHEIGHT - 1; y++)
			advect_cell(tvx, tvy, ts, dt, x, y);
	memcpy(VX,    tvx, GWIDTH * GHEIGHT * sizeof(float));
	memcpy(VY,    tvy, GWIDTH * GHEIGHT * sizeof(float));
	memcpy(SMOKE, ts,  GWIDTH * GHEIGHT * sizeof(float));
}

void project_cell(double dt, size_t x, size_t y) {
	float vup    = VY[(y + 1) * GWIDTH + (x + 0)];
	float vdown  = VY[(y + 0) * GWIDTH + (x + 0)];
	float vleft  = VX[(y + 0) * GWIDTH + (x + 0)];
	float vright = VX[(y + 0) * GWIDTH + (x + 1)];
	float eup    = !SOLID[(y + 1) * GWIDTH + (x + 0)];
	float edown  = !SOLID[(y - 1) * GWIDTH + (x + 0)];
	float eleft  = !SOLID[(y + 0) * GWIDTH + (x - 1)];
	float eright = !SOLID[(y + 0) * GWIDTH + (x + 1)];
	float v      = vup - vdown - vleft + vright;
	float e      = eup + edown + eleft + eright;
	if (SOLID[y * GWIDTH + x] || e == 0) return;
	float m      = -v / e * 1.5;
	VY[(y + 1) * GWIDTH + (x + 0)] += m * eup;
	VY[(y + 0) * GWIDTH + (x + 0)] -= m * edown;
	VX[(y + 0) * GWIDTH + (x + 0)] -= m * eleft;
	VX[(y + 0) * GWIDTH + (x + 1)] += m * eright;
}

void project(double dt) {
	for (size_t i = 0; i < 50; i++)
		for (size_t x = 1; x < GWIDTH - 1; x++)
			for (size_t y = 1; y < GHEIGHT - 1; y++)
				project_cell(dt, x, y);
}

int main(void) {
	WIN     = RGFW_createWindow("", 0, 0, WIDTH, HEIGHT, RGFW_windowCenter);
	BUF     = malloc(WIDTH * HEIGHT * 4);
	SURFACE = RGFW_createSurface(BUF, WIDTH, HEIGHT, RGFW_formatRGBA8);

	for (size_t x = 0; x < GWIDTH; x++)
		for (size_t y = 0; y < GHEIGHT; y++)
			if (x == 0 || y == 0 || x == GWIDTH - 1 || y == GHEIGHT - 1) {
				SOLID[x] = true;
				SOLID[(GHEIGHT - 1) * GWIDTH + x] = true;
			}

	while (RGFW_window_shouldClose(WIN) == RGFW_FALSE) {
		const float dt = dt_get();
		RGFW_pollEvents();
		if (RGFW_isMouseDown(RGFW_mouseLeft)) {
			float vx, vy;
			RGFW_getMouseVector(&vx, &vy);
			i32 x, y;
			RGFW_window_getMouse(WIN, &x, &y);
			x = (float)x / WIDTH * GWIDTH;
			y = (float)y / HEIGHT * GHEIGHT;
			for (i32 ox = -5; ox <= 5; ox++)
				for (i32 oy = -5; oy <= 5; oy++) {
					if (SOLID[y * GWIDTH + x]) continue;
					VX[(y + oy) * GWIDTH + x + ox] += 100 * vx * dt;
					VY[(y + oy) * GWIDTH + x + ox] += 100 * vy * dt;
					SMOKE[(y + oy) * GWIDTH + x + ox] += 1;
				}
		}

		render_rect(0, 0, WIDTH, HEIGHT, BLACK);
		//render_grid(WHITE);
		render_smoke();
		//render_vel(WHITE);
		//render_vel_split(WHITE);
		RGFW_window_blitSurface(WIN, SURFACE);

		advect(dt);
		project(dt);
	}
}
