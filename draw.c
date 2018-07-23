
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include "png_util.h"
#include "draw.h"

#define FILLED 1.0
#define EMPTY 0.0

image* init_image_empty(int w, int h) {
    image* img = malloc(sizeof(image));
    img->width = w;
    img->height = h;
    img->pixels = calloc(w * h, sizeof(float));
    return img;
}    

image* init_image(int w, int h, float* px) {
    image* img = malloc(sizeof(image));
    img->width = w;
    img->height = h;
    img->pixels = px;
    return img;
}    

void draw_pixel(image* img, int x, int y) {
    if (x >= 0 && x < img->width && y >= 0 && y < img->height) {
        img->pixels[y*(img->width) + x] = FILLED;
    }
}

float get_pixel(image* img, int x, int y) {
    if (x < 0 || x >= img->width || y < 0 || y >= img->height) {
        return -1.0;
    }
    return img->pixels[y*(img->width) + x];
}

void draw_circle(image* img, int center_x, int center_y, int radius, bool fill) {
    /*
     * A near-direct copy-paste from the Wikipedia implementation of the Midpoint
     * Circle Algorithm.  Idea is to start at 0 degrees and choose pixel path
     * to 45 degrees which maximizes x^2 + y^2 without exceeding r^2.  Then 
     * symmetry is exploited to duplicate this path 8 times.
     *
     */

    int x = radius-1;
    int y = 0;
    int dx = 1;
    int dy = 1;
    int err = dx - (radius << 1);
    printf("center (%d, %d)\n", center_x, center_y);
    while (x >= y)
    {
        draw_pixel(img, center_x + x, center_y + y);
        draw_pixel(img, center_x + y, center_y + x);
        draw_pixel(img, center_x - y, center_y + x);
        draw_pixel(img, center_x - x, center_y + y);
        draw_pixel(img, center_x - x, center_y - y);
        draw_pixel(img, center_x - y, center_y - x);
        draw_pixel(img, center_x + y, center_y - x);
        draw_pixel(img, center_x + x, center_y - y);

        if (err <= 0)
        {
            y++;
            err += dy;
            dy += 2;
        }
        if (err > 0)
        {
            x--;
            dx += 2;
            err += (-radius << 1) + dx;
        }
    }
    
    if (fill) {
        fill_shape(img, center_x, center_y);
    }
}

void draw_line(image* img, int x0, int y0, int x1, int y1) {
    
    // Trivial one-pixel line
    if ((x0 == x1) && (y0 == y1)) {
        draw_pixel(img, x0, y0);
        return;
    }
    int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
    int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1; 
    int err = dx+dy, e2; /* error value e_xy */

    if (dx == 0) {
        for (int i = 0; i < abs(dy); i++) {
            draw_pixel(img, x0, y0 + (sy*i));
        }
        return;
    }
    else if (dy == 0) {
        for (int i = 0; i < abs(dx); i++) {
            draw_pixel(img, x0 + (sx*i), y0);
        }
        return;
    }
 
    for(;;){  /* loop */
        draw_pixel(img, x0, y0);
        if (x0==x1 && y0==y1) break;
        e2 = 2*err;
        if (e2 >= dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
        if (e2 <= dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
    }
}

void draw_polygon(image* img, int* points_x, int*points_y, int num_vertices, bool fill) {
    int x0,y0,x1,y1;
    
    // Draw the border one line at a time
    for (int ind = 1; ind <= num_vertices; ind++) {
        
        // Fetch current endpoints
        x0 = *(points_x+ind-1);
        y0 = *(points_y+ind-1);
        x1 = *(points_x+(ind%num_vertices));
        y1 = *(points_y+(ind%num_vertices));
        
        // Fill in the line pixel by pixel
        draw_line(img, x0, y0, x1, y1);
    }        
    
    if (fill) {
        float center_x = 0;
        float center_y = 0;
        int next_x, next_y;
        for (int i = 0; i < num_vertices; i++) {

            // These checks should make the fill-in work when the center of the 
            // object is outside the bounds of img.
            next_x = *(points_x + i);
            if (next_x < 0)
                next_x = -1;
            if (next_x >= img->width)
                next_x = img->width;
            next_y = *(points_y + i);
            if (next_y < 0)
                next_y = -1;
            if (next_y >= img->height)
                next_y = img->height;

            center_x = center_x + next_x;
            center_y = center_y + next_y;
        }
        int center_px_x = (int) (center_x / num_vertices);
        int center_px_y = (int) (center_y / num_vertices);
        fill_shape(img, center_px_x, center_px_y);
    }
}

void fill_shape(image* img, int x, int y) {
    /*
     * Classic recursive implementation of flood fill.
     * If this ends up being too slow, there are some
     * easy optimizations to be made.
     */

    if (get_pixel(img, x, y) == EMPTY) {
        draw_pixel(img, x, y);
    }
    else return;
    if (get_pixel(img, x+1, y) == EMPTY) {
        fill_shape(img, x+1, y);
    }
    if (get_pixel(img, x-1, y) == EMPTY) {
        fill_shape(img, x-1, y);
    }
    if (get_pixel(img, x, y+1) == EMPTY) {
        fill_shape(img, x, y+1);
    }
    if (get_pixel(img, x, y-1) == EMPTY) {
        fill_shape(img, x, y-1);
    }
}

image* spec_to_image(char* fname) {
    FILE* fp = fopen(fname, "r");
    int width, height, x, y, rad;
    char* shape = malloc(40 * sizeof(char));

    fscanf(fp, "%d,%d\n", &width, &height);
    image* img = init_image_empty(width, height);
    
    
    
    while (!feof(fp)) {
        fgets(shape, 10, fp);
        if (strcmp(shape, "circle\n") == 0) {
            fscanf(fp, "%d,%d,%d\n", &x, &y, &rad);
            draw_circle(img, x, y, rad, 0);
        }
        else if (strcmp(shape, "polygon\n") == 0) {
            int sides;
            fscanf(fp, "sides,%d\n", &sides);
            int* xs = malloc(sides * sizeof(int));
            int* ys = malloc(sides * sizeof(int));
            for (int i = 0; i < sides; i++) {
                fscanf(fp, "point,%d,%d\n", &(xs[i]), &(ys[i]));
            }
            draw_polygon(img, xs, ys, sides, 1);
        }
        else
            return img;
    }
    return img;
}

void save_image(image* img, char* fname) {
    float minI = EMPTY;
    float maxI = FILLED;
    FILE* fp = fopen(fname, "w");
    write_gray_png(fp, img->width, img->height, img->pixels, minI, maxI);
}

int main(int argc, char* argv[]) {
    uint16_t num_files_start = atoi(argv[1]);
    uint16_t num_files_end = atoi(argv[2]);
    uint16_t num_files = num_files_end - num_files_start;
    
    char* fname_in;
    char* fname_out;
    image* img;
    
    fname_in = malloc(30 * sizeof(char));
    fname_out = malloc(30 * sizeof(char));
    for (uint16_t ind = num_files_start; ind < num_files_end; ind++) {
        sprintf(fname_in, "./simulations/plane_%03"PRId16".txt", ind);
        sprintf(fname_out, "./simulations/plane_%03"PRId16".png", ind);
        img = spec_to_image(fname_in);
        save_image(img, fname_out);
    }
    printf("processed images %03"PRId16" through %03"PRId16"\n", num_files_start, num_files_end);
    free(fname_in);
    free(fname_out);

	return 0;
}

