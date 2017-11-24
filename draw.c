
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
        fgets(shape, 7, fp);
        if (strcmp(shape, "circle")) {
            fscanf(fp, "%d,%d,%d", &x, &y, &rad);
            draw_circle(img, x, y, rad, 1);
        }
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
    int num_files = atoi(argv[1]);
    //printf("%d\n", num_files);
    char* fname_in = malloc(11 * sizeof(char));
    char* fname_out = malloc(11 * sizeof(char));
    image* img;    
    for (int ind = 0; ind < num_files; ind++) {
        sprintf(fname_in, "plane_%d.txt", ind);
        sprintf(fname_out, "plane_%d.png", ind);
        img = spec_to_image(fname_in);
        save_image(img, fname_out);
    }
	return 0;
}

