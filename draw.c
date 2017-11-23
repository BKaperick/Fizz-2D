
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include "png_util.h"
#include "draw.h"

#define FILLED 1.0
#define EMPTY 0.0

void init_image_empty(image* img, int w, int h) {
    img->width = w;
    img->height = h;
    img->pixels = calloc(w * h, sizeof(float));
}    

void init_image(image* img, int w, int h, float* px) {
    img->width = w;
    img->height = h;
    img->pixels = px;

}    

void draw_pixel(image* img, int x, int y) {
    img->pixels[y*(img->width) + x] = FILLED;
}

float get_pixel(image* img, int x, int y) {
    if (x < 0 || x > img->width || y < 0 || y > img->height) {
        return -1.0;
    }
    return img->pixels[y*(img->width) + x];

}

void draw_circle(image* img, int center_x, int center_y, int radius, bool fill) {
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


int main() {
    FILE* fp = fopen("test.png", "w");
    int width = 64;
    int height = 64;
    
    image img;// = malloc(sizeof(image));
    init_image_empty(&img, width, height);
    draw_circle(&img, 32, 32, 12, 1);


//    float *img = malloc(width * height * sizeof(float));
//    for (int i = 0; i < width * height; i++) {
//        img[i] = (float) i / (width * height);
//    }
    float minI = EMPTY;
    float maxI = FILLED;
    write_gray_png(fp, img.width, img.height, img.pixels, minI, maxI);
	return 0;
}

