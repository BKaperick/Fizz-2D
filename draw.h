struct image {
    int width;
    int height;
    float* pixels;
};
typedef struct image image;

void init_image_empty(image* img, int w, int h);

void init_image(image* img, int w, int h, float* px);

void draw_pixel(image* img, int x, int y);

void draw_circle(image* img, int center_x, int center_y, int radius, bool fill);

void fill_shape(image* img, int x, int y);

