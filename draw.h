struct image {
    int width;
    int height;
    float* pixels;
};
typedef struct image image;

image* init_image_empty(int w, int h);

image* init_image(int w, int h, float* px);

void draw_pixel(image* img, int x, int y);

void draw_circle(image* img, int center_x, int center_y, int radius, bool fill);

void fill_shape(image* img, int x, int y);

void save_img(image* img, char* fname);

image* spec_to_image(char* fname);

void draw_line(image* img, int x0, int y0, int x1, int y1);

void draw_polygon(image* img, int* points_x, int*points_y, int num_vertices, bool fill);
