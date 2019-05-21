Fizz is a simple 2D physics engine designed for use in experimentation with some genetic algorithms to optimize shapes in an evolutionary process.

To compile, from the parent directory, run

    $gcc -lpng -o draw ./src/draw.c ./src/png_util.c

Currently, I ran into some issues with `libpng 1.6.34` but had success installing `libpng 1.2.59` and changing the `LIBPNG_VER_STRING` variable in `png.h` to equal `"1.6.34"`.

To run the simulation, run

    $./simulate.sh [INPUTFILE] [NUM_TIMESTEPS] [VERBOSITY]

Note: [VERBOSITY] is a nonnegative integer indicating how much debugging info you wish to be displayed.  VERBOSITY = 0 does not print anything.  If you do not include VERBOSITY, it defaults to 0.

This script generates `.png` files of the form `plane_%d.png` in the `simulations` directory as well as `simul.gif`, which is the stitching of all these files into a single GIF.

Here are some examples from the latest simulation run.

![Example 1](./simulations/simul.gif)
![Example 2](./simulations/simul2.gif)
![Example 3](./simulations/simul3.gif)
