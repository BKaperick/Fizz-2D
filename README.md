Fizz is a simple 2D physics engine designed for use in experimentation with some genetic algorithms to optimize shapes in an evolutionary process.

To compile, run

    $gcc -lpng -o draw draw.c png_util.c

Currently, I ran into some issues with `libpng 1.6.34` but had success installing `libpng 1.2.59` and changing the `LIBPNG_VER_STRING` variable in `png.h` to equal `"1.6.34"`.

To run the simulation, run

    $./simulate.sh [INPUTFILE.in] [NUM_TIMESTEPS]

Note: in `simulate.sh`, you will need to change `py` to the path to your Python 3.x installation.

This script generates `.png` files of the form `plane_%d.png` in the current directory as well as `simul.gif` which is the stitching of all these files into a single GIF.

Here is an example from the latest simulation run.
![Example](simul.gif)
