###############################################################################
#                                                                             #
#       MAIN.PY-                                                              #
#       THIS FILE CONTAINS THE FILE I/O AND DATA SAVING                       #
#                                                                             #
###############################################################################


# ARGV[1] - INPUT FILE
# ARGV[2] - NUMBER OF TIME STEPS
# ARGV[3] - VERBOSITY (OPTIONAL)

import physics
from sys import argv
import subprocess
import os
import time

def read_input(fname, verbosity = 0):
    with open(fname, "r") as f:

        # First line of file is the dimension of the frame in pixels
        w, h = [int(x) for x in f.readline().strip().split(",")]
        world = physics.World(width=w, height = h, verbosity = verbosity)

        current_shape = ""
        for line in f.readlines():
            line = line.strip()

            # Ignore commented out lines and empty lines
            if not line or line[0] == "#":
                continue
            
            words = line.split()
            if words[0] == "PARAM":
               setattr(physics, words[1], float(words[2]))
               continue
                

            if current_shape == "circle":
                data = line.split(",")
                if data[0] == "mass":
                    mass = float(data[1])
                else:
                    center_x = float(data[0])
                    center_y = float(data[1])
                    rad = float(data[2])
                    obj = physics.Ball(world, pos = [center_x, center_y], radius = rad)
                    current_shape = ""
            
            elif current_shape == "polygon" or current_shape == "fixedpolygon":
                data = line.split(",")
                if data[0] == "mass":
                    mass = float(data[1])
                elif data[0] == "sides":
                    sides = int(data[1])
                    npoints = sides
                    points = []
                elif data[0] == "point":
                    npoints -= 1
                    pp = [float(_) for _ in data[1:]]
                    point = physics.Point(world, mass = mass / sides, pos = pp)
                    points.append(point)
                    if npoints == 0:
                        if current_shape == "polygon":
                            obj = physics.Polygon(world, points = points, mass = mass)
                        else:
                            obj = physics.FixedPolygon(world, points = points)
                        current_shape = ""
            else:
                current_shape = line
                mass = physics.np.inf
    return world


if __name__ == '__main__':
        
    num_iters = int(argv[2])
    verbosity = int(argv[3]) if len(argv) >= 4 else 0
    plane = read_input(argv[1], verbosity)
    processes = []
    for t in range(num_iters):
        plane.update()
        #print(t,", ")
        #for i,obj in enumerate(plane.objs):
        #    print(obj.name,i, obj.pos, obj.vel, obj.acc)
        #    pass
        with open("./simulations/plane_%03d.txt" % t, "w") as f:
            f.write(str(plane))
        if t and t % 126 == 0:
            start = time.time()

            # Does not wait for this to finish running
            p = subprocess.Popen(["./draw", str(t - 126), str(t), "&"] , stdin=None, stdout=None, stderr=None,)
            processes.append(p)
    
    # Does wait for this to finish running
    subprocess.run(["./draw", str(t - (t%126)), str(t)])

    # If any processes are still running, wait for them since ffmpeg will fail 
    # if all pngs are not created 
    for p in processes:
        p.wait()
