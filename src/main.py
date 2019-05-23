###############################################################################
#                                                                             #
#       MAIN.PY-                                                              #
#       THIS FILE CONTAINS THE FILE I/O AND DATA SAVING                       #
#                                                                             #
###############################################################################


# ARGV[1] - INPUT FILE
# ARGV[2] - NUMBER OF TIME STEPS
# ARGV[3] - VERBOSITY (OPTIONAL)
# ARGV[4] - ENERGY LOG (OPTIONAL)

import physics
from sys import argv
import subprocess
import os
import time
from plot import plot_energy
import config
from config import SIMULATION_DIR

def read_input(fname, energylog = 0):
    with open(fname, "r") as f:

        # First line of file is the dimension of the frame in pixels
        w, h = [int(x) for x in f.readline().strip().split(",")]
        world = physics.World(width=w, height = h)

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
                if data[0] == "density":
                    density = float(data[1])
                else:
                    center_x = float(data[0])
                    center_y = float(data[1])
                    rad = float(data[2])
                    obj = physics.Ball(world, pos = [center_x, center_y], radius = rad)
                    current_shape = ""
            
            elif current_shape == "polygon" or current_shape == "fixedpolygon":
                data = line.split(",")
                if data[0] == "density":
                    density = float(data[1])
                elif data[0] == "sides":
                    sides = int(data[1])
                    npoints = sides
                    points = []
                elif data[0] == "velocity":
                    vel = physics.np.array([float(data[1]), float(data[2])])
                elif data[0] == "point":
                    npoints -= 1
                    pp = [float(_) for _ in data[1:]]
                    point = physics.Point(world, pos = pp, speed = vel)
                    points.append(point)
                    if npoints == 0:
                        if current_shape == "polygon":
                            if density == physics.np.inf:
                                print("Warning: density never set for polygon, set to 1")
                                density = 1
                            obj = physics.Polygon(world, points = points, density = density, speed = vel)
                        else:
                            obj = physics.FixedPolygon(world, points = points)
                        current_shape = ""
            else:
                vel = physics.np.array([0.0,0.0])
                current_shape = line
                density = physics.np.inf
    return world


if __name__ == '__main__':
        
    num_iters = int(argv[2])

    # Override defaults set in config.py
    physics.verbosity = int(argv[3]) if len(argv) >= 4 else config.VERBOSITY
    physics.energylog = int(argv[4]) if len(argv) >= 5 else config.ENERGYLOG
    plane = read_input(argv[1])
    processes = []
    for t in range(num_iters):
        plane.update()
        #print(t,", ")
        #for i,obj in enumerate(plane.objs):
        #    print(obj.name,i, obj.pos, obj.vel, obj.acc)
        #    pass
        with open(SIMULATION_DIR + "plane_%03d.txt" % t, "w") as f:
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
    
    # Lastly, display energy
    if physics.energylog:
        plane.log.close()
        plot_energy(SIMULATION_DIR+'energy.txt')
