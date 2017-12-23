###############################################################################
#                                                                             #
#       MAIN.PY-                                                              #
#       THIS FILE CONTAINS THE FILE I/O AND DATA SAVING                       #
#                                                                             #
###############################################################################


import physics
from sys import argv

def read_input(fname):
    with open(fname, "r") as f:
        w, h = [int(x) for x in f.readline().strip().split(",")]
        world = physics.World(width=w, height = h)
        current_shape = ""
        for line in f.readlines():
            if line.strip()[0] == "#":
                continue
            if current_shape == "circle":
                data = line.strip().split(",")
                if data[0] == "mass":
                    mass = float(data[1])
                else:
                    center_x = float(data[0])
                    center_y = float(data[1])
                    rad = float(data[2])
                    obj = physics.Ball(world, pos = [center_x, center_y], radius = rad)
                    current_shape = ""
            
            elif current_shape == "polygon" or current_shape == "fixedpolygon":
                data = line.strip().split(",")
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
                current_shape = line.strip()
                mass = physics.np.inf
    return world


if __name__ == '__main__':
        
    plane = read_input(argv[1])
    num_iters = int(argv[2])
    for t in range(num_iters):
        plane.update()
        for i,obj in enumerate(plane.objs):
            print(obj.name,i, obj.pos, obj.vel, obj.acc)
        with open("plane_{0}.txt".format(t), "w") as f:
            f.write(str(plane))
