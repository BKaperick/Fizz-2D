from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

GRAVITY = np.array([0,.98])

class World:
    def __init__(self, width, height, objs = set(), global_force = GRAVITY, time_disc = 1, gamma = []):
        self.width = width
        self.height = height
        self.objs = set()
        self.fixed_objs = set()
        self.time_disc = time_disc
        self.state = 0
        
        self.global_force = global_force
        self.global_damping_force = np.array([])
        if len(gamma) > 0:
            #self.global_damping_force = lambda v : gamma * v^2
            self.global_damping_force = gamma

        for obj in objs:
            self.init_obj(obj)
    
    def init_obj(self, obj):
        self.objs.add(obj)
        #forces.append(self.global_force)
        if self.global_damping_force:
            obj.damping_force = True

    def init_fixed_obj(self, obj):
        self.fixed_objs.add(obj)

    def update(self):

        # Apply linear damping force if exists
        if len(self.global_damping_force) == 0:
            gdf = []
        else:
            gdf = [self.global_damping_force]
        
        # Do first pass of position, velocity and acceleration updates for each object in the world
        for obj in self.objs:
            obj.pre_update([self.global_force], gdf, self.time_disc)
        
        # Do second (final) pass for each object in the world
        for obj in self.objs:
            obj.finish_update()

        # Advance time
        self.state += self.time_disc

    def __str__(self):
        outstr = str(self.width) + "," + str(self.height) + "\n"
        for obj in self.objs:
            outstr += str(obj)
        for obj in self.fixed_objs:
            outstr += str(obj)
        return outstr


class Obj:
    def __init__(self, points = [], world = None, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        
        self.points = points
        self.com = point(world, pos = sum([pt.pos for pt in points]) / len(points))
        
        self.mass = mass
        
        # Attributes of motion
        self.pos = pos
        self.new_pos = pos
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rot_ang = rotation_angle
        self.rot_spd = rotation_speed
        
        self.world = world
        world.init_obj(self)
    
    def pre_update(self, force, damping_force, dt):
        if damping_force:
            update_x, update_v, new_acc = self.com.linear_damping_move([self.world.global_force], dt)
        
            for i in range(len(self.points)):
                self.points[i].pos += update_x
                self.points[i].vel += update_v
                self.points[i].acc = new_acc

        if force:
            update_x, update_v, new_acc = self.com.move([self.world.global_force], dt)
        
            for i in range(len(self.points)):
                self.points[i].pos += update_x
                self.points[i].vel += update_v
                self.points[i].acc = new_acc

    def finish_update(self):
        self.pos = self.com.new_pos
        self.com.pos = self.com.new_pos
        self.vel = self.com.vel
        self.acc = self.com.acc
        

class Point:
    def __init__(self, world, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0])):

        self.mass = mass
        self.pos = pos
        self.new_pos = pos
        self.vel = speed

        self.acc = np.array([0.0,0.0])
        self.world = world


    def move(self, forces, dt):
        '''
        Velocity Verlet Algorithm to update x,v,a with global error of order 2.
        '''
        update_v = .5 * self.acc * dt
        v_avg = self.vel + update_v

        # Update position
        update_x = v_avg * dt
        self.new_pos = self.pos + update_x
        
        # Update acceleration
        new_acc = sum(forces) / self.mass
        self.acc = new_acc
        
        # Update velocity
        update_v += .5 * self.acc * dt
        self.vel = v_avg + (.5 * self.acc * dt)
        
        return update_x, update_v, new_acc


    def linear_damping_move(self, forces, dt):
        update_x = (dt * self.vel) + (.5 * (dt**2) * self.acc)
        self.new_pos = self.pos + update_x
        
        tmp = dt * self.world.global_damping_force / (2 * self.mass)
        Fprev = self.acc * self.mass
        F = sum(forces)
        update_v = (1 / (1 + tmp)) * (self.vel * (1 -  tmp) + (dt / (2 * self.mass)) * (F + Fprev )) - self.vel
        self.vel += update_v
        
        new_acc = F / self.mass
        self.acc = F / self.mass

        return update_x, update_v, new_acc


class Polygon(Obj):
    def __init__(self, world = None, mass = 1, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        
        self.points = points
        self.com = Point(world, pos = sum([pt.pos for pt in points]) / len(points))
        
        self.mass = mass
        
        # Attributes of motion
        self.pos = self.com.pos
        self.new_pos = self.pos
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rot_ang = rotation_angle
        self.rot_spd = rotation_speed
        
        self.world = world
        if world:
            world.init_obj(self)
    
    def __str__(self):
        # Note we leave out mass here since that is irrelevant to graphical representation at a fixed point in time.
        return "polygon\nsides," + str(len(self.points))+ "\n" + "\n".join(["point," + str(int(pt.pos[0])) + "," + str(int(pt.pos[1])) for pt in self.points]) + "\n"


class Ball(Obj):
    def __init__(self, world = None, mass = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        
        self.points = []
        self.com = Point(world, pos = pos)
        self.radius = radius
        
        self.mass = mass
        
        # Attributes of motion
        self.pos = pos
        self.new_pos = pos
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rot_ang = rotation_angle
        self.rot_spd = rotation_speed
        
        self.world = world
        if world:
            world.init_obj(self)

    def __str__(self):
        return "circle\n" + str(int(self.pos[0])) + "," + str(int(self.pos[1])) + "," + str(int(self.radius)) + "\n"



class FixedPolygon(Polygon):
    def __init__(self, world = None, points = []):
        
        self.points = []
        self.com = Point(world, pos = sum([pt.pos for pt in self.points] / len(self.points)))
        
        self.mass = np.inf
        
        # Attributes of motion
        self.pos = pos
        self.new_pos = pos
        self.vel = 0
        self.acc = np.array([0.0,0.0])
        
        self.world = world
        if world:
            world.init_fixed_obj(self)
        
def read_input(fname):
    with open(fname, "r") as f:
        w, h = [int(x) for x in f.readline().strip().split(",")]
        world = World(width=w, height = h)
        current_shape = ""
        for line in f.readlines():
            if current_shape == "circle":
                data = line.strip().split(",")
                if data[0] == "mass":
                    mass = float(data[1])
                else:
                    center_x = float(data[0])
                    center_y = float(data[1])
                    rad = float(data[2])
                    obj = Ball(world, pos = np.array([center_x, center_y]), radius = rad)
                    world.init_obj(obj)
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
                    pp = np.array([float(_) for _ in data[1:]])
                    point = Point(world, mass = mass / sides, pos = pp)
                    points.append(point)
                    if npoints == 0:
                        if current_shape == "polygon":
                            obj = Polygon(world, points = points, mass = mass)
                        else:
                            obj = FixedPolygon(world, points = points)
                        world.init_obj(obj)
                        current_shape = ""
            else:
                current_shape = line.strip()
    return world
            

if __name__ == '__main__':
    
#    plane = World(width=100, height = 500)#gamma = np.array([0,1]))
#    ball = Ball(plane, pos=np.array([50.0, 100.0]), radius = 15)
#    floor = Polygon(plane, points = [
#        Point(plane, pos=np.array([0.0, plane.height])),
#        Point(plane, pos = np.array([100.0, plane.height])),
#        Point(plane, pos = np.array([100.0, plane.height-1])),
#        Point(plane, pos = np.array([0.0, plane.height-1]))])
        
    plane = read_input(argv[1])
    num_iters = int(argv[2])
    for t in range(num_iters):
        plane.update()
        for obj in plane.objs:
            print(obj.pos)
        with open("plane_{0}.txt".format(t), "w") as f:
            f.write(str(plane))
