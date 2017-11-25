from sys import argv
import numpy as np
from math import acos
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

GRAVITY = np.array([0.1,.2])

class World:
    def __init__(self, width, height, objs = set(), global_accel = GRAVITY, time_disc = 1, gamma = []):
        self.width = width
        self.height = height
        self.objs = []
        self.fixed_objs = []
        self.time_disc = time_disc
        self.state = 0
        
        self.global_accel = global_accel
        self.global_damping_force = np.array([])
        if len(gamma) > 0:
            #self.global_damping_force = lambda v : gamma * v^2
            self.global_damping_force = gamma

        for obj in objs:
            self.init_obj(obj)
    
    def init_obj(self, obj):
        self.objs.append(obj)
        #forces.append(self.global_force)
        if self.global_damping_force:
            obj.damping_force = True

    def init_fixed_obj(self, obj):
        self.fixed_objs.append(obj)

    def update(self):

        # Apply linear damping force if exists
        if len(self.global_damping_force) == 0:
            gdf = []
        else:
            gdf = [self.global_damping_force]
        
        # Do first pass of position, velocity and acceleration updates for each object in the world
        for obj in self.objs:
            obj.pre_update([], gdf, self.time_disc)
        
        collisions = self.check_collisions()
        while collisions:
            for obj1, obj2, c_vec_normal, c_vec_mag, f1, f2 in collisions:
                if obj1.is_fixed:
                    obj1,obj2 = obj2,obj1
                    #c_vec_normal = -1 * c_vec_normal
                c_vec = c_vec_normal * (c_vec_mag + 1e-7)
                print("COLLISION", c_vec, f1, f2)
                if f2 == 1:
                    c_vec *= -1
                    c_vec_normal *= -1
#                
                mass_prop = obj2.mass / (obj1.mass + obj2.mass) 
                if np.isnan(mass_prop):
                    mass_prop = 0
                
                if obj2.is_fixed:
                    #self.forces.append(obj1.mass * 
                    for i in range(len(obj1.points)):
                        #print("updates: ", c_vec, np.dot(obj1.points[i].vel, c_vec))
                        obj1.points[i].pos += c_vec
                        
                        vdotN = np.dot(obj1.points[i].vel, c_vec_normal)
                        obj1.points[i].vel -= 2*vdotN * c_vec_normal
                        
                        #adotN = np.dot(obj1.points[i].acc, c_vec_normal)
                        #obj1.points[i].acc -= 2*adotN * c_vec_normal

                    obj1.com.new_pos += c_vec
                    
                    vdotN = np.dot(obj1.com.vel, c_vec_normal)
                    obj1.com.vel -= 2*vdotN * c_vec_normal
                    
                    #adotN = np.dot(obj1.com.acc, c_vec_normal)
                    #obj1.com.acc -= 2*adotN * c_vec_normal
                    
                    obj1.finish_update()
                    obj2.finish_update()

                else:
                    obj1.com.new_pos += (1-mass_prop) * (c_vec)
                    for i in range(len(obj1.points)):
                        obj1.points[i].pos += (1-mass_prop) * (c_vec)
                        
                    obj1.finish_update()
                    
                    obj2.com.new_pos -= mass_prop * c_vec
                    for i in range(len(obj2.points)):
                        obj2.points[i].pos -= mass_prop * (c_vec)
                    obj2.finish_update()
            collisions = self.check_collisions()
            #exit()

        # Do second (final) pass for each object in the world
        for obj in self.objs:
            obj.finish_update()

        # Advance time
        self.state += self.time_disc

    def check_collisions(self):
        collisions = []
        objects = [o for o in self.objs + self.fixed_objs if o.name == "polygon" or o.name == "fixedpolygon"]
        for i,obj in enumerate(objects):
            for other_obj in objects[i+1:]:
                correction, correct_mag, f1, f2 = polypoly_collision(obj, other_obj)
                if len(correction) > 0:
                    collisions.append((obj, other_obj, correction, correct_mag, f1, f2))
        return collisions

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
        self.com = Point(world, pos = sum([pt.pos for pt in points]) / len(points), mass = mass)
        
        self.mass = mass
        
        # Attributes of motion
        self.pos = pos
        self.new_pos = pos
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rot_ang = rotation_angle
        self.rot_spd = rotation_speed

        self.update_x = 0.0
        self.update_v = 0.0
        self.new_acc = 0.0
        
        self.world = world
        if world:
            world.init_obj(self)

        self.is_fixed = False
    
    def pre_update(self, force, damping_force, dt):
        if damping_force:
            update_x, update_v, new_acc = self.com.linear_damping_move([], dt)
        
#            for i in range(len(self.points)):
#                self.points[i].pos += update_x
#                self.points[i].vel += update_v
#                self.points[i].acc = new_acc
        self.update_x, self.update_v, self.new_acc = self.com.move(force, dt)
        for i in range(len(self.points)):
            self.points[i].pos += self.update_x
            self.points[i].vel += self.update_v
            self.points[i].acc = np.array(self.new_acc, copy=True)

    def finish_update(self):
        self.pos = np.array(self.com.new_pos, copy=True)
        self.com.pos = np.array(self.com.new_pos, copy=True)
        self.vel = np.array(self.com.vel, copy=True)
        self.acc = np.array(self.com.acc, copy=True)
#        for i in range(len(self.points)):
#            self.points[i].pos += self.update_x
#            self.points[i].vel += self.update_v
#            self.points[i].acc = self.new_acc
        

class Point:
    def __init__(self, world, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0])):
        
        self.name = "point"

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
        new_acc = (sum(forces) / self.mass) + self.world.global_accel
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
        
        super().__init__(world = world, mass = mass, points = points, speed = speed, rotation_angle = rotation_angle, rotation_speed = rotation_speed)
        self.name = "polygon"
        
    
    def __str__(self):
        # Note we leave out mass here since that is irrelevant to graphical representation at a fixed point in time.
        return "polygon\nsides," + str(len(self.points))+ "\n" + "\n".join(["point," + str(int(pt.pos[0])) + "," + str(int(pt.pos[1])) for pt in self.points]) + "\n"


class Ball(Obj):
    def __init__(self, world = None, mass = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        
        self.name = "circle"

        self.points = []
        self.com = Point(world, pos = pos, mass = mass)
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
    def __init__(self, world = None, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        super().__init__(world = None, mass = np.inf, points = points, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
        self.name = "fixedpolygon"        
        self.world = world
        if world:
            world.init_fixed_obj(self)
        self.is_fixed = True
        
def read_input(fname):
    with open(fname, "r") as f:
        w, h = [int(x) for x in f.readline().strip().split(",")]
        world = World(width=w, height = h)
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
                    obj = Ball(world, pos = np.array([center_x, center_y]), radius = rad)
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
                        current_shape = ""
            else:
                current_shape = line.strip()
                mass = np.inf
    return world

def polypoly_collision(poly1, poly2):
    poly1_edges = list(zip(poly1.points[:-1], poly1.points[1:])) + [(poly1.points[-1], poly1.points[0])]
    poly1_normals = [np.array([p2.pos[1] - p1.pos[1], p1.pos[0] - p2.pos[0]]) for p1,p2 in poly1_edges]
    
    poly2_edges = list(zip(poly2.points[:-1], poly2.points[1:])) + [(poly2.points[-1], poly2.points[0])]
    poly2_normals = [np.array([p2.pos[1] - p1.pos[1], p1.pos[0] - p2.pos[0]]) for p1,p2 in poly2_edges]
    

    overlap = np.inf
    for axis,flag in list(zip(poly1_normals,[1]*len(poly1_normals))) + list(zip(poly2_normals, [2]*len(poly2_normals))):
        
        # Normalize
        normal_axis = axis / np.linalg.norm(axis)
        
        poly1_bounds = [np.dot(normal_axis, poly1.points[0].pos)]*2
        poly2_bounds = [np.dot(normal_axis, poly2.points[0].pos)]*2

        for point in poly1.points:
            curr = np.dot(normal_axis, point.pos)
            if curr < poly1_bounds[0]:
                poly1_bounds[0] = curr
            if curr > poly1_bounds[1]:
                poly1_bounds[1] = curr
        
        for point in poly2.points:
            curr = np.dot(normal_axis, point.pos)
            if curr < poly2_bounds[0]:
                poly2_bounds[0] = curr
            if curr > poly2_bounds[1]:
                poly2_bounds[1] = curr
        
        mag_flag = 1
        
        if poly1_bounds[0] > poly2_bounds[0]:
            mag_flag = 2
            poly1_bounds,poly2_bounds = poly2_bounds, poly1_bounds
            
        #If a single check fails, then there is no collision
        if poly1_bounds[1] < poly2_bounds[0]:
            return [],None,None,None
        
        # Save axis and overlap magnitude if this is the smallest gap
        current_overlap = poly1_bounds[1] - poly2_bounds[0]
        if current_overlap < overlap:
            fix_axis = normal_axis
            fix_flag = flag
            fix_mag_flag = mag_flag
            overlap = current_overlap

    #bounceback_displacement = fix_axis * overlap
    return fix_axis, overlap, fix_flag, fix_mag_flag
    #return bounceback_displacement


if __name__ == '__main__':
        
    plane = read_input(argv[1])
    num_iters = int(argv[2])
    for t in range(num_iters):
        plane.update()
        for i,obj in enumerate(plane.objs):
            print(obj.name,i, obj.pos, obj.vel, obj.acc)
        with open("plane_{0}.txt".format(t), "w") as f:
            f.write(str(plane))
