from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

GRAVITY = np.array([0,.98])

class World:
    def __init__(self, width, height, objs = set(), global_force = GRAVITY, time_disc = 1, gamma = []):
        self.objs = set()
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


class obj:
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
        

class Ball(obj):
    def __init__(self, world = None, mass = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0):
        
        self.points = []
        self.com = point(world, pos = pos)
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
        


class point:
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

    
        



if __name__ == '__main__':
    plane = World(width=100, height = 500)#gamma = np.array([0,1]))
    pt = point(plane)
    #ball = obj([pt], plane, pos=np.array([0.0, 100]))
    ball = Ball(plane, pos=np.array([50.0, 100.0]), radius = 15)
    num_iters = int(argv[1])
    for t in range(num_iters):
        #print(ball.pos)
        with open("plane_{0}.txt".format(t), "w") as f:
            f.write("100,500\n")
            f.write("circle\n{0},{1},{2}\n".format(int(ball.pos[0]), int(ball.pos[1]), int(ball.radius)))
        plane.update()
        
