from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

GRAVITY = np.array([0,-9.8])

class world:
    def __init__(self, objs = set(), global_force = GRAVITY, time_disc = 1, gamma = []):
        self.objs = set()
        self.time_disc = time_disc
        self.state = 0
        
        self.global_force = global_force
        self.global_damping_force = None
        if len(gamma) > 0:
            #self.global_damping_force = lambda v : gamma * v^2
            self.global_damping_force = gamma

        for obj in objs:
            self.init_obj(obj)
    
    def init_obj(self, obj):
        self.objs.add(obj)
        #forces.append(self.global_force)
        if len(self.global_damping_force) > 0:
            obj.damping_force = True

    def update(self):
        for obj in self.objs:
            obj.pre_update([self.global_force], self.time_disc)
#        update = self.check_collisions()
#        for obj1,obj2 in update:
#            self.collide_gracefully(obj1, obj2)
        for obj in self.objs:
            obj.finish_update()
        self.state += 1

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

    
    def pre_update(self, dt, damping_force):
        if damping_force:
            update = self.com.linear_damping_move([world.global_force], dt)
        else:
            update = self.com.move([world.global_force], dt)

        for i in range(len(self.points)):
            self.points[i] += update

    def finish_update(self):
        self.pos = self.new_pos
        

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
        # Update position
        v_avg = self.vel + (.5 * self.acc * dt)
        self.new_pos = self.pos + (v_avg * dt)
        
        # Update acceleration
        self.acc = sum(forces) / self.mass
        
        # Update velocity
        self.vel = v_avg + (.5 * self.acc * dt)

    def linear_damping_move(self, forces, dt):
        self.new_pos = self.pos + (dt * self.vel) + (.5 * (dt**2) * self.acc)
        
        tmp = dt * self.world.global_damping_force / (2 * self.mass)
        Fprev = self.acc * self.mass
        F = sum(forces)
        self.vel = (1 / (1 + tmp)) * (self.vel * (1 -  tmp) + (dt / (2 * self.mass)) * (F + Fprev ))

        self.acc = F / self.mass

    
        



if __name__ == '__main__':
    plane = world(gamma = np.array([0,1]))
    pt = point(plane)
    ball = obj([pt], plane)
    num_iters = int(argv[1])
    
    
    plt.ion()
    fig, ax = plt.subplots()

    plot = ax.scatter([], [])
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)


    for t in range(num_iters):
        plt.figure(1)
        plt.plot(*ball.pos,'-')
        print(ball.pos)

        # get the current points as numpy array with shape  (N, 2)
        array = plot.get_offsets()

        # add the points to the plot
        array = np.append(array, ball.pos)
        plot.set_offsets(array)

        # update x and ylim to show all points:
        ax.set_xlim(array[:, 0].min() - 0.5, array[:,0].max() + 0.5)
        ax.set_ylim(array[:, 1].min() - 0.5, array[:, 1].max() + 0.5)
        # update the figure
        fig.canvas.draw()       
        plane.update()
