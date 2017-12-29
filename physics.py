###############################################################################
#                                                                             #
#       PHYSICS.PY-                                                           #
#       THIS FILE CONTAINS THE CLASSES FOR ALL OBJECTS IN THE WORLD           #
#       AND THEIR PHYSICAL PROPERTIES                                         #
#                                                                             #
###############################################################################

import numpy as np
from   math import acos

# World constants
GRAVITY = np.array([0.1,.2])
EPSILON = 1e-7
ELASTICITY = 0.5
COLLISION_TOL = 1

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
        '''
        This function is called once each time step.  All position and velocity
        updates, and all collisions are handled within.
        '''
        
        dt = self.time_disc
        max_collision_overlap = COLLISION_TOL + 1
        first_iter = True
        while max_collision_overlap > COLLISION_TOL:       
            if first_iter:
                first_iter = False
            else:
                print("collision")
                dt /= 2
                for obj in self.objs:
                    obj.reverse_update()
            
            # Apply linear damping force if exists
            if len(self.global_damping_force) == 0:
                gdf = []
            else:
                gdf = [self.global_damping_force]
            
            # Do first pass of position, velocity and acceleration updates for each object in the world
            for obj in self.objs:
                obj.pre_update([], gdf, dt)

            # check_collisions() compiles a list of collisions with information:
            #   (1,2) the two objects contained in the collision
            #   (3) the unit vector normal to the surface of collision
            #   (4) magnitude of normal vector
            #   (5) a flag indicating a sign change for the normal vector
            collisions = self.check_collisions()
            magnitudes = [abs(magnitude) for _,__,___,magnitude,____ in collisions]
            #print(len(collisions), max([0] + magnitudes), dt)
            if collisions:
                max_collision_overlap = max(magnitudes)
            else:
                max_collision_overlap = 0
            
        
        while collisions:
            for obj1, obj2, unit_normal_vec, normal_vec_mag, flag in collisions:
                
                # Makes things easier if the first object is never fixed
                if obj1.is_fixed:
                    obj1,obj2 = obj2,obj1
                
                # Reassemble normal vector with some added spacing so two 
                # objects are definitely non-overlapping after collision 
                # resolution
                normal_vec = unit_normal_vec * (normal_vec_mag + EPSILON)
                
                # If flag, the normal vector needs to be flipped
                if flag == 1:
                    normal_vec *= -1
                    unit_normal_vec *= -1
                
                mass_prop = obj2.mass / (obj1.mass + obj2.mass) 
                if np.isnan(mass_prop):
                    mass_prop = 0
                
                # obj1 is NEVER fixed, so this branch is for a collision with a 
                # free object and a fixed object
                #
                # 1. The free object updates its position by the normal vector
                # 2. The free object updates its velocity by twice the (negated)
                #    component of its velocity normal to the surface
                #
                if obj2.is_fixed:
                    
                    for i in range(len(obj1.points)):
                        obj1.points[i].pos += normal_vec
                        
                        vdotN = np.dot(obj1.points[i].vel, unit_normal_vec)
                        obj1.points[i].vel -= 2*vdotN * unit_normal_vec

                    obj1.com.pos += normal_vec
                    
                    vdotN = np.dot(obj1.com.vel, unit_normal_vec)
                    obj1.com.vel -= 2*vdotN * unit_normal_vec
                    
                    obj1.finish_update()
                    obj2.finish_update()

                # Only other possible case is that both objects are free
                else:
                    obj1.com.pos += (1-mass_prop) * (normal_vec)
                    for point in obj1.points:
                        point.pos += (1-mass_prop) * (normal_vec)
                        
                    obj1.finish_update()
                    
                    obj2.com.pos -= mass_prop * normal_vec
                    for point in obj2.points:
                        point.pos -= mass_prop * (normal_vec)
                    obj2.finish_update()

            collisions = self.check_collisions()

        # Do second (final) pass for each object in the world
        for obj in self.objs:
            obj.finish_update()

        # Advance time
        print("dt=", self.time_disc, dt)
        self.state += dt#self.time_disc

    def check_collisions(self):
        '''
        Iterates through each pair of objects in the world, and determines if any are currently overlapping.
        Currently, only polygons and fixed polygons are supported.
        TODO: Implement collision detection for circles.
        '''
        collisions = []
        objects = [o for o in self.objs + self.fixed_objs if o.name == "polygon" or o.name == "fixedpolygon"]
        for i,obj in enumerate(objects):
            for other_obj in objects[i+1:]:
                correction, correct_mag, flag = polypoly_collision(obj, other_obj)
                if len(correction) > 0:
                    collisions.append((obj, other_obj, correction, correct_mag, flag))
        return collisions

    def __str__(self):
        '''
        Output the current state of the world as a string to be read by png 
        script.
        '''
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
        '''
        All points in the object are updated one step in accordance with the 
        calculations done on the center of mass point.
        The object's official position, velocity and acceleration are not updated
        until self.finish_update().
        '''
        if damping_force:
            update_x, update_v, new_acc = self.com.linear_damping_move([], dt)
        
#            for i in range(len(self.points)):
#                self.points[i].pos += update_x
#                self.points[i].vel += update_v
#                self.points[i].acc = new_acc
        
        self.update_x, self.update_v, self.new_acc = self.com.move(force, dt)
        for point in self.points:

            # Retain copy of previous step
            point.oldpos = np.copy(point.pos)
            point.oldvel = np.copy(point.vel)
            point.oldacc = np.copy(point.acc)
            
            point.update_with_object(self, self.update_x, self.update_v, self.new_acc)

    def reverse_update(self):
        '''
        Can be used to take one step backward in time, since each point stores one previous
        position, velocity and acceleration.  
        Currently, this is used in collision resolution to test out various time steps before
        choosing one sufficiently close to the time of collision.
        '''
        self.com.pos = self.com.oldpos
        self.com.vel = self.com.oldvel
        self.com.acc = self.com.oldacc

        for point in self.points:
            point.pos = point.oldpos
            point.vel = point.oldvel
            point.acc = point.oldacc 

    def finish_update(self):
        '''
        The object's position, velocity and acceleration are updated to be 
        consistent with its center of mass.
        '''
        self.pos = np.array(self.com.pos, copy=True)
        self.vel = np.array(self.com.vel, copy=True)
        self.acc = np.array(self.com.acc, copy=True)
        

class Point:
    def __init__(self, world, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0])):
        
        self.name = "point"

        self.mass = mass

        self.oldpos = np.array(pos, dtype=float)
        self.pos = np.array(pos, dtype=float)

        self.oldvel = speed
        self.vel = speed
        
        self.oldacc = np.array([0.0,0.0])
        self.acc = np.array([0.0,0.0])

        self.world = world


    def move(self, forces, dt):
        '''
        Velocity Verlet Algorithm to update x,v,a with global error of order 2.
        '''

        halfdt = .5*dt
        update_v = self.acc * halfdt
        v_avg = self.vel + update_v
        
        # Update position
        update_x = v_avg * dt
        self.pos += update_x
        
        # Update acceleration
        new_acc = (sum(forces) / self.mass) + self.world.global_accel
        self.acc = new_acc
        
        # Update velocity
        update_v += self.acc * halfdt
        self.vel = v_avg + (self.acc * halfdt)
        
        return update_x, update_v, new_acc


    def linear_damping_move(self, forces, dt):
        update_x = (dt * self.vel) + (.5 * (dt**2) * self.acc)
        self.pos += update_x
        
        tmp = dt * self.world.global_damping_force / (2 * self.mass)
        Fprev = self.acc * self.mass
        F = sum(forces)
        update_v = (1 / (1 + tmp)) * (self.vel * (1 -  tmp) + (dt / (2 * self.mass)) * (F + Fprev )) - self.vel
        self.vel += update_v
        
        new_acc = F / self.mass
        self.acc = F / self.mass

        return update_x, update_v, new_acc
    
    def update_with_object(self, obj, update_x, update_v, new_acc):
        '''
        Eventually this update will incorporate the current rotation parameters of the object
        The update values are the updates applied to the object's center of mass.  With rotation,
        this update will be more complex than a simple addition
        '''
        self.pos += update_x
        self.vel += update_v
        self.acc = np.array(new_acc, copy=True)


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
        self.pos = np.array(pos, dtype=float)
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
            return [],None,None
        
        # Save axis and overlap magnitude if this is the smallest gap
        current_overlap = poly1_bounds[1] - poly2_bounds[0]
        if current_overlap < overlap:
            fix_axis = normal_axis
            #fix_flag = flag
            fix_mag_flag = mag_flag
            overlap = current_overlap

    #bounceback_displacement = fix_axis * overlap
    #return fix_axis, overlap, fix_flag, fix_mag_flag
    return fix_axis, overlap, fix_mag_flag
    #return bounceback_displacement

