###############################################################################
#                                                                             #
#       PHYSICS.PY-                                                           #
#       THIS FILE CONTAINS THE CLASSES FOR ALL OBJECTS IN THE WORLD           #
#       AND THEIR PHYSICAL PROPERTIES                                         #
#                                                                             #
###############################################################################

import numpy as np
from numpy.linalg import norm 
from   math import acos,pi
import config

# World constants
GRAVITY = np.array([0.0,9.8])
TIME_DISC = .05

DENSITY = .1
EPSILON = 1e-12
# 1 == perfectly elastic 
# 0 == perfectly inelastic
ELASTICITY = 1.0
TIME_TOL = 1e-10
COLLISION_TOL = .1
VIBRATE_TOL = 1e-5

ANGLE_TOL = 1e-4
CONTACT_TOL = 1e-2

class World:
    
    def __init__(self, width, height, objs = set(), global_accel = GRAVITY, 
                 time_disc = TIME_DISC, gamma = []):
        self.width = width
        self.height = height
        self.objs = []
        self.fixed_objs = []
        self.time_disc = time_disc
        self.state = 0
        if energylog:
            self.log = open(config.SIMULATION_DIR + 'energy.txt', 'a+')
        
        self.global_accel = global_accel
        self.global_damping_force = np.array([])
        if len(gamma) > 0:
            #self.global_damping_force = lambda v : gamma * v^2
            self.global_damping_force = gamma

        for obj in objs:
            self.init_obj(obj)

        # Heat contained in all objects
        self.heat = 0

    
    def init_obj(self, obj):
        if self.global_damping_force:
            obj.damping_force = True
        # Only force initially acting is gravity
        obj.forces.append(obj.mass*self.global_accel)
        self.objs.append(obj)

    def init_fixed_obj(self, obj):
        self.fixed_objs.append(obj)
    
    def energy(self):
        # total, kinetic, potential, heat
        E = np.array([0.0,0.0,0.0,0.0])
        for obj in self.objs:
            E[1] += obj.k_energy()
            E[2] += obj.p_energy()
        E[3] = self.heat    
        E[0] = sum(E[1:4])
        return E
    
    def update(self):
        '''
        This function is called once each time step.  All position and velocity
        updates, and all collisions are handled within.
        '''
        if energylog:
            self.log.write(','.join([str(x) for x in self.energy()]) + '\n')
        # Apply linear damping force if exists
        if len(self.global_damping_force) == 0:
            gdf = []
        else:
            gdf = [self.global_damping_force]
    
        dt = self.time_disc
        max_collision_overlap = COLLISION_TOL + 1
        first_iter = True

        # Check there is at least one 'true' collision, AND time is stable
        if verbosity>2:
            print('in:',self.objs[0].com.pos[1])
        if verbosity:
            hashes = [hash(obj) for obj in self.objs]
        while max_collision_overlap > COLLISION_TOL and dt > TIME_TOL:       
            if first_iter:
                first_iter = False
            else:
                dt /= 2
                for i,obj in enumerate(self.objs):
                    obj.reverse_update()
                    if verbosity:
                        if hash(obj) != hashes[i]:
                            print('WRONG, dt = ',dt)
                            print(obj)
            
            # Do first pass of position, velocity and acceleration updates for 
            # each object in the world
            for obj in self.objs:
                #obj.forces = []
                obj.pre_update(gdf, dt)

            # check_collisions() compiles a list of collisions with information:
            #   (1,2) the two objects contained in the collision
            #   (3) the unit vector normal to the surface of collision
            #   (4) magnitude of normal vector
            #   (5) a flag indicating a sign change for the normal vector
            collisions = self.check_collisions()
            magnitudes = [abs(magnitude) for _,__,___,magnitude,____ in collisions]
            if collisions:
                max_collision_overlap = max(magnitudes)
            else:
                max_collision_overlap = 0
            if max_collision_overlap > 0:# and dt < TIME_TOL:
                if verbosity == 1:
                    print(dt,max_collision_overlap, self.objs[0].com.pos[1])
        if dt < TIME_TOL:
            print("Warning: time became too small while still in collision.  ")
        
        while collisions:
           
            for obj1, obj2, unit_normal_vec, normal_vec_mag, flag in collisions:
                
                # Makes things easier if the first object is never fixed
                if obj1.is_fixed:
                    obj1,obj2 = obj2,obj1
                
                [side_length, left, right, normal] = side_contact(obj1,obj2)
                if side_length:
                    # Check if torque is induced 
                    # (TODO: this is not gravity-direction neutral)
                    if obj1.com.pos[0] < left[0] or obj1.com.pos[0] > right[0]:
                        pass # induce torque
                
                # Reassemble normal vector with some added spacing so two 
                # objects are definitely non-overlapping after collision 
                # resolution
                normal_vec = unit_normal_vec * (normal_vec_mag + EPSILON)
                if verbosity>1:
                    print("normal vec: ", normal_vec)
                
                # If flag, the normal vector needs to be flipped
                if flag == 1:
                    normal_vec *= -1
                    unit_normal_vec *= -1


                # obj1 is NEVER fixed, so this branch is for a collision with a 
                # free object and a fixed object
                #
                # 1. The free object updates its position by the normal vector
                # 2. The free object updates its velocity by twice the (negated)
                #    component of its velocity normal to the surface
                #
                #print(obj1.is_fixed, obj2.is_fixed)
                if obj2.is_fixed:
                    
                    #Jpart = (obj1.mass*obj2.mass/(obj1.mass + obj2.mass))*(1+CR)
                    #obj1.com.vel += (Jpart/obj1.mass + Jpart/obj2.mass)*(v2n - v1n)*unit_normal_vec
                    for i in range(len(obj1.points)):
                        obj1.points[i].pos += normal_vec

                        ## This would work for inelastic collisions with others
                        #obj1.points[i].vel += (Jpart/obj1.mass + Jpart/obj2.mass)*(unit_normal_vec)*(v2n - np.dot(obj1.points[i].vel,unit_normal_vec))
                        
                        vdotN = np.dot(obj1.points[i].vel, unit_normal_vec)
                        obj1.points[i].vel -= (1+ELASTICITY)*vdotN*unit_normal_vec

                    # TODO: why is this here?
                    if norm(normal_vec) > VIBRATE_TOL:
                        obj1.com.pos += normal_vec
                    
                    vdotN = np.dot(obj1.com.vel, unit_normal_vec)
                    #if norm(vdotN) > VIBRATE_TOL:
                    ##obj1.com.vel -= (1+ELASTICITY)*vdotN * unit_normal_vec
                    #J = obj1.mass*
                    K0 = .5*obj1.mass*norm(obj1.com.vel)**2
                    obj1.com.vel -= (1+ELASTICITY)*vdotN * unit_normal_vec
                    
                    # Update heat energy, dispersed equally to each object
                    H = .5*obj1.mass*(1-ELASTICITY**2)*(vdotN**2)
                    K1 = .5*obj1.mass*norm(obj1.com.vel)**2
                    obj1.heat += H/2
                    obj2.heat += H/2
                    self.heat += H

                    obj1.finish_update()
                    obj2.finish_update()

                # Only other possible case is that both objects are free
                else:
                    v1 = obj1.com.vel
                    mtotal = obj1.mass + obj2.mass
                    mprop1 = obj1.mass / mtotal 
                    mprop2 = 1-mprop1 
                    obj1.com.pos += mprop1 * (normal_vec)
                    vn1 = np.dot(obj1.com.vel,unit_normal_vec)
                    vn2 = np.dot(obj2.com.vel,unit_normal_vec)
                    v1_update = (mprop1 - mprop2)*vn1 + 2*mprop2*vn2
                    if verbosity>1:
                        v0a = np.copy(obj1.com.vel)
                        v0b = np.copy(obj2.com.vel)
                        print("inits", v0a,v0b)
                        K0a = .5*obj1.mass*norm(obj1.com.vel)**2
                        K0b = .5*obj2.mass*norm(obj2.com.vel)**2
                        print('K0a= ', K0a)
                        print('K0b= ', K0b)
                        #print('mbefore:',obj1.momentum(normal_vec) , obj2.momentum(normal_vec))
                    obj1.com.vel += (v1_update - np.dot(obj1.com.vel,unit_normal_vec))*unit_normal_vec
                    for point in obj1.points:
                        point.pos += mprop1 * (normal_vec)
                        point.vel += (v1_update - np.dot(point.vel,unit_normal_vec))*unit_normal_vec
                        
                    obj1.finish_update()
                    
                    obj2.com.pos -= mprop2 * normal_vec
                    v2_update = 2*mprop1*vn1 + (mprop2 - mprop1)*vn2
                    obj2.com.vel += (v2_update - np.dot(obj2.com.vel,unit_normal_vec))*unit_normal_vec
                    for point in obj2.points:
                        point.pos -= mprop2 * (normal_vec)
                        point.vel += (v2_update - np.dot(point.vel,unit_normal_vec))*unit_normal_vec
                    obj2.finish_update()
                    if verbosity>1:
                        print("inits", v0a,v0b)
                        print("a_update=",vn2, v1_update)
                        print("b_update=",vn1, v2_update)
                        print("va1 =", obj1.com.vel, v0a+vn1)
                        print("vb1 =", obj2.com.vel, v0b+vn2)
                        print('K1a=',K0a + .5*obj1.mass*(vn1**2 - vn2**2), .5*obj1.mass*norm(obj1.com.vel)**2)
                        print('K1b=',K0b + .5*obj2.mass*(vn2**2 - vn1**2), .5*obj2.mass*norm(obj2.com.vel)**2)
                        #print('mafter:',obj1.momentum(normal_vec) , obj2.momentum(normal_vec))

            collisions = self.check_collisions()

        # Do second (final) pass for each object in the world
        for obj in self.objs:
            # Obj.com x,v,a attrs are given to Obj itself
            obj.finish_update()
        
        # Update the remainder of the time step
        # ensures each frame transition is consistent time width
        for obj in self.objs:
            obj.pre_update(gdf, self.time_disc - dt)
            obj.finish_update()
        
        # Advance time
        self.state += self.time_disc


    def check_collisions(self):
        '''
        Iterates through each pair of objects in the world, and determines if 
        any are currently overlapping.
        Currently, only polygons and fixed polygons are supported.
        *TODO:* Implement collision detection for circles.
        '''
        collisions = []
        objects = [o for o in self.objs + self.fixed_objs 
                if o.name == "polygon" or o.name == "fixedpolygon"]
        for i,obj in enumerate(objects):
            for other_obj in objects[i+1:]:
                
                # Don't even check collisions if the two objects are not close 
                # to each other
                com_dist = norm(obj.pos - other_obj.pos)
                if com_dist > max(obj.radius, other_obj.radius):
                    continue

                if verbosity >= 2:
                    print([p.pos for p in obj.points])
                    print([p.pos for p in other_obj.points])
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
    def __init__(self, points = [], world = None, mass = 1, density = 1, 
                pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0]), 
                rotation_speed = .00):
        
        self.cached_I = None
        self.mass = mass
        self.density = density
        self.points = points
        
        com_pos = sum([pt.pos for pt in points]) / len(points)
        self.com = Point(
                world, 
                pos = com_pos, 
                mass = mass,
                speed = speed,
                rotvel = rotation_speed
                )
        self.radius = max([norm(pt.pos - com_pos) for pt in points])
        
        # Attributes of motion
        self.pos = pos
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rotpos = 0.0
        self.rotvel = rotation_speed
        self.rotacc = 0.0

        self.update_x = 0.0
        self.update_v = 0.0
        self.new_acc = 0.0
        
        self.update_rotx = 0.0
        self.update_rotv = 0.0
        self.new_rotacc = 0.0
        
        self.world = world
        self.forces = []
        self.torques = []
        if world:
            world.init_obj(self)

        self.is_fixed = False
        
        self.heat = 0.0

    def rotate(self, angle,dontupdatecom=True):
        '''
        Rotate Object angle radians
        '''
        for point in self.points:
            radius = (point.pos - self.com.pos)
            radius = radius*norm(point.radius)/norm(radius)
            th = point.rotpos + angle
            point.pos[0] = norm(radius)*np.cos(point.rotpos+angle) + self.com.pos[0]
            point.pos[1] = norm(radius)*np.sin(point.rotpos+angle) + self.com.pos[1]
            point.rotpos += angle
        if not dontupdatecom:
            self.rotpos += angle
            self.com.rotpos += angle

    def momentum(self, direction):
        return self.mass * np.dot(self.com.vel, direction)/norm(direction)
    
    def k_energy(self):
        return .5*self.mass * norm(self.vel)**2
    
    def p_energy(self):
        # TODO: Fix to always be parallel to gravity in arbitrary direction 
        p = GRAVITY[1] * self.mass * (self.world.height-self.pos[1])
        return p

    def side_pos_iter(self):
        '''
        Simple generator to yield a common need:
        endpoint positions of each side, and a vector centered at zero
        '''
        for p1,p2 in zip(self.points, self.points[1:]+[self.points[0]]):
            yield p1.pos,p2.pos,p2.pos-p1.pos

    def pre_update(self, damping_force, dt):
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
        
        # Retain copy of previous step
        self.com.oldpos = np.copy(self.com.pos)
        self.com.oldvel = np.copy(self.com.vel)
        self.com.oldacc = np.copy(self.com.acc)
        
        # Retain copy of previous step
        self.com.oldrotpos = np.copy(self.com.rotpos)
        self.com.oldrotvel = np.copy(self.com.rotvel)
        self.com.oldrotacc = np.copy(self.com.rotacc)
        #print(self.com.rotvel)
        # Use time-integrator of choice to find new x,v,a
        self.update_x, self.update_v, self.new_acc = self.com.move(self.forces, dt)
        
        # Use time-integrator of choice to find new x,v,a
        self.update_rotx, self.update_rotv, self.new_rotacc = self.com.move_angular(self.torques, dt)
        self.rotate(self.update_rotx)
        for point in self.points:

            # Retain copy of previous step
            point.oldpos = np.copy(point.pos)
            point.oldvel = np.copy(point.vel)
            point.oldacc = np.copy(point.acc)
            
            # Retain copy of previous step
            point.oldrotpos = np.copy(point.rotpos)
            point.oldrotvel = np.copy(point.rotvel)
            point.oldrotacc = np.copy(point.rotacc)
            
            # update point x,v,a values
            point.update_with_object(self, self.update_x, self.update_v, self.new_acc, self.update_rotx, self.update_rotv, self.new_rotacc)
            #print('\t',point.rotvel,self.update_rotv)


    def reverse_update(self):
        '''
        Can be used to take one step backward in time and undo a pre_update, 
        since each point stores one previous position, velocity and acceleration.  
        Currently, this is used in collision resolution to test out various time steps before
        choosing one sufficiently close to the time of collision.
        '''
        self.com.pos = self.com.oldpos
        self.com.vel = self.com.oldvel
        self.com.acc = self.com.oldacc
        
        self.com.rotpos = self.com.oldrotpos
        self.com.rotvel = self.com.oldrotvel
        self.com.rotacc = self.com.oldrotacc

        for point in self.points:
            point.pos = point.oldpos
            point.vel = point.oldvel
            point.acc = point.oldacc 
            point.rotpos = point.oldrotpos
            point.rotvel = point.oldrotvel
            point.rotacc = point.oldrotacc 

    def finish_update(self):
        '''
        The object's position, velocity and acceleration are updated to be 
        consistent with its center of mass.
        '''
        self.pos = np.array(self.com.pos, copy=True)
        self.vel = np.array(self.com.vel, copy=True)
        self.acc = np.array(self.com.acc, copy=True)
        
        self.rotpos = np.array(self.com.rotpos, copy=True)
        self.rotvel = np.array(self.com.rotvel, copy=True)
        self.rotacc = np.array(self.com.rotacc, copy=True)
    
    def __str__(self):
        return str(self)
        
        

class Point:
    def __init__(self, 
                world, 
                mass = 1, 
                pos = np.array([0.0,0.0]),
                speed = np.array([0.0,0.0]),
#                forces = np.array([0.0,0.0])
                I = 0.0,
                rotvel = 0.0
                ):
        
        self.name = "point"

        self.oldpos = np.array(pos, dtype=float)
        self.pos = np.array(pos, dtype=float)
        
        self.rotpos = 0.0
        self.rotvel = rotvel
        self.rotacc = 0.0

        self.oldvel = speed
        self.vel = speed
        self.oldacc = np.array([0.0,0.0])
        self.acc = np.array([0.0,0.0])
#        self.forces = np.array(mass*world.global_accel)

        self.world = world
        self.mass = mass
#        self.moment_of_inertia = I


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
        new_acc = (sum(forces) / self.mass)
        self.acc = new_acc
        
        # Update velocity
        update_v += self.acc * halfdt
        self.vel = v_avg + (self.acc * halfdt)
        if verbosity > 1:
            print('update', update_x)
        if verbosity==1:
            print('dt:',dt)
            print('updatex:', update_x)
            print('updatev:', update_v)
            print('new acc:', new_acc)
        return update_x, update_v, new_acc
    
    def move_angular(self, torques, dt):
        '''
        Velocity Verlet Algorithm to update x,v,a with global error of order 2.
        '''

        halfdt = .5*dt
        update_v = self.rotacc * halfdt
        v_avg = self.rotvel + update_v
        
        # Update position
        update_x = v_avg * dt
        self.rotpos = update_x
        
        # Update acceleration
        new_acc = (sum(torques) / self.moment_of_inertia)
        self.rotacc = new_acc
        
        # Update velocity
        update_v += self.rotacc * halfdt
        self.rotvel = v_avg + (self.rotacc * halfdt)
        if verbosity > 1:
            print('update', update_x)
        if verbosity>1:
            print('updatex:', update_x)
            print('updatev:', update_v)
            print('new acc:', new_acc)
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
    
    def update_with_object(self, obj, update_x, update_v, new_acc,
                                      update_rotx,update_rotv,new_rotacc):
        '''
        Eventually this update will incorporate the current rotation parameters of the object
        The update values are the updates applied to the object's center of mass.  With rotation,
        this update will be more complex than a simple addition
        '''
        if verbosity > 1:
            print(update_x, update_v, new_acc)
        
        radius = norm(self.pos - obj.com.pos)
        self.pos += update_x
        obj.rotate(update_rotx)
        #self.pos[0] += radius*(np.sin(update_rotx+self.rotpos) - np.sin(self.rotpos))
        #self.pos[1] += radius*(np.cos(update_rotx+self.rotpos) - np.cos(self.rotpos))
        self.vel += update_v + radius*update_rotv
        self.acc = np.array(new_acc, copy=True) + radius*new_rotacc
        
        self.rotpos += update_rotx
        self.rotvel += update_rotv
        self.rotacc = np.array(new_rotacc, copy=True)

        com = obj.com

    def __str__(self):
        int_pos_x = int(self.pos[0])
        int_pos_y = int(self.pos[1])
        bnded_pos_x = max(min(int_pos_x, self.world.width-1), 0)
        bnded_pos_y = max(min(int_pos_y, self.world.height-1), 0)
        #return "point," + str(int(self.pos[0])) + "," + str(int(self.pos[1]))
        return "point," + str(bnded_pos_x) + "," + str(bnded_pos_y)


class Polygon(Obj):
    def __init__(self, world = None, density = 1, mass = 1, points = [], speed = np.array([0.0,0.0]), rotation_speed = .0):
        super().__init__(world = world, mass = mass, points = points, speed = speed, rotation_speed = rotation_speed)
        self.name = "polygon"
        self.num_edges = len(points)
        
        self.pos = self.com.pos
        self.vel = self.com.vel
        self.rotvel = rotation_speed

        if density == np.inf:
            self.com.moment_of_inertia = np.inf
            self.mass = np.inf
        else:
            self.com.rotvel = rotation_speed
            self.com.moment_of_inertia = self.moment_of_inertia()
            self.mass = self.area*density
            if verbosity>1:
                print('mass:',self.mass)
                print('height:',self.pos[1])
                print('potent:',self.p_energy())
        
        # second pass to fix masses  for unfixed polygons
        for point in self.points:
            point.radius = point.pos - self.com.pos
            point.rotpos = np.arccos(point.radius[0]/norm(point.radius))
            if point.radius[1] > 0:
                point.rotpos = 2*np.pi - point.rotpos
            point.rotvel = rotation_speed    
            point.mass = self.mass
            # TODO: this is a kind of hacky way to differentiate 
            # fixed and unfixed polygons
            if density != np.inf:
                point.rotvel = 0.0
                point.moment_of_inertia = self.moment_of_inertia()
            else:
                point.moment_of_inertia = np.inf

#    def area(self):
#        if self.area:
#            return self.area
#        if len(self.points) == 3:
#            p0,p1,p2 = [p.pos for p in self.points]
#            self.area = triangle_area(p0,p1,p2)
#        else:
#            self.area = 0.0
#            p0 = self.com.pos
#            for point1,point2 in zip(self.points[:-1], self.points[1:]):
#                
#                p1 = point1.pos
#                p2 = point2.pos
#                area = triangle_area(p0,p1,p2)
#                self.area += area
#        return self.area


    def moment_of_inertia(self):#, cached_I = None):
        '''
        Computes the moment of inertia by splitting polygon into triangles 
        cornered at self.com.
        '''
        if self.cached_I:
            return self.cached_I
        if len(self.points) == 3:
            p0,p1,p2 = [p.pos for p in self.points]
            #area = triangle_area(p0,p1,p2)
            #com = (p0 + p1 + p2) / 3
            I = triangle_moment_of_inertia(p0,p1,p2)
            self.cached_I = I
            self.area = triangle_area(p0,p1,p2)
        else:
            I = 0.0
            weighted_masses = 0.0
            self.area = 0.0
            p0 = self.com.pos
            #for point1,point2 in zip(self.points,deque(self.points).rotate(-1)):
            for point1,point2 in zip(self.points, self.points[1:]+[self.points[0]]):
            #for point1,point2 in zip(self.points[:-1], self.points[1:]):
                
                p1 = point1.pos
                p2 = point2.pos
                I_tri = triangle_moment_of_inertia(p0,p1,p2)

                # Weight the mass of each triangle by its distance from the COM of the polygon
                com = (p0 + p1 + p2) / 3
                dsquared = norm(com - p0)

                
                # Compute area of polygon "for free"
                area = triangle_area(p0,p1,p2)
                self.area += area

                # Once weighted_masses is scaled by M/A we add it into I.
                # However, we don't know A until the end of the loop
                I += I_tri
                weighted_masses += area * dsquared
            
            # Now that self.area is fully computed, we add components together
            weighted_masses *= (self.mass / self.area) 
            I += weighted_masses

            # Save value for next time
            self.cached_I = I
        return I
                
    def __str__(self):
        '''
        This prints the polygon object as a string compatible with draw.c.

        Note: we leave out mass here since that is irrelevant to graphical 
        representation at a fixed point in time.
        '''
        
        #return "polygon\nmass," + str(self.mass) + "\nsides," + str(len(self.points))+ "\n" + "\n".join([str(pt) for pt in self.points]) + "\n"
        return "polygon\nsides," + str(len(self.points))+ "\n" + "\n".join([str(pt) for pt in self.points]) + "\n"


class Ball(Obj):
    def __init__(self, world = None, density = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_speed = 0.0):
        
        self.name = "circle"

        self.points = []
        self.radius = radius
        self.mass = (self.radius**2) * pi * density
        self.com = Point(world, pos = pos, mass = self.mass, speed = speed)
        
        self.mass = mass
        
        # Attributes of motion
        self.pos = np.array(pos, dtype=float)
        self.vel = speed
        self.acc = np.array([0.0,0.0])
        self.rotvel = rotation_speed
        
        self.world = world
        if world:
            world.init_obj(self)

    def __str__(self):
        return "circle\n" + str(int(self.pos[0])) + "," + str(int(self.pos[1])) + "," + str(int(self.radius)) + "\n"



class FixedPolygon(Polygon):
    def __init__(self, world = None, points = [], speed = np.array([0.0,0.0]), rotation_speed = 0.0):
        super().__init__(world = None, density = np.inf, mass = np.inf, points = points, speed = np.array([0.0,0.0]), rotation_speed = 0.0)
        self.name = "fixedpolygon"        
        self.world = world
        self.mass = np.inf
        if world:
            world.init_fixed_obj(self)
        self.is_fixed = True

def moment_of_inertia(self):
    return np.inf

def polypoly_collision(poly1, poly2):
    '''
    Determine whether two polygons are currently intersecting each other
    '''
    
    # array of edges in the form
    # poly1_edges = [(Point1, Point2), (Point2, Point3), ..., (Point[n], Point1)]
    poly1_edges = zip(poly1.points, poly1.points[1:] + [poly1.points[0]])
    
    # array of outward-facing normals for each of the previous edges
    poly1_normals = [np.array([p2.pos[1] - p1.pos[1], p1.pos[0] - p2.pos[0]]) for p1,p2 in poly1_edges]
    
    poly2_edges = zip(poly2.points, poly2.points[1:] + [poly2.points[0]])
    poly2_normals = [np.array([p2.pos[1] - p1.pos[1], p1.pos[0] - p2.pos[0]]) for p1,p2 in poly2_edges]
    
    if verbosity >= 2:
        print('IN COLLISION')
        if verbosity >= 3:
            print([p.pos for p in poly1.points[:-1]])
            print([p.pos for p in poly2.points[:-1]])
    
    overlap = np.inf
    for axis,flag in list(zip(poly1_normals,[1]*len(poly1_normals))) + list(zip(poly2_normals, [2]*len(poly2_normals))):
        
        # Normalize
        normal_axis = axis / norm(axis)
        
        # Determine projection bounds onto separating axis (perpindicular to separating line

        # initialize bounds with current edge
        poly1_bounds = [np.dot(normal_axis, poly1.points[0].pos)]*2
        poly2_bounds = [np.dot(normal_axis, poly2.points[0].pos)]*2
        
        # increase bounds as necessary if other edges increase projection line length
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
        if verbosity>1:
            print('co:',current_overlap,poly1_bounds,poly2_bounds)
        if current_overlap < overlap:
            fix_axis = normal_axis
            fix_mag_flag = mag_flag
            overlap = current_overlap
    
    '''
    Sanity Check: This whole fix_axis direction makes physical sense.  It 
    is resolving the collision in a direction which *must* be orthogonal to one
    of the two objects' sides, which should be the only way a collision could
    happen, aside from some point-point intersection which is a.s. not the case.  
    '''    
    return fix_axis, overlap, fix_mag_flag
    

def triangle_area(pA,pB,pC):
    if verbosity>1:
        print(pA,pB,pC)
    area = .5 * abs(pA[0]*pB[1] + pB[0]*pC[1] + pC[0]*pA[1] - pA[0]*pC[1] - pC[0]*pB[1] - pB[0]*pA[1])
    return area

def triangle_moment_of_inertia(p0,p1,p2):

    # If the line from p0,p2 is close to vertical, the slope calculation is 
    # unstable - We resolve this by simply cycling the point values
    if abs(p2[0] - p0[0]) < EPSILON:
        p1,p2,p0 = p0,p1,p2

    slope = (p2[1] - p0[1]) / (p2[0] - p0[0])
    b = norm(p2 - p0)

    # Create new point where argmin of distance from p1 to the line 
    # between p0 and p2 occurs
    xmin = (slope**2 * p0[0] + slope * p1[1] - slope * p0[1] + p1[0]) / (slope**2 + 1)
    pmin = np.array([xmin, slope*(xmin - p0[0]) + p0[1]])
    
    a = norm(p0 - pmin)
    h = norm(p1 - pmin)

    I = ((b**3)*h - (b**2)*h*a + b*h*(a**2) + b*(h**3)) / 36
    return I


def side_contact(poly1, poly2):
    '''
    Returns (length,normal) giving length of "face-to-face" contact, and if so
    left and right endpoints, and a unit vector facing poly1 normal to the 
    contact.
    '''
    
    # first check if the polygons can possibly be in contact
    # note: we don't need to factor in collision tol, since the radius is a 
    # maximum reached at a corner point, never at a side point, so this is
    # always a conservative check
    com_dist = norm(poly1.pos - poly2.pos)
    if com_dist > max(poly1.radius, poly2.radius):
        return 0, None, None, None
    best_side = None
    best_score = .01
    for i1,(p1p1,p1p2,s1) in enumerate(poly1.side_pos_iter()):
        for i2,(p2p1,p2p2,s2) in enumerate(poly2.side_pos_iter()):
            anglediff = abs(np.dot(s1,s2))/(norm(s1)*norm(s2))
            if  1-anglediff < ANGLE_TOL:
                
                # I believe this always has to be minimized with p1p2 
                # assuming coords entered consistently cw or ccw
                # Lines should be in reverse orientation to each other
                z = p1p1-p2p2
                score = norm(np.cross(z/norm(z),s2/norm(s2)))
                if score < best_score:
                    best_side = (i1,i2)
                    best_score = score
                    crosslen1 = norm(p1p1-p2p1)
                    crosslen2 = norm(p2p2-p1p2)
                    minlen = min(norm(s1),norm(s2))
                    contact_length = min(minlen, crosslen1, crosslen2)
                    
                    # Again, with CW assumption, this turns toward object 1
                    normal = np.array([p2p2[1] - p2p1[1], p2p1[0] - p2p2[0]])

                    a,b,c,d = p1p1,p1p2,p2p1,p2p2
                    if a[0] > b[0]:
                        a,b = b,a
                    if c[0] > d[0]:
                        c,d = d,c
                    if a[0] < c[0] and c[0] < b[0]:
                        left = c
                    else:
                        left = a
                    if a[0] < d[0] and d[0] < b[0]:
                        right = d
                    else:
                        right = b
    #assert(abs(contact_length - norm(right - left)) / abs(contact_length) < EPSILON)
    if not best_side:
        return 0,None, None, None
    unit_normal = -normal/norm(normal)
    return contact_length, left, right, unit_normal
