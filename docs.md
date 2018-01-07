
# ./main_old.py
## read_input(fname)

# ./physics.py
## World
###     __init__(self, width, height, objs = set(), global_accel = GRAVITY,
###     init_obj(self, obj)
###     init_fixed_obj(self, obj)
###     update(self)
>         This function is called once each time step.  All position and velocity
>         updates, and all collisions are handled within.

###     check_collisions(self)
>         Iterates through each pair of objects in the world, and determines if 
>         any are currently overlapping.
>         Currently, only polygons and fixed polygons are supported.
>         *TODO:* Implement collision detection for circles.

###     __str__(self)
>         Output the current state of the world as a string to be read by png 
>         script.

## Obj
###     __init__(self, points = [], world = None, mass = 1,
###     pre_update(self, force, damping_force, dt)
>         All points in the object are updated one step in accordance with the 
>         calculations done on the center of mass point.
>         The object's official position, velocity and acceleration are not updated
>         until self.finish_update().

###     reverse_update(self)
>         Can be used to take one step backward in time and undo a pre_update, 
>         since each point stores one previous position, velocity and acceleration.  
>         Currently, this is used in collision resolution to test out various time steps before
>         choosing one sufficiently close to the time of collision.

###     finish_update(self)
>         The object's position, velocity and acceleration are updated to be 
>         consistent with its center of mass.

## Point
###     __init__(self,
###     move(self, forces, dt)
>         Velocity Verlet Algorithm to update x,v,a with global error of order 2.

###     linear_damping_move(self, forces, dt)
###     update_with_object(self, obj, update_x, update_v, new_acc)
>         Eventually this update will incorporate the current rotation parameters of the object
>         The update values are the updates applied to the object's center of mass.  With rotation,
>         this update will be more complex than a simple addition

###     __str__(self)
## Polygon(Obj)
###     __init__(self, world = None, mass = 1, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
###     moment_of_inertia(self, cached_I = None)
>         Computes the moment of inertia by splitting polygon into triangles 
>         cornered at self.com.

###     __str__(self)
>         This prints the polygon object as a string compatible with draw.c.
>         Note: we leave out mass here since that is irrelevant to graphical 
>         representation at a fixed point in time.

## Ball(Obj)
###     __init__(self, world = None, mass = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
###     __str__(self)
## FixedPolygon(Polygon)
###     __init__(self, world = None, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
## polypoly_collision(poly1, poly2)
## triangle_area(p0,p1,p2)
## triangle_moment_of_inertia(p0,p1,p2)

# ./main.py
## read_input(fname)
