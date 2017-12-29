
./physics.py--------------------------------------------------------------------
World
    __init__(self, width, height, objs = set(), global_accel = GRAVITY, time_disc = 1, gamma = [])
    init_obj(self, obj)
    init_fixed_obj(self, obj)
    update(self)
        This function is called once each time step.  All position and velocity
        updates, and all collisions are handled within.

    check_collisions(self)
    __str__(self)
Obj
    __init__(self, points = [], world = None, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
    pre_update(self, force, damping_force, dt)
    reverse_update(self)
    finish_update(self)
Point
    __init__(self, world, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0]))
    move(self, forces, dt)
        Velocity Verlet Algorithm to update x,v,a with global error of order 2.

    linear_damping_move(self, forces, dt)
    update_with_object(self, obj, update_x, update_v, new_acc)
        Eventually this update will incorporate the current rotation parameters of the object
        The update values are the updates applied to the object's center of mass.  With rotation,
        this update will be more complex than a simple addition

Polygon(Obj)
    __init__(self, world = None, mass = 1, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
    __str__(self)
Ball(Obj)
    __init__(self, world = None, mass = 1, pos = np.array([0.0,0.0]), radius = 1, speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
    __str__(self)
FixedPolygon(Polygon)
    __init__(self, world = None, points = [], speed = np.array([0.0,0.0]), rotation_angle = 0.0, rotation_speed = 0.0)
polypoly_collision(poly1, poly2)

./main.py-----------------------------------------------------------------------
read_input(fname)
