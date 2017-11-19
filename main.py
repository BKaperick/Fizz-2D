import numpy as np

GRAVITY = np.array([0,-9.8])

class world:
    def __init__(self, objs = set(), global_force = GRAVITY, time_disc = 1):
        self.objs = set()
        self.time_disc = time_disc
        self.state = 0
        self.global_force = global_force

        for obj in objs:
            self.init_obj(obj)
    
    def init_obj(self, obj):
        self.objs.add(obj)
        obj.forces.append(self.global_force)

    def update(self):
        for obj in self.objs:
            obj.pre_update(self.time_disc)
#        update = self.check_collisions()
#        for obj1,obj2 in update:
#            self.collide_gracefully(obj1, obj2)
        for obj in self.objs:
            obj.finish_update()
        self.state += 1
            

class obj:
    def __init__(self, world, mass = 1, pos = np.array([0.0,0.0]), speed = np.array([0.0,0.0]), forces = []):

        self.mass = mass
        self.pos = pos
        self.new_pos = pos
        self.vel = speed
        self.forces = forces
        self.force_update = True
        self.acc = np.array([0.0,0.0])
        
        self.world = world
        world.init_obj(self)

    def pre_update(self, dt):
        self.move(dt)

    def finish_update(self):
        self.pos = self.new_pos
    
    def get_forces(self, world):
        return False

    def move(self,dt,):
        '''
        Velocity Verlet Algorithm to update x,v,a with global error of order 2.
        '''
        # Update position
        v_avg = self.vel + (.5 * self.acc * dt)
        self.new_pos += v_avg * dt
        
        # Update acceleration
        self.acc = sum(self.forces) / self.mass
        
        # Update velocity
        self.vel = v_avg + (.5 * self.acc * dt)


if __name__ == '__main__':
    plane = world()
    ball = obj(plane)

    for t in range(10):
        print(ball.pos)
        plane.update()
    print(ball.pos)
    
