import matplotlib.pyplot as plt

def plot_energy(fname):
    data = open(fname).read().splitlines()
    data = data[4:]
    totals = [float(x[0]) for x in [y.split(',') for y in data]]
    kin = [float(x[1]) for x in [y.split(',') for y in data]]
    pot = [float(x[2]) for x in [y.split(',') for y in data]]
    heat = [float(x[3]) for x in [y.split(',') for y in data]]

    plt.plot(totals)
    plt.plot(kin)
    plt.plot(pot)
    plt.plot(heat)
    plt.legend(['total','kinetic','potential','heat'])
    plt.show()

