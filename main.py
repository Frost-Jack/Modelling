import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter
import tkinter as tk  # Python 3.x
from time import *
from tqdm import tqdm

master = tk.Tk()

data_file = open('data.txt', 'r')

data_array = [float(i) for i in data_file.readline().split()]
m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0 = data_array
data_vars = [m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0]

days, hours, minutes, seconds, iters = [2000, 0, 0, 0, 1000]

m1_text = tk.StringVar()
m2_text = tk.StringVar()
m3_text = tk.StringVar()

x1_0_text = tk.StringVar()
y1_0_text = tk.StringVar()
x2_0_text = tk.StringVar()
y2_0_text = tk.StringVar()
x3_0_text = tk.StringVar()
y3_0_text = tk.StringVar()

vx1_0_text = tk.StringVar()
vy1_0_text = tk.StringVar()
vx2_0_text = tk.StringVar()
vy2_0_text = tk.StringVar()
vx3_0_text = tk.StringVar()
vy3_0_text = tk.StringVar()

text_vars = [m1_text, m2_text, m3_text, x1_0_text, y1_0_text, x2_0_text, y2_0_text, x3_0_text, y3_0_text,
             vx1_0_text, vy1_0_text, vx2_0_text, vy2_0_text, vx3_0_text, vy3_0_text]

print(m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0)


def save_data():
    global data_file, data_fields, message_label, flag
    data_file.close()
    data_file = open('data.txt', 'w')
    for data_field in data_fields:
        data_file.write(data_field.get().strip())
        data_file.write(' ')


def update_quit():
    global data_array, master, data_file
    global m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0, \
        days__, hours__, minutes__, seconds__, iters__, days, hours, minutes, seconds, iters

    days = int(days__.get().strip())
    hours = int(hours__.get().strip())
    minutes = int(minutes__.get().strip())
    seconds = int(seconds__.get().strip())
    iters = int(iters__.get().strip())

    data_file.close()
    data_file = open('data.txt', 'r')
    data_array = [float(i) for i in data_file.readline().split()]
    m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0 = data_array
    data_file.close()
    master.destroy()


for i in range(len(text_vars)):
    text_vars[i].set(data_vars[i])

tk.Label(master, text="Masses").grid(row=1, column=3)
tk.Label(master, text="Coordinates").grid(row=3, column=3)
tk.Label(master, text="Speeds").grid(row=6, column=3)

m__1 = tk.Entry(master, textvariable=m1_text)
m__2 = tk.Entry(master, textvariable=m2_text)
m__3 = tk.Entry(master, textvariable=m3_text)

x__1_0 = tk.Entry(master, textvariable=x1_0_text)
y__1_0 = tk.Entry(master, textvariable=y1_0_text)
x__2_0 = tk.Entry(master, textvariable=x2_0_text)
y__2_0 = tk.Entry(master, textvariable=y2_0_text)
x__3_0 = tk.Entry(master, textvariable=x3_0_text)
y__3_0 = tk.Entry(master, textvariable=y3_0_text)

vx__1_0 = tk.Entry(master, textvariable=vx1_0_text)
vy__1_0 = tk.Entry(master, textvariable=vy1_0_text)
vx__2_0 = tk.Entry(master, textvariable=vx2_0_text)
vy__2_0 = tk.Entry(master, textvariable=vy2_0_text)
vx__3_0 = tk.Entry(master, textvariable=vx3_0_text)
vy__3_0 = tk.Entry(master, textvariable=vy3_0_text)

days__ = tk.Entry(master)
hours__ = tk.Entry(master)
minutes__ = tk.Entry(master)
seconds__ = tk.Entry(master)
iters__ = tk.Entry(master)

m__1.grid(row=2, column=1)
tk.Label(master, text="m1").grid(row=2, column=0)
m__2.grid(row=2, column=3)
tk.Label(master, text="m2").grid(row=2, column=2)
m__3.grid(row=2, column=5)
tk.Label(master, text="m3").grid(row=2, column=4)

x__1_0.grid(row=4, column=1)
tk.Label(master, text="x1").grid(row=4, column=0)
y__1_0.grid(row=5, column=1)
tk.Label(master, text="y1").grid(row=5, column=0)
x__2_0.grid(row=4, column=3)
tk.Label(master, text="x2").grid(row=4, column=2)
y__2_0.grid(row=5, column=3)
tk.Label(master, text="y2").grid(row=5, column=2)
x__3_0.grid(row=4, column=5)
tk.Label(master, text="x3").grid(row=4, column=4)
y__3_0.grid(row=5, column=5)
tk.Label(master, text="y3").grid(row=5, column=4)

vx__1_0.grid(row=7, column=1)
tk.Label(master, text="vx1").grid(row=7, column=0)
vy__1_0.grid(row=8, column=1)
tk.Label(master, text="vy1").grid(row=8, column=0)
vx__2_0.grid(row=7, column=3)
tk.Label(master, text="vx2").grid(row=7, column=2)
vy__2_0.grid(row=8, column=3)
tk.Label(master, text="vy2").grid(row=8, column=2)
vx__3_0.grid(row=7, column=5)
tk.Label(master, text="vx3").grid(row=7, column=4)
vy__3_0.grid(row=8, column=5)
tk.Label(master, text="vy2").grid(row=8, column=4)

tk.Label(master, text="Days, hours, minutes, seconds").grid(row=9, column=2, columnspan=3)
days__.grid(row=10, column=1)
hours__.grid(row=10, column=2)
minutes__.grid(row=10, column=3)
seconds__.grid(row=10, column=4)

tk.Label(master, text="Iters").grid(row=11, column=3)
iters__.grid(row=12, column=3)

data_fields = [m__1, m__2, m__3, x__1_0, y__1_0, x__2_0, y__2_0, x__3_0, y__3_0,
               vx__1_0, vy__1_0, vx__2_0, vy__2_0, vx__3_0, vy__3_0, ]

tk.Button(master, text='Save', command=save_data).grid(row=14, column=1)
tk.Button(master, text='Quit', command=update_quit).grid(row=14, column=5)

master.mainloop()

print(m1, m2, m3, x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0)

k = 1 / np.sqrt(6.67e-11  # Gravitational constant
                * 1.99e30  # Mass of Sun
                / 1.5e11 ** 3  # Astronomical unit
                )

time_to_simulate = seconds + minutes * 60 + hours * 60 * 60 + days * 24 * 60 * 60
dt = time_to_simulate / iters
years_per_step = dt / (60 * 60 * 24 * 365.25)
timing = dt / k * iters
print(dt)

sleep(2)


def dSdt(t, S):
    x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3 = S
    r12 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    r13 = np.sqrt((x3 - x1) ** 2 + (y3 - y1) ** 2)
    r23 = np.sqrt((x2 - x3) ** 2 + (y2 - y3) ** 2)
    return [vx1,
            vy1,
            vx2,
            vy2,
            vx3,
            vy3,
            m2 / r12 ** 3 * (x2 - x1) + m3 / r13 ** 3 * (x3 - x1),  # mass 1
            m2 / r12 ** 3 * (y2 - y1) + m3 / r13 ** 3 * (y3 - y1),
            m1 / r12 ** 3 * (x1 - x2) + m3 / r23 ** 3 * (x3 - x2),  # mass 2
            m1 / r12 ** 3 * (y1 - y2) + m3 / r23 ** 3 * (y3 - y2),
            m1 / r13 ** 3 * (x1 - x3) + m2 / r23 ** 3 * (x2 - x3),  # mass 3
            m1 / r13 ** 3 * (y1 - y3) + m2 / r23 ** 3 * (y2 - y3)
            ]


t = np.linspace(0, timing, iters)

sol = solve_ivp(dSdt, (0, timing), y0=[x1_0, y1_0, x2_0, y2_0, x3_0, y3_0, vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0],
                method='DOP853', t_eval=t, rtol=1e-10, atol=1e-13)

t = sol.t
x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[2]
y2 = sol.y[3]
x3 = sol.y[4]
y3 = sol.y[5]
Ep = np.ndarray(iters)
Ek = np.ndarray(iters)
#print(k)

for i in range(iters):
    Ep[i] = (6.67e-11 * m1 * m2 / np.sqrt((x2[i] - x1[i]) ** 2 + (y2[i] - y1[i]) ** 2) / 1.5e11 ** 2 * 1.99e30 ** 2 +
             6.67e-11 * m1 * m3 / np.sqrt((x3[i] - x1[i]) ** 2 + (y3[i] - y1[i]) ** 2) / 1.5e11 ** 2 * 1.99e30 ** 2 +
             6.67e-11 * m2 * m3 / np.sqrt((x3[i] - x2[i]) ** 2 + (y3[i] - y2[i]) ** 2) / 1.5e11 ** 2 * 1.99e30 ** 2)
    Ek[i] = (m1 / 2 * (sol.y[6][i] ** 2 + sol.y[7][i] ** 2) + m2 / 2 * (
            sol.y[8][i] ** 2 + sol.y[9][i] ** 2) + m3 / 2 * (
                     sol.y[10][i] ** 2 + sol.y[11][i] ** 2)) * 1.99e30 * k ** 2 * 1.5e11 ** 2
    #print(Ek[i] - Ep[i])

# tt = 1 / np.sqrt(6.67e-11  # Gravitational constant
#                  * 1.99e30  # Mass of Sun
#                  / 1.5e11 ** 3  # Astronomical unit
#                  )
# print(tt)
# tt = tt / (60 * 60 * 24 * 365.25) * np.diff(t)[0]  # per time step (in years)
# print(np.diff(t)[0])


times = []
Ep_array = []
Ek_array = []
max_times = time_to_simulate / (60 * 60 * 24 * 365.25)
max_E = max(Ep.max(), Ek.max())


def animate(i):
    # ln1.set_data([x1[i], x2[i], x3[i]], [y1[i], y2[i], y3[i]])
    ln1.set_data([x1[i]], [y1[i]])
    ln2.set_data([x2[i]], [y2[i]])
    ln3.set_data([x3[i]], [y3[i]])
    times.append(i * years_per_step)
    Ep_array.append(Ep[i])
    Ek_array.append(Ek[i])
    ax2.clear()
    ax2.set_xlim(0, max_times)
    ax2.set_ylim(0, max_E)
    ax2.plot(times, Ep_array)
    ax2.plot(times, Ek_array)
    text1.set_text('Time = {:.3f} Years'.format(i * years_per_step))


max_x = max(x1.max(), x2.max(), x3.max(), abs(x1.min()), abs(x2.min()), abs(x3.min()))
max_y = max(y1.max(), y2.max(), y3.max(), abs(y1.min()), abs(y2.min()), abs(y3.min()))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
ax1.grid()
ln1, = ax1.plot([], [], 'ro', lw=3, markersize=3)
ln2, = ax1.plot([], [], 'bo', lw=3, markersize=3)
ln3, = ax1.plot([], [], 'yo', lw=3, markersize=3)
ln4, = ax2.plot([], [], 'ro', lw=3, markersize=2)
text1 = ax1.text(0, max_y * 1.1, '', fontsize=15, backgroundcolor='white', ha='center')
plt.subplots_adjust(top=0.9)
ax1.set_xlim(-max_x, max_x)
ax1.set_ylim(-max_y, max_y)
ax2.set_xlim(-0.05, max_times)
ax2.set_ylim(-0.05, max_E)
ani = animation.FuncAnimation(fig, animate, frames=tqdm(range(iters), colour="green"), interval=50)
ani.save('plan.gif', writer='pillow', fps=30)
6