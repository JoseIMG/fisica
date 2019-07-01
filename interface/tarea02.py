import tkinter as tk
import matplotlib.animation as ani
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# Hélice Conica
def helice_conica():
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    def gen(n):
        t = -6 * np.pi

        a1 = 30

        a = 15
        for t in np.linspace(-5 * np.pi, 5 * np.pi, 200):
            yield np.array([a1 * np.sin(a) * t * np.cos(t), a1 * np.sin(a) * t * np.sin(t), a1 * np.cos(a) * t])
        theta = 0
        e = 2.718281
        a = 0.5
        an1 = 45
        an2 = 30

        while theta > -25 * np.pi:
            x = a * (e ** (np.sin(an1) * (1 / np.tan(an2) * theta))) * np.cos(theta)
            y = a * (e ** (np.sin(an1) * (1 / np.tan(an2) * theta))) * np.sin(theta)
            z = a * (e ** (np.sin(an1) * (1 / np.tan(an2) * theta))) * (1 / np.tan(an1))
            yield np.array([x, y, z])
            theta += -8 * np.pi / n

    def update(num, data, line):
        line.set_data(data[:2, :num])
        line.set_3d_properties(data[2, :num])

    N = 200
    N = 150
    data = np.array(list(gen(N))).T
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

    # Setting the axes properties
    ax.set_xlim3d([-400.0, 400.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([-400.0, 400.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-400.0, 400.0])
    ax.set_xlim3d([-10.0, 10.0])
    ax.set_xlabel('X')
    ax.set_ylim3d([-10.0, 10.0])
    ax.set_ylabel('Y')
    ax.set_zlim3d([0.0, 10.0])
    ax.set_zlabel('Z')

    ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=1000 / N, blit=False, repeat=False)

    plt.show()


# Hélice Circular
def helice_circular_1():
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    def gen(n):
        t = 0
        r = 5
        while t < 8 * np.pi:
            yield np.array([r * np.cos(t), r * np.sin(t), t])
            t += 5.5 * np.pi / n

    def update(num, data, line):
        line.set_data(data[:2, :num])
        line.set_3d_properties(data[2, :num])

    N = 100
    data = np.array(list(gen(N))).T
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], label='Curva Helice Circular')
    ax.legend()

    # Setting the axes properties

    ax.set_xlim3d([-8.0, 8.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([-8.0, 8.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, 16.0])
    ax.set_zlabel('Z')
    ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=16, repeat=False)
    plt.show()

    pass


# Corona Sinusoidal

# Curva de Viviani
def anima(num, ejes, linea):
    linea.set_data(ejes[:2, :num])
    linea.set_3d_properties(ejes[2, :num])


def curva_de_viviani():
    # Creamos una figura
    fig = plt.figure()
    # Establecemos que es una figura de tipo 3D
    grafico = fig.gca(projection='3d')
    # Titulo del grafico
    grafico.set_title("Curva Viviani")
    # limites del grafico
    grafico.set_xlim3d([-1.0, 2.0])  # el eje x entre 0 y 2 (se uso corchete por la cantidad de datos)
    grafico.set_xlabel('X')  # lo que se esta escrito en el eje x
    grafico.set_ylim3d([-1.0, 2.0])  # el eje y entre -1 y 2 (se uso corchete por la cantidad de datos)
    grafico.set_ylabel('Y')  # lo que se esta escrito en el eje y
    grafico.set_zlim3d([-2.0, 2.0])  # el eje z entre -2 y 2 (se uso corchete por la cantidad de datos)
    grafico.set_zlabel('Z')  # lo que se esta escrito en el eje z
    # El valor recuperado
    a = 1
    # Los valores de T
    t = np.linspace(-4, 4 * np.pi, 100)
    # Ejes
    x = a * (1 + np.cos(t))
    y = a * np.sin(t)
    z = 2 * a * np.sin(t / 2)
    eje = [x, y, z]
    # agregamos la lista en un array para dibujar la linea
    ejes = np.array(eje)

    linea, = grafico.plot(ejes[0, 0:1], ejes[1, 0:1], ejes[2, 0:1], 'm', label="Curva de Viviani", lw=5)
    # figura,funcion,fargs:(datos,linea dibujada),FPS,ms,blit
    animacion = animation.FuncAnimation(fig, anima, fargs=(ejes, linea), frames=80, interval=16, blit=False,
                                        repeat=False)
    # mostrar el grafico
    plt.show()


# Hipopoda
def Hipopoda():
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    def gen():
        for theta in np.linspace(0, 4 * np.pi, 99):
            yield np.array([20 + (50 - 20) * np.cos(theta), (50 - 20) * np.sin(theta),
                            2 * (5 * (50 - 20)) ** (1 / 2) * np.sin(theta / 2)])

    def update(num, data, line):
        line.set_data(data[:2, :num])
        line.set_3d_properties(data[2, :num])

    N = 100
    plt.rcParams['legend.fontsize'] = 12
    data = np.array(list(gen())).T
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], label='Hipopede de Eudoxo')

    ax.set_xlim3d([-50.0, 50.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([-50.0, 50.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-50.0, 50.0])
    ax.set_zlabel('Z')
    anim = ani.FuncAnimation(fig, update, N, fargs=(data, line), interval=1700 / N, repeat=False)
    ax.legend()
    plt.show()
    pass


# Espiral Cónica de Papus
def conico_papus():
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    def gen(n):
        t = -6 * np.pi

        a1 = 30

        a = 15
        for t in np.linspace(-5 * np.pi, 5 * np.pi, 200):
            yield np.array([a1 * np.sin(a) * t * np.cos(t), a1 * np.sin(a) * t * np.sin(t), a1 * np.cos(a) * t])

    def update(num, data, line):
        line.set_data(data[:2, :num])
        line.set_3d_properties(data[2, :num])

    N = 200
    data = np.array(list(gen(N))).T
    line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

    # Setting the axes properties
    ax.set_xlim3d([-400.0, 400.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([-400.0, 400.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-400.0, 400.0])
    ax.set_zlabel('Z')

    ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=1000 / N, blit=False, repeat=False)

    plt.show()


# Curva de Arquitas

# Horóptera

def sig(i, XYZ, horoptera):
    horoptera.set_data(XYZ[:2, :i])
    horoptera.set_3d_properties(XYZ[2, :i])
def horoptera():
    plt.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = Axes3D(fig)
    r = 3   #corresponde al radio
    c = 1   #es una constante
    t = np.linspace(-2.*np.pi, 2.0*np.pi, 100)
    #se crean los ejes
    z = t
    x = 2 * r / (1 + ((c ** 2) * (t ** 2)))
    y = 2 * r * c * t / (1 + ((c ** 2) * (t ** 2)))
    #se crea lista de ejes
    XYZ = [x, y, z]
    XYZ = np.array(list(XYZ))
    horoptera, = ax.plot(x[0:1], y[0:1], z[0:1])
    #limites de los ejes
    ax.set_xlim3d([-8.0, 8.0])
    ax.set_xlabel('X')
    ax.set_ylim3d([-8.0, 8.0])
    ax.set_ylabel('Y')
    ax.set_zlim3d([-8.0, 8.0])
    ax.set_zlabel('Z')
    ax.set_title("Horoptera")
    animacion = animation.FuncAnimation(fig, sig, fargs=(XYZ, horoptera), frames=100, blit=False, interval=16,
                                        repeat=False)
    plt.show()
# Curva Bicilindrica

if __name__ == '__main__':
    # Creación de Ventanas
    root = tk.Tk()
    root.wm_title("funciones")
    root.geometry("800x600")

    # Crear frame contenedor de los elementos
    frame = tk.Frame(root)
    frame.pack(padx=20, pady=20)
    horoptera_img = tk.PhotoImage(file="horoptera.png")
    curva_horoptera = tk.Button(master=frame, image=horoptera_img, command=horoptera, height="100", width="100")
    curva_horoptera.pack(side=tk.BOTTOM, padx=10, pady=10)
    # Añadir titulo
    label = tk.Label(frame, text="Curvas Paramétricas Famosas", height="2")
    label.pack(fill=tk.X, expand=1)
    hipopoda_im = tk.PhotoImage(file="Figure_1.gif")
    hipopoda_button = tk.Button(master=frame, text="Hipopoda", command=Hipopoda, image=hipopoda_im)
    hipopoda_button.pack(side=tk.BOTTOM, padx=10, pady=10)

    imagen_viviani = tk.PhotoImage(file="Viviani_curve.png")
    viviani = tk.Button(master=frame, image=imagen_viviani, width=100, height=100, command=curva_de_viviani)
    viviani.pack()
    Conica_im = tk.PhotoImage(file="helice_conica.gif")
    conica_button = tk.Button(master=frame, text="Hélice Cónica", command=helice_conica, image=Conica_im)
    conica_button.pack(side=tk.BOTTOM, padx=10, pady=10)

    helice_circular_im = tk.PhotoImage(file="helice_circular1.gif")
    helice_circular_button = tk.Button(master=frame, text="helice circular", command=helice_circular_1,
                                       image=helice_circular_im)
    helice_circular_button.pack(side=tk.BOTTOM, padx=10, pady=10)

    curva_de_ejemplo = tk.Button(master=frame, text="conica de papus", command=conico_papus)
    curva_de_ejemplo.pack(side=tk.BOTTOM, padx=10, pady=10)
    # Conica_im = tk.PhotoImage(file="conica_papu.gif")
    conica_button = tk.Button(master=frame, text="conica de papus", command=conico_papus, image=Conica_im)
    conica_button.pack(side=tk.BOTTOM, padx=10, pady=10)

    tk.mainloop()
