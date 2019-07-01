import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def Atras(ventana):
    ventana.quit()
    ventana.destroy()

def sig(i, XYZ, horoptera):
    horoptera.set_data(XYZ[:2, :i])
    horoptera.set_3d_properties(XYZ[2, :i])
def aceleracion():
    Tiempo = int(Entrada_ti.get())
    ti=t[Tiempo]
    x_a = -a*np.cos(ti)
    y_a = (b*np.tan(ti/2))/(np.cos(ti/2)**2)
    z_a = -a*np.sin(ti)
    ace = "(" + str(round(x_a, 3)) + "," + str(round(y_a, 3)) + "," + str(round(z_a, 8)) + ")"
    Var.set((ace))

def velocidad():
    Tiempo = int(Entrada_ti.get())
    ti=t[Tiempo]
    x_v = -a*np.sin(ti)
    y_v = (b*(1/(np.cos(ti)**2)))/2
    z_v = a*np.cos(ti)
    vel = "("+str(round(x_v,3))+","+str(round(y_v,3))+","+str(round(z_v,8))+")"
    Var.set(vel)

def velocidad_media():
    TiempoI = int(Entrada_ti.get())
    TiempoF = int(Entrada_tf.get())
    ti = t[TiempoI]
    tf = t[TiempoF]
    x_vi = -a*np.sin(ti)
    y_vi = b/(2*(np.cos(ti/2)**2))
    z_vi = a*np.cos(ti)

    x_vf = -a*np.sin(tf)
    y_vf = b/(2*(np.cos(tf/2)**2))
    z_vf = a*np.cos(tf)
    x_v = (x_vf - x_vi)/(tf-ti)
    y_v = (y_vf - y_vi)/(tf-ti)
    z_v = (z_vf - z_vi)/(tf-ti)
    vel = "(" + str(round(x_v, 8)) + "," + str(round(y_v, 8)) + "," + str(round(z_v, 8)) + ")"
    Var.set((vel))

def aceleracion_media():
    TiempoI = int(Entrada_ti.get())
    TiempoF = int(Entrada_tf.get())
    ti = t[TiempoI]
    tf = t[TiempoF]
    x_ai = -a*np.cos(ti)
    y_ai = (b*np.tan(ti/2))/(np.cos(ti/2)**2)
    z_ai = -a*np.sin(ti)

    x_af = -a*np.cos(tf)
    y_af = (b*np.tan(tf/2))/(np.cos(tf/2)**2)
    z_af = -a*np.sin(tf)

    x_a = (x_af - x_ai)/(tf-ti)
    y_a = (y_af - y_ai)/(tf-ti)
    z_a = (z_af - z_ai)/(tf-ti)
    ace = "(" + str(round(x_a, 8)) + "," + str(round(y_a, 8)) + "," + str(round(z_a, 8)) + ")"
    Var.set((ace))

def curvatura():
    Tiempo = int(Entrada_ti.get())
    ti = t[Tiempo]
    x_pder = -a*np.sin(ti)
    y_pder = b/(2*(np.cos(ti/2)**2))
    z_pder = a*np.cos(ti)

    x_sder = -a*np.cos(ti)
    y_sder = (b*np.tan(ti/2))/(np.cos(ti/2)**2)
    z_sder = -a*np.sin(ti)
    cruz_x = [(y_pder*z_sder)-(y_sder*z_pder),-((x_pder*z_sder)-(x_sder*z_pder)),
                      (x_pder*y_sder)-(x_sder*y_pder)]
    r = np.sqrt((cruz_x[0])**2+(cruz_x[1])**2+(cruz_x[2])**2)
    s = (np.sqrt((x_pder**2)+(y_pder**2)+(z_pder**2)))**3
    curv = round(r/s,8)
    Var.set(curv)

def radio_curvatura():
    Tiempo = int(Entrada_ti.get())
    ti = t[Tiempo]
    #derivadas
    x_pder = -a*np.sin(ti)
    y_pder = b/(2*(np.cos(ti/2)**2))
    z_pder = a*np.cos(ti)

    x_sder = -a*np.cos(ti)
    y_sder = (b*np.tan(ti/2))/(np.cos(ti/2)**2)
    z_sder = -a*np.sin(ti)
    cruz_x = [(y_pder * z_sder) - (y_sder * z_pder), -((x_pder * z_sder) - (x_sder * z_pder)),
              (x_pder * y_sder) - (x_sder * y_pder)]
    r = np.sqrt((cruz_x[0])**2+(cruz_x[1])**2+(cruz_x[2])**2)
    s = (np.sqrt((x_pder**2)+(y_pder**2)+(z_pder**2)))**3
    curv = r/s
    radio_curv = round(1/curv,8)
    Var.set(radio_curv)

def torsion():
    Tiempo = int(Entrada_ti.get())
    ti = t[Tiempo]
    #derivadas
    x_pder = -a * np.sin(ti)
    y_pder = b / (2 * (np.cos(ti / 2) ** 2))
    z_pder = a * np.cos(ti)

    x_sder = -a * np.cos(ti)
    y_sder = (b * np.tan(ti / 2)) / (np.cos(ti / 2) ** 2)
    z_sder = -a * np.sin(ti)

    x_tder = a*np.sin(ti)
    y_tder = (b*(np.tan(ti/2)**2)*(1/np.cos(ti/2))+b*(1/(np.cos(ti/2)**3)))/4
    z_tder = -a*np.cos(ti)
    cruz_x = [(y_pder * z_sder) - (y_sder * z_pder), -((x_pder * z_sder) - (x_sder * z_pder)),
              (x_pder * y_sder) - (x_sder * y_pder)]
    r = np.sqrt((cruz_x[0])**2 + (cruz_x[1])**2 + (cruz_x[2])**2)**2
    cruz_numerador=[(y_sder * z_tder) - (y_tder * z_sder), -((x_sder * z_tder) - (x_tder * z_sder)),
              (x_sder * y_tder) - (x_tder * y_sder)]
    r_torsion= round((cruz_numerador[0]*x_pder+cruz_numerador[1]*y_pder+cruz_numerador[2]*z_pder)/r,3)
    Var.set(r_torsion)

def radio_torsion():
    Tiempo = int(Entrada_ti.get())
    ti = t[Tiempo]
    #derivadas
    x_pder = -a * np.sin(ti)
    y_pder = b / (2 * (np.cos(ti / 2) ** 2))
    z_pder = a * np.cos(ti)

    x_sder = -a * np.cos(ti)
    y_sder = (b * np.tan(ti / 2)) / (np.cos(ti / 2) ** 2)
    z_sder = -a * np.sin(ti)

    x_tder = a*np.sin(ti)
    y_tder = (b*(np.tan(ti/2)**2)*(1/np.cos(ti/2))+b*(1/(np.cos(ti/2)**3)))/4
    z_tder = -a*np.cos(ti)
    cruz_x = [(y_pder * z_sder) - (y_sder * z_pder), -((x_pder * z_sder) - (x_sder * z_pder)),
              (x_pder * y_sder) - (x_sder * y_pder)]
    r = np.sqrt((cruz_x[0])**2 + (cruz_x[1])**2 + (cruz_x[2])**2)**2
    cruz_numerador=[(y_sder * z_tder) - (y_tder * z_sder), -((x_sder * z_tder) - (x_tder * z_sder)),
              (x_sder * y_tder) - (x_tder * y_sder)]
    r_torsion= round((cruz_numerador[0]*x_pder+cruz_numerador[1]*y_pder+cruz_numerador[2]*z_pder)/r,8)
    radio_t = 1/r_torsion
    Var.set(radio_t)

def posicion():
    Tiempo = int(Entrada_ti.get())

    X = round(XYZ[0][Tiempo],8)
    Y = round(XYZ[1][Tiempo],8)
    Z = round(XYZ[2][Tiempo],8)
    Axes = "("+str(X)+","+str(Y)+","+str(Z)+")"
    Var.set(Axes)

def desbloquear_botones():
    boton_posicion.config(state="normal")
    boton_velocidad.config(state="normal")
    boton_aceleracion.config(state="normal")
    boton_velocidad_media.config(state="normal")
    boton_aceleracion_media.config(state="normal")
    boton_curvatura.config(state="normal")
    boton_radio_curvatura.config(state="normal")
    boton_torsion.config(state="normal")
    boton_radio_torsion.config(state="normal")
    boton_longitud_arco.config(state="normal")

def presionado(boton):
    desbloquear_botones()
    if(boton=="posicion"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_posicion.config(state="disabled")
        Boton_R.config(state="normal",command=lambda : posicion())
    elif (boton == "velocidad"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_velocidad.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: velocidad())
    elif (boton == "aceleracion"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_aceleracion.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: aceleracion())
    elif (boton == "velocidad media"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="normal")
        boton_velocidad_media.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: velocidad_media())
    elif (boton == "aceleracion media"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="normal")
        boton_aceleracion_media.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: aceleracion_media())
    elif (boton == "curvatura"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_curvatura.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: curvatura())
    elif (boton == "radio curvatura"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_radio_curvatura.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: radio_curvatura())
    elif (boton == "torsion"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_torsion.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: torsion())
    elif (boton == "radio torsion"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="disabled")
        boton_radio_torsion.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: radio_curvatura())
    elif (boton == "longitud de arco"):
        Entrada_ti.config(state="normal")
        Entrada_tf.config(state="normal")
        boton_longitud_arco.config(state="disabled")
        Boton_R.config(state="normal", command=lambda: radio_curvatura())

if __name__ == '__main__':

    ventana = tk.Tk()
    ventana.geometry("900x900")
    ventana.title("Horoptera")
    #creamos el frame contenedor de los elementos
    frame = tk.Frame(ventana)
    frame.grid(column=0, row=0)

    #creamos el grafico
    plt.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    frame_grafico = FigureCanvasTkAgg(fig, master=frame)
    frame_grafico.get_tk_widget().grid(column=1, row=1, columnspan=6, rowspan=9)

    ax = fig.gca(projection='3d')
    a = 3  # corresponde al radio
    b = 1  # es una constante
    t = np.linspace(0.3 - np.pi, -0.3 + np.pi, 102)
    x = a + a * np.cos(t)
    y = b * np.tan(t / 2)
    z = a * np.sin(t)
    XYZ = [x, y, z]
    XYZ = np.array(list(XYZ))
    horoptera, = ax.plot(x[0:1], y[0:1], z[0:1])

    ax.set_xlim3d([-8.0, 8.0])
    ax.set_xlabel('Eje X')
    ax.set_ylim3d([-8.0, 8.0])
    ax.set_ylabel('Eje Y')
    ax.set_zlim3d([-8.0, 8.0])
    ax.set_zlabel('Eje Z')
    animacion = animation.FuncAnimation(fig, sig, fargs=(XYZ, horoptera), frames=100, blit=False, interval=1,
                                        repeat=False)
    #creando Entradas para los numeros
    frame_entradas_texto = tk.Frame(frame)
    frame_entradas_texto.grid(column=1,row=10)
    texto_ti = tk.Label(frame_entradas_texto,text="Ti")
    texto_ti.grid(column=0,row=0)
    Entrada_ti = tk.Entry(frame_entradas_texto,state="disabled")
    Entrada_ti.grid(column=1,row=0,padx=(0, 50))
    texto_tf = tk.Label(frame_entradas_texto, text="Tf")
    texto_tf.grid(column=2, row=0)
    Entrada_tf = tk.Entry(frame_entradas_texto, state="disabled")
    Entrada_tf.grid(column=3, row=0,padx=(0, 50))
    Var = tk.StringVar()
    Boton_R = tk.Button(frame_entradas_texto, text="Resp =",state="disabled")
    Boton_R.grid(column=4, row=0)
    Res = tk.Label(frame_entradas_texto,textvariable=Var)
    Res.grid(column=5, row=0)

    #BOTONES--------------------------------------

    titulo = tk.Label(frame, text="Horóptera",font=("letter case", 50))
    titulo.grid(column=1, row=0, pady=40,columnspan=5,)

    boton_posicion = tk.Button(frame, text="Posicion", width=20, height=2, command=lambda : presionado("posicion"))
    boton_posicion.grid(column=0, row=1)

    boton_velocidad = tk.Button(frame, text="Velocidad", width=20, height=2, command=lambda : presionado("velocidad"))
    boton_velocidad.grid(column=0, row=2)

    boton_velocidad_media = tk.Button(frame, text="Velocidad Media", width=20, height=2,
                                      command=lambda : presionado("velocidad media"))
    boton_velocidad_media.grid(column=0, row=3)

    boton_aceleracion = tk.Button(frame, text="Aceleracion", width=20, height=2,
                                  command=lambda : presionado("aceleracion"))
    boton_aceleracion.grid(column=0, row=4)

    boton_aceleracion_media = tk.Button(frame, text="Aceleracion Media", width=20, height=2,
                                        command=lambda : presionado("aceleracion media"))
    boton_aceleracion_media.grid(column=0, row=5)

    boton_curvatura = tk.Button(frame, text="Curvatura", width=20, height=2,
                                command=lambda : presionado("curvatura"))
    boton_curvatura.grid(column=0, row=6)

    boton_radio_curvatura = tk.Button(frame, text="Radio de Curvatura", width=20, height=2,
                                      command=lambda : presionado("radio curvatura"))
    boton_radio_curvatura.grid(column=0, row=7)

    boton_torsion = tk.Button(frame, text="Torsion",width=20, height=2,
                              command=lambda : presionado("torsion"))
    boton_torsion.grid(column=0, row=8)

    boton_radio_torsion = tk.Button(frame, text="Radio de Torsion", width=20, height=2,
                                    command=lambda : presionado("radio torsion"))
    boton_radio_torsion.grid(column=0, row=9)

    boton_longitud_arco = tk.Button(frame, text="Longitud de Arco",width=20, height=2,
                                    command=lambda : presionado("longitud de arco"))
    boton_longitud_arco.grid(column=0, row=10)

    boton_ani = tk.Button(frame, text="Ani",width=20, height=2)
    boton_ani.grid(column=6, row=0)
    boton_atras = tk.Button(frame, text="<",width=20, height=2,command=lambda : Atras(ventana))
    boton_atras.grid(column=0, row=0)
    tk.mainloop()
