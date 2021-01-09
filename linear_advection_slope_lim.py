import numpy as np
import matplotlib.pyplot as plt
import time



    
def get_user_opinion(question):
    while True:
        method = input(question)
        if method == "1":
            return "cent_dif"
        if method == "2":
            return "upwind"
        if method == "3":
            use_method = "lax"
            return "lax"
        if method == "4":
            use_method = "lax_min"
            return "lax_min"
        if method == "5":
            use_method = "lax_sup"
            return "lax_sup"
        else:
            print("please enter 1,2,3,4 or 5")
    
        
def minmod(a,b):
    if abs(a)<abs(b) and a*b>0:
        return a
    if abs(b)<=abs(a) and a*b>0:
        return b
    if a*b<=0:
        return 0
    
def maxmod(a,b):
    if abs(a)<abs(b) and a*b>0:
        return b
    if abs(b)<=abs(a) and a*b>0:
        return a
    if a*b<=0:
        return 0
    
    


def make_step(use_method,u,timestep):
    #Randbedingungen Anwenden
    u[0] = u[N]
    u[N+1] = u[1]
    
    
    #u zum neuen Zeitschritt ausrechnen
    u_new=np.zeros_like(u)
    
    if use_method == "cent_dif":
        for i in range (1,N+1):
            u_new[i] = u[i] - timestep/(delta_x*2)*(u[i+1]-u[i-1])
        return u_new
        
    
    if use_method == "upwind":  
        for i in range (1,N+1):
            u_new[i] = u[i] - timestep/delta_x*(u[i]-u[i-1])       
        return u_new
    
    
    if use_method == "lax":
        for i in range(1,N+1):
            u_new[i] = u[i] - sigma/2*(u[i+1]-u[i-1])+sigma**2/2*(u[i+1]-2*u[i]+u[i-1])
        return u_new
    
    if use_method == "lax_min":
        for i in range(1,N+1):
            u_new[i] = u[i] - sigma*(u[i]-u[i-1])-sigma/2*(delta_x-timestep)*((minmod((u[i]-u[i-1])/delta_x,(u[i+1]-u[i])/delta_x)-minmod((u[i-1]-u[i-2])/delta_x,(u[i]-u[i-1]))/delta_x))
        return u_new
    
    if use_method == "lax_sup":
        for i in range(1,N+1):
            u_new[i] = u[i] - sigma*(u[i]-u[i-1])-sigma/2*(delta_x-timestep)*(maxmod(minmod((u[i+1]-u[i])/delta_x,2*(u[i]-u[i-1])/delta_x),(minmod(2*(u[i+1]-u[i])/delta_x,(u[i]-u[i-1])/delta_x)))-maxmod(minmod((u[i]-u[i-1])/delta_x,2*(u[i-1]-u[i-2])/delta_x),minmod(2*(u[i]-u[i-1])/delta_x,(u[i-1]-u[i-2])/delta_x)))
        return u_new     
    
    
        
def make_steps(use_method, u,timestep,max_time):
    time_passed = 0
    while time_passed < max_time:
        time_passed += timestep
        u = make_step(use_method,u,timestep)
    return u

def plot_grid(u,x,label):
    """Plottet ein Grid im physikalisch relevanten Teil
    grid: das Gitter
    val_name: Name der Größe als plotting label"""
    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(x[2:len(x)-4],u[2:len(u)-4],label=label)

def analytical_sol(t,n):
    """Gibt ein array der Länge N+4 mit der Analytischen Lösung zurück"""
    arr = np.linspace(-1,1,n)
    arr -= t
    arr = (arr+1) % 2 -1
    return abs(arr) <= 1/3

print("Program for solving the linear advection using different methods")
N= int(input("How Many Gridpoints? > "))
delta_x = 2/N
sigma = 0.8
timestep = sigma*delta_x
max_time = int(input("What shall be the maximum time?: "))
method = get_user_opinion("Which method shall be used?\n\n[1] Centered Differencing \n[2] Upwind \n[3] Lax-Wendrof\n[4] Lax-Wendrof with minmod\n[5] Lax-Wendrof with superbee\n")
u = np.zeros(N+2)
x = np.zeros(N+2)
for i in range(1,N+1):
    x[i]=2*i/N-1
    if abs(x[i])<=1/3:
        u[i]=1
begin_time = time.time()
u_new=make_steps(method,u,timestep,max_time)
print("The computation took {0:7.2f} ms".format((time.time()-begin_time)*1000))
plot_grid(u_new,x,"Simulation")
n = 1000
plot_grid(analytical_sol(max_time,n),np.linspace(-1,1,n),label="Analytical")
plt.legend()
plt.show()
