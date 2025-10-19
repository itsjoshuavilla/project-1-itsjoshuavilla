[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/0km5Hvvm)
# Project1

# Project 1 Details
Important: This project is due 03/23/2025 at 11:59 PM.  No late submissions will be accepted.  

This program will be written in the C++ language, and recap pointers, structures, and dynamic memory allocation.

Create a program that will solve for voltage and current across a capacitor as a function of time.  

# Constant current source power supply
For charging a capacitor with a constant current power supply, the resulting circuit looks like

![cc](https://user-images.githubusercontent.com/123007812/221931749-494ce0b3-fc32-403b-9437-77c5e375d0f7.png)


The current through a capacitor as a function of time can be written as

```math
I(t) = C \frac{dV(t)}{dt}    
```

Using a finite-difference method, the voltage across the capacitor can be found as 

```math
V(t+1) = V(t) + I(t) \Delta t  \frac{1}{C}   
```

At time t = 0, the voltage across the capacitor is 0.

# Constant voltage source power supply
In this case, a resistor is used to protect the power supply.  The resulting circuit looks as

![rc](https://user-images.githubusercontent.com/123007812/221923282-443cbb78-7e5c-4881-9b94-7a3d7dae092e.PNG)
(source wikipedia)

Using KVL, the voltage, $V_0$ has to be equal to the voltage drop across $R$, $V_r$, and $C$, $V_c$.  This can be written as, 

```math
V_0(t) = V_r(t) + V_c(t)   
```

Integrating equation (1) and plugging into equation (3) with Ohm's law gi9ves
```math
V_0(t) = I(t)R + \frac{1}{C} \int I(\tau) d\tau
```

Taking the derivative of this equation and multiplying by C gives

```math
0 = RC\frac{dI(t)}{dt} + I(t)
```

Using finite differences, the current in the circuit, and through the capacitor is then

```math
I(t+1) = I(t) - \frac{I(t)}{RC} \Delta t
```

At time t = 0, $I(t=0)$ is equal to $\frac{V_0}{R}$.

# Data Format
Use a structure that looks like
```
struct _capacitor
{
  double *time;       // time array
  double *voltage;    // voltage array
  double *current;    // current array
  double C;           // capacitance value
};
typedef struct _capacitor Capacitor
```

Solve for the voltage across the capacitor in the constant current case using 

```math
V(t+1) = V(t) + I(t) \Delta t  \frac{1}{C}   
```

In words, this means that the voltage at the next time step (t+1) is equal to the voltage and the previous time step (t) plus the current at the previous timestep (t) multiplied by $\Delta t$ and inverse capacitance.

and the current through the capacitor using 

```math
I(t+1) = I(t) - \frac{I(t)}{RC} \Delta t
```
Again, in words, this means that the current at the next timestep (t+1) is equal to the current at the previous timestep (t-1) minus the current at the previous timestep (t) multiplied by resistance, capacitance and $\Delta t$.

Use the values of 
  - $\Delta t$ = 1e-10 s
  - Final time = 5e-6 s
  - Number of timesteps = 5e-6/1e-10 = 50000
  - $R$ = 1 $k\Omega$
  - $C$ = 100e-12 F
  - $I(t)$ = 1e-2 A (for constant current supply)
  - $V_0$ = 10.0 V  (for constant voltage supply)
  
You will need to dyanmically allocate the time, voltage, and current arrays in the structure.  Output the resulting current and voltage for the capacitor for the different charging circuits every 200 timesteps to the screen.

A loop over the number of timesteps will be needed to iterate the voltage and current in time.  Start at timestep 1, with timestep 0 being filled with the initial conditions given in the equation details.
  
