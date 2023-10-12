[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/SlWp5koN)
# Physics 427 Homework 6

__Due 11:59pm Wednesday 10/11/2023__

There is only one problem in this assignment.

## 1. The Quantum Anharmonic Oscillator

The time-independent Schrödinger equation in 1D reads:

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
$$

where $\psi$ is the wave function with $|\psi(x)|^2$ describing the probability of finding the particle at position $x$, $\hbar$ is the reduced Planck constant, $m$ is the mass of the particle, and $E$ is the energy eigenvalue. The function $V(x)$ specifies a potential energy, and its shape determines the allowable energies $E$ and the corresponding wave functions. In general, the allowable energies are discrete, hence the name "Quantum Mechanics".

Let's first make the equation dimensionless. In quantum mechanics, a natural length scale to use is the reduced Compton wavelength $ƛ = \hbar/mc$. For an electron, $ƛ = 3.86\times 10^{-13}\,\mathrm{m}$. We can then define the _dimensionless length_ $\hat{x}$ using $x = ƛ\hat{x}$. The 1D time-independent Schrödinger equation now reads (after dividing both sides by $mc^2$):

$$
-\frac{1}{2}\frac{d^2\psi}{d\hat{x}^2} + \frac{V}{mc^2}\psi = \frac{E}{mc^2}\psi
$$

You probably recognize that $mc^2$ is simply the rest mass energy of the particle. It is convenient to normalize energies (both $E$ and the potential $V$) to $mc^2$. So we define the dimensionless $\hat{V}$ and $\hat{E}$ using:

$$
V = \hat{V}mc^2,\quad E = \hat{E}mc^2
$$

Then we recover a dimensionless time-indepenent Schrödinger equation that is much simpler:

$$
-\frac{1}{2}\frac{d^2\psi}{d\hat{x}^2} + \hat{V}\psi = \hat{E}\psi
$$

We would like to look at the following potential which defines the 1D quantum anharmonic oscillator (written in a dimensionless way):

$$
\hat{V}(\hat{x}) = \frac{1}{2}\hat{k}_0\hat{x}^2 + \hat{k}_1 \hat{x}^4
$$

This is a variation on the ordinary quantum harmonic oscillator that you might have seen, where $V(x) = kx^2/2$. Unlike the harmonic oscillator, this variation does not have a closed analytic solution. We can use numerical methods to find the energy eigenvalues and the corresponding wave functions.

Let's first consider $\hat{k}_1 = 0$ which reduces to the harmonic oscillator case. The frequency of the harmonic oscillator is $\omega = \sqrt{k_0/m}$, and the ground state energy is simply $\hbar\omega/2$. If we normalize the ground state energy by $mc^2$, it can be shown that the dimensionless ground state energy is:

$$
\hat{E}_0 = \frac{1}{2}\frac{\hbar\omega}{mc^2} = \frac{1}{2}\sqrt{\hat{k}_0}
$$

When $0 <\hat{k}_1 \ll \hat{k}_0$, the solution should be close enough to the ordinary harmonic oscillator, so that the ground state energy is not too different. This gives you a sense of what kind of energy to expect.

In a C++ file `problem1.cpp`, use the shooting method to solve the anharmonic oscillator with $\hat{k}_0 = 10^{-4}$ and $\hat{k}_1 = 10^{-7}$. Technically, the wave function $\psi(x)$ is defined from $-\infty$ to $\infty$. However, we can simply use a reduced interval such as $[-a, a]$, and let $\psi(-a) = \psi(a) = 0$ be your boundary conditions. For example, $a = 30$ is sufficient for our problem. You can introduce a new function $\phi$ to reduce the order of the ODE to 1:

$$
\frac{d\psi}{d\hat{x}} = \phi,\quad \frac{d\phi}{d\hat{x}} = 2(\hat{V} - \hat{E})\psi
$$

Even though we have 2 equations, there is still only one condition we need to satisfy at $\hat{x} = a$, which is $\psi(a) = 0$. Use RK4 to solve the initial condition problem and use Newton's method to find the ground state energy $\hat{E}$. The goal of Newton's method here is to find the eigenvalue $\hat{E}$ such that $\psi(a) = 0$, and you can use a small increment $h_a$ on $\hat{E}$ to compute the numerical derivative required in Newton's method. Due to the use of numerical derivative, you won't be able to directly reuse what you wrote in Homework 1. To make sure you get to the ground state, you can start your Newton's method using a trial solution of $\hat{E} = 0$. Your program should print this energy:

``` sh
E = [the ground-state energy you find]
```

You will probably need $\phi(-a)$ as part of your initial condition. Use a small value such as $\phi(-a) \sim 10^{-3}$. The exact value shouldn't matter since it only determines the overall normalization of the wave function. Use a small increment size $h_a$ for your Newton's method, such as $h_a \sim 10^{-6}$. Again the exact value doesn't matter but you need to make sure it's much smaller than the energy you are looking for.

Write the output of your program (one-liner containing the ground-state energy) to `problem1.txt` and commit it to the repo. In addition, plot the ground state wave function in Python and include the plot in the repo as `problem1.png`.
