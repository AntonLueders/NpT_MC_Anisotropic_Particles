Here, we present a simple example on how to calculate the equation of state (EOS) for hexagons that are horizontally deformed as characterized by $\delta = 0.5$. As in [1], we want to study a system where the area of the particles $A_k(\delta)$ is set to the area of a reference system of circles described by a diameter $\sigma = 1$.

**1) Calculate the parameters**

First, we have to determine the parameters of the system. For this, we can check [Parameters.pdf](../Parameters.pdf) to find the input needed for a simulation of hexagons. We can identify that we have to specify the number of vertices $m$, the length of the edges $a$, and the scaling parameter $\lambda$. Since we analyze hexagons, we have $m = 6$. Also, since we want to stretch the particles horizontally, we know that $\lambda > 1.0$.

From [Parameters.pdf](../Parameters.pdf), we find $\lambda = \delta + 1 = 1.5$. The area of an arbitrary regular hexagon that is stretched by the factor $\lambda = 1.5$ is given by $A_{\textrm{hex}} = 3\sqrt{3} a^2 \lambda / 2$. Setting this relation equal to the particle area of the reference system $\pi/4$, we find $a \approx 0.4489$.

**2) Downloand and compile the code**

Since we want to study hexagons, we download the code for Polygons. Then, one can compile the code using
``
make
``
in a console opened at the path where the code is saved.

**3) Run the code**

Now, we have to prepare the input file. For the horizontally stretched hexagons with $\delta = 0.5$, we can download the example given in this folder. Before running the code, we must modify it by adding the derived values for $a$, $m$, and $\lambda$ (they are currently left blank in the example input file). The meaning of the other parameters is omitted here. Please check the README and the code for more information.

To obtain the EOS, we have to run the program multiple times. For each run, we modify the predetermined pressure (since we run $NpT$ simulations). This can be done via the **p** parameter line. Note that the system must be equilibrated to obtain meaningful data.

**4) Extract the data**

After equilibration, the total area of the simulation box should fluctuate around a fixed mean value. The last column of the output file *Position\_For\_Ovito.xyz* is the box length $L$ of the simulation box at the time of the particular simulation step. Note that this value should be the same for all particles. We can use this to determine the packing fraction for the studied pressure value via $\eta = \pi \langle N/L^2 \rangle / 4$, where $N$ is the number of simulated particles (which is 2000 in the example input file). Since the simulation runs in units of $k_{\mathrm{B}}T$, the pressure normalized by the corresponding ideal gas pressure can be calculated by $p / \langle N/L^2 \rangle$, where $p$ is the number as set in the input file. Plotting this quotient over $\eta$ gives the EOS.
