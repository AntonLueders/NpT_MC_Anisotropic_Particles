Here, we present a simple example on how to calculate the EOS for hexagons that are horizontally deformed as characterized by $\delta = 0.5$. As in [1], we want to study a system, where the area of the particles $A_k(\delta)$ is set to the area of a reference system of circles described by a diameter $\sigma=1$. 

**1) Calculate the parameters**

First, we have to determine the parameters of the system. For this, we can check [Parameters.pdf](../Parameters.pdf) to find the input needed for a simulation of hexagons. We can identify that we have to specify the number of vertices $m$, the length of the edges $a$, and the scaling parameter $\lambda$. Since we analyze hexagons, we have $m=6$. Also, since we want to stracht the partically horizontally, we know that $\lambda > 1.0$.

From [Parameters.pdf](../Parameters.pdf), we find $\lambda=\delta+1=1.5$. The area of a regular hexagon that is streched by the factor $\lambda=1.5$ is given by $A_{\textrm{hex}}=3\sqrt{3}a^2\lambda/2$. Setting this relation equals to particle area of the reference $\pi/4$, we find $a \approx 0.4489$.

**2) Downloand and compile the code**

Since we want to study hexagons, we download the code for Polygons. Then, one can compile the code using
``
make
``
in a console opened at the bath where the code is saved.

**3) Run the code**

Now, we have to prepare the inputfile. For the horizontally streched hexagons with $\delta=0.5$, we can downlad the example given in this folder. 
