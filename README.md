# The frame of Finite Element Method 

- Model problem:

$$
-\nabla\cdot(p~\nabla u)+\mathbf{q}\cdot\nabla u+ru=f,~u\in \Omega,
$$

and its Dirichlet boundary condition: $u=0,~u\in\partial \Omega$.

- Variational form: Find $u\in H_0^1(\Omega)=\{\omega\in H^1(\Omega),\omega|_{\partial \Omega}=0\}$, s.t.

$$
a(u,v)=f(v),\quad \forall v\in H_0^1(\Omega),
$$

where $a(u,v)=(p\nabla u,\nabla v)+(\mathbf{q}\cdot\nabla u,v)+(ru,v),\quad f(v)=(f,v).$

- And we use Galerkin method to discrete the infinite dimensional functional space.



[TOC]



## Main frame

Here, we mainly focus on the standard model problem, just like $\Omega =[0,1]^n$, and $n$ is the dimension.



The programming of the FEM follows:

### PDE information input

we will input the **interval** of PDE model, **degree** of basis function and **number of division** in each axis.

- take the two dimension as an example:

```matlab
pde.left = 0;
pde.right = 1;
pde.bottom = 0;
pde.top = 1;

deg = 1;% degree of basis function

N1 = 2;% number of division on X-axis
N2 = 2;% number of division on Y-axis
```



### Mesh matrix generation

Here we use $P,T,P_b,T_b$ four matrix to record mesh information.

And $P-n\times nofnode$: record the coordinates of each node in real mesh.

​		 $T-nofeletype\times nofelement$: record the node number of each element in real mesh.

​		 $P_b-n\times nofelmnode$: record the coordinates of each node in FEM mesh.

​		 $T_b-nofeletype\times nofemelement$: record the node number of each element in FEM mesh.

- take the two dimension and one degree of basis function as an example:

  - triangle mesh: ==$h$ is the radius of circumcircle==

    ![TriangleMesh](https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/T_triangle.jpg)

    And its $P$ and $T$ matrix:<img src="https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/P_triangle.jpg" alt="P_triangle" style="zoom:80%;" />

    ***

    ![T_triangle](https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/T_triangle.jpg)

  - rectangle mesh:

    ![RectangleMesh](https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/RectangleMesh.jpg)

    And its $P$ and $T$ matrix:

    <img src="https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/P_rectangle.jpg" alt="P_rectangle" style="zoom:80%;" />

    ***

    <img src="https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/T_rectangle.jpg" alt="T_rectangle" style="zoom:80%;" />

here $P_b=P,~T_b=T$.

==Attention: $P_b$ and $T_b$ is different from $P$ and $T$ when $deg\neq 1$.==

- 3D cubic mesh:

  <img src="https://github.com/SuperMiffy/FEM_Poisson/blob/main/figure/CubeMesh.png" alt="CubeMesh" style="zoom: 80%;" />

### Stiffness matrix

In this part, we mainly focus on the different part of stiffness matrix and we use local mass matrix to assemble the all stiffness matrix by **element loop**.

```matlab
%% Generate the stiffness matrix
A = StiffnessMatrix(P,Pb,T,Tb,'p',deg,1,0,1,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'p',deg,0,1,0,1)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',deg,1,0,0,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',deg,0,1,0,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'r',deg,0,0,0,0);
```



### Load vector

Assemble the load vector is easier compared to the stiffness matrix.

```matlab
%% Generate the load vector
b = LoadVector(P,Pb,T,Tb,deg,0,0);
```



### Boundary condition

Here, we use ==**Boundary Matrix($3\times nofboundarynode$)**== to record the node information on the boundary. 

- for each column:
  - BoundaryMatrix(1,i) -- boundary type(-1: Dirichlet,  0: Neumann,  1: Robin);
  - BoundaryMatrix(2,i) -- current boundary node number;
  - BoundaryMatrix(3,i) -- **boundary value**.(==attention==)

And indeed, here is the standard model problem, and we can easily get the boundary node by **for loop**.

#### Dirichlet Boundary condition

- make the BoundaryMatrix(3,i) to be equal to zero and treat linear systems.

#### Neumann Boundary condition

- calculate the right handside value and put it on BoundaryMatrix(3,i).

#### Robin Boundary condition

- treat the corresponding stiffness matrix and RHS(RightHandSide) load vector at the same time.



### Solve the linear system

we need to pay attention to the **characteristic of stiffness matrix** and choose appropriate numerical method to do solution.

here the stiffness matrix is always tridiagonal matrix and sparse.

- A/b   bigstable



### Error estimate

#### $L^2$-norm error

$$
\left\|u-u_{h}\right\|_{0}=\sqrt{\int_\Omega(u-u_h)^2\mathrm{d}x\mathrm{d}y}=\sqrt{\sum_{n=1}^{N} \int_{E_{n}}\left(u-\sum_{k=1}^{N_{l b}} u_{T_{b}(k, n)} \psi_{n k} .\right)^{2} \mathrm{d}x\mathrm{d}y}.
$$

#### $H^1$-semi-norm error

$$
\left|u-u_{h}\right|_{1, x}=\sqrt{\int_{\Omega}\left(\frac{\partial\left(u-u_{h}\right)}{\partial x}\right)^{2}\mathrm{d}x\mathrm{d}y}=\sqrt{\sum_{n=1}^{N} \int_{E_{n}}\left(\frac{\partial u}{\partial x}-\sum_{k=1}^{N_{l b}} u_{T_{b}(k, n)} \frac{\partial \psi_{n k}}{\partial x}\right)^{2} \mathrm{d} x \mathrm{d} y}.
$$

And the remaining partial derivatives is the same way.

#### Error order

error order is the powerful verification tool. 

here we always use $\log$ function to estimate the error order.

- e.g.: we have the theoretical result: $\left\|u-u_h\right\|_a = error_h\leqslant Ch^p\left\|u\right\|_b.$ 

​				and we do the bisection: $error_{h/2}\leqslant C(h/2)^p \left\|u\right\|_b.$

​				thus we can get:
$$
p = \log_2(error_h/error_{h/2}).
$$

- we need to notice that we always do many of bisection and take the average.



## Preprocessing function

### Basis function

Here, we main use affine transformation to generate the basis function from reference element onto real element.

#### One dimension

- affine transformation: [-1,1] -> [a,b]
  $$
  \hat{x} = f(x)=\frac{b-a}{2}x+\frac{b+a}{2},\quad x\in [-1,1].
  $$

- Basis function -- use **Lagrange Interpolation** to generate the basis polynomial function.

#### Two dimension

- affine transformation:==pay attention to the chain rule for partial derivatives==
  - triangulation mesh:
  
    reference element:  △[(0,0),(1,0),(0,1)]   -->  real element: △[$(x_1,y_1),(x_2,y_2),(x_3,y_3)$].
  
    - affine transformation: $\mathbf{x} = T\mathbf{\hat{x}}+R$,(from reference to real) and here
  
    $$
    T=\left[\begin{array}{ll}
    x_2-x_1 & x_3-x_1 \\
    y_2-y_1 & y_3-y_1
    \end{array}\right],\quad R=\left[\begin{array}{ll}
    x_1 \\
    y_1
    \end{array}\right].
    $$
  
    - ==chain rule for derivatives:==
  
      here we want to calculate the derivatives on real element:
      $$
      \frac{\partial F(x,y)}{\partial x} = \frac{\partial F(\hat{x},\hat{y})}{\partial\hat{x}}\frac{\partial\hat{x}(x,y)}{\partial x}+\frac{\partial F(\hat{x},\hat{y})}{\partial\hat{y}}\frac{\partial\hat{y}(x,y)}{\partial x}.
      $$
      
  
      - hence we calculate the affine inverse transformation:
        $$
        \mathbf{\hat{x}}=\left[\begin{array}{ll}
        \hat{x} \\
        \hat{y}
        \end{array}\right]\xlongequal{\Delta}T^{-1}(\mathbf{x}-R)=\frac{1}{|T|}\left[\begin{array}{ll}
        y_3-y_1 & x_1-x_3 \\
        y_1-y_2 & x_2-x_1
        \end{array}\right]\left[\begin{array}{ll}
        x-x_1 \\
        y-y_1
        \end{array}\right].
        $$
        so we note: $\hat{x}\xlongequal{\Delta}T_1(x,y),~\hat{y}\xlongequal{\Delta}T_2(x,y)$.
  
      - Thus:
  
        - $$
          \frac{\partial F(x,y)}{\partial x} = \frac{\partial F(\hat{x},\hat{y})}{\partial\hat{x}}\frac{\partial T_1(x,y)}{\partial x}+\frac{\partial F(\hat{x},\hat{y})}{\partial\hat{y}}\frac{\partial T_2(x,y)}{\partial x}.
          $$
  
        - $$
          \frac{\partial F(x,y)}{\partial y} = \frac{\partial F(\hat{x},\hat{y})}{\partial\hat{x}}\frac{\partial T_1(x,y)}{\partial y}+\frac{\partial F(\hat{x},\hat{y})}{\partial\hat{y}}\frac{\partial T_2(x,y)}{\partial y}.
          $$
  
- 



- Basis function: ==here the index order of basis function should be corresponding with the order in $T$ matrix.==
  - triangle reference element △[(0,0),(1,0),(0,1)]
    - node (0,0): $\phi_1(x,y)=1-x-y$ ;
    - node (1,0): $\phi_2(x,y)=x$ ;
    - node (0,1): $\phi_3(x,y)=y$ .
  - rectangle reference element □[(0,0),(1,0),(1,1),(0,1)]
    - node (0,0): $\phi_1(x,y)=(1-x)(1-y)$ ;
    - node (1,0): $\phi_2(x,y)=x(1-y)$ ;
    - node (1,1): $\phi_3(x,y)=xy$ ;
    - node (0,1): $\phi_4(x,y)=(1-x)y$ .



### Numerical integral

Here, we mainly use Gauss integral to calculate the integral. ==pay attention to the affine transformation==

# Numerical result

## 2D Rectangle

- Model problem:
  $$
  \begin{aligned}
  -\nabla \cdot(\nabla u)&=f, \quad\text { in } (0,1)^2 \\
  u&=0, \quad \text { on } \partial (0,1)^2 .
  \end{aligned}
  $$
  
- True solution: $u(x,y)=\sin(\pi x)\sin(\pi y)$ and $f=2\pi^2 \sin(\pi x)\sin(\pi y)$.

- NURBS information:

```matlab
nurbs_info.knotU = [0,0,1,1]; % row vector 1*m
nurbs_info.knotV = [0,0,1,1];
nurbs_info.pu = 1;
nurbs_info.pv = 1;
nurbs_info.ContrPoints = [0,0;1,0;0,1;1,1]; % control points n*dim
nurbs_info.weights = [1,1,1,1]'; %column vector n*1
```



- Mesh: Three h-refinement

![Rectangle3Refinement](C:\IGA\IGA_Poisson\figure\Rectangle3Refinement.png)

- error estimate:

  - Octave

  | n(h-refinement)   | 1       | 2       | 3       | 4       | 5       | 6       |
  | ----------------- | ------- | ------- | ------- | ------- | ------- | ------- |
  | $L^2$-norm        | 2.60e-2 | 2.03e-3 | 2.18e-4 | 2.61e-5 | 3.23e-6 | 4.03e-7 |
  | convergence order |         | 3.6769  | 3.2211  | 3.0611  | 3.0157  | 3.0040  |
  | $H^1$-seminorm    | 2.79e-1 | 5.53e-2 | 1.30e-2 | 3.21e-3 | 7.99e-4 | 2.00e-4 |
  | convergence order |         | 2.3364  | 2.0858  | 2.0215  | 2.0054  | 2.0013  |



- Model problem:

$$
\begin{aligned}
-\nabla \cdot(p\nabla u)+\mathbf{q}\cdot\nabla u+ru&=f, \quad\text { in } (0,1)^2 \\
u&=0, \quad \text { on } \partial (0,1)^2 .
\end{aligned}
$$



## 2D Circle

- Model problem:
  $$
  \begin{aligned}
  -\nabla \cdot(\nabla u)&=f, \quad\text { in } x^2+y^2\leqslant 1 \\
  u&=0, \quad \text { on } x^2+y^2=1 .
  \end{aligned}
  $$
  
- True solution: $u(x,y)=\sin(\frac14-x^2-y^2)$ and $f=4\cos(\frac14-x^2-y^2)+4(x^2+y^2)\sin(\frac14-x^2-y^2)$.

- NURBS information:

```matlab
nurbs_info.knotU = [0,0,0,1,1,1];
nurbs_info.knotV = [0,0,0,1,1,1];
nurbs_info.pu = 2;
nurbs_info.pv = 2;
a = sqrt(2)/4;
nurbs_info.ContrPoints = [-a,a;-2*a,0;-a,-a;0,2*a;0,0;0,-2*a;a,a;2*a,0;a,-a];
nurbs_info.weights = [1,2*a,1,2*a,1,2*a,1,2*a,1]';
```



- Mesh: Three times h-refinement

![Disk3Refinement](C:\IGA\IGA_Poisson\figure\Disk3Refinement.png)

- error estimate:

  - Octave(Improved IGA)

  | n(h-refinement)   |    1    |    2    |    3    |    4    |    5    |    6     |
  | ----------------- | :-----: | :-----: | :-----: | :-----: | :-----: | :------: |
  | $L^2$-norm        | 1.74e-4 | 2.42e-5 | 1.27e-6 | 7.90e-8 | 5.00e-9 | 3.16e-10 |
  | convergence order |         |  2.85   |  4.26   |  4.01   |  3.98   |   3.99   |
  | $H^1$-seminorm    | 3.53e-3 | 6.72e-4 | 7.28e-5 | 9.13e-6 | 1.16e-6 | 1.46e-7  |
  | convergence order |         |  2.39   |  3.21   |  3.00   |  2.98   |   2.99   |
  
  - Octave(IGA)
  
  | n(h-refinement)   |    1    |    2    |    3    |    4    |    5    |    6     |
  | ----------------- | :-----: | :-----: | :-----: | :-----: | :-----: | :------: |
  | $L^2$-norm        | 8.85e-5 | 1.58e-5 | 8.29e-7 | 5.34e-8 | 3.46e-9 | 2.21e-10 |
  | convergence order |         |  2.49   |  4.25   |  3.96   |  3.95   |   3.97   |
  | $H^1$-seminorm    | 1.86e-3 | 4.20e-4 | 4.73e-5 | 6.09e-6 | 7.86e-7 | 1.00e-7  |
  | convergence order |         |  2.15   |  3.15   |  2.96   |  2.95   |   2.97   |



## 3D Cubic

- Model problem:
  $$
  \begin{aligned}
  -\nabla \cdot(\nabla u)&=f, \quad\text { in } (-1,1)^3 \\
  u&=0, \quad \text { on } \partial (-1,1)^3 .
  \end{aligned}
  $$

- True solution: $u(x,y,z)=\sin(\pi x)\sin(\pi y)\sin(\pi z)$  and $f=3\pi^2 \sin(\pi x)\sin(\pi y)\sin(\pi z)$ .

- NURBS information

```matlab
nurbs_info.knotU = [0,0,1,1]; % row vector 1*m
nurbs_info.knotV = [0,0,1,1];
nurbs_info.knotW = [0,0,1,1];
nurbs_info.pu = 1;
nurbs_info.pv = 1;
nurbs_info.pw = 1;
a = -1; b = 1;% \Omega = [a,b]^3
nurbs_info.ContrPoints = zeros(8,3); % control points n*dim
nurbs_info.ContrPoints(:,1) = repmat([a;b],4,1);% x-axis
nurbs_info.ContrPoints(:,2) = repmat([a;a;b;b],2,1);% y-axis
nurbs_info.ContrPoints(:,3) = [a;a;a;a;b;b;b;b];% z-axis
% nurbs_info.ContrPoints = zeros(2,2,2);
% nurbs_info.ContrPoints(:,:,1) = [0,0;1,1];
% nurbs_info.ContrPoints(:,:,2) = [0,1;0,1];
nurbs_info.weights = [1,1,1,1,1,1,1,1]'; %column vector n*1
```



- error estimate:

  - Matlab

  | n(h-refinement)   |    1    |    2    |    3    |    4    |    5    |    6    |
  | ----------------- | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: |
  | $L^2$-norm        | 9.17e-2 | 6.30e-2 | 4.96e-3 | 5.34e-4 | 6.40e-5 | 7.91e-6 |
  | convergence order |         | 0.5421  | 3.6653  | 3.2172  | 3.06000 | 3.0154  |
  | $H^1$-seminorm    | 4.02e-1 | 7.10e-1 | 1.36e-1 | 3.19e-2 | 7.86e-3 | 1.96e-3 |
  | convergence order |         | -0.8199 | 2.3818  | 2.0931  | 2.0230  | 2.0057  |




## 3D Cyclinder

- Model problem:
  $$
  \begin{aligned}
  -\nabla \cdot(\nabla u)&=f, \quad\text { in } (-1,1)^3 \\
  u&=0, \quad \text { on } \partial (-1,1)^3 .
  \end{aligned}
  $$

- True solution: $u(x,y,z)=\sin(x^2+y^2-1)z(1-z)$ and 
  $$
  f=4(x^2+y^2)\sin(x^2+y^2-1)z(1-z)-4\cos(x^2+y^2-1)z(1-z)+2\sin(x^2+y^2-1).
  $$
  
- NURBS information:

```matlab
nurbs_info.knotU = [0,0,0,1,1,1];
nurbs_info.knotV = [0,0,0,1,1,1];
nurbs_info.knotW = [0,0,0,1,1,1];
nurbs_info.pu = 2;
nurbs_info.pv = 2;
nurbs_info.pw = 2;
a = sqrt(2)/2;

nu = length(nurbs_info.knotU) - nurbs_info.pu - 1;
nv = length(nurbs_info.knotV) - nurbs_info.pv - 1;
nw = length(nurbs_info.knotW) - nurbs_info.pw - 1;
ConPts=zeros(nu,nv,nw,3);
weights=zeros(nu,nv,nw);

X=[-a,0,a; -2*a   0   2*a;-a, 0, a];
Y=[ a,2*a,  a; 0,0, 0;-a, -2*a, -a];

w=[1,a,1;a,1,a;1,a,1];


for i=1:nu
    for j=1:nv
        ConPts(i,j,1,:)=[X(i,j),Y(i,j),0];
        weights(i,j,1)=w(i,j);
    end
end

for i=1:nu
    for j=1:nv
        ConPts(i,j,2,:)=[X(i,j),Y(i,j),0.5];
        weights(i,j,2)=w(i,j);
    end
end

for i=1:nu
    for j=1:nv
        ConPts(i,j,3,:)=[X(i,j),Y(i,j),1];
        weights(i,j,3)=w(i,j);
    end
end
nurbs_info.ContrPoints = reshape(ConPts,nu*nv*nw,3);
nurbs_info.weights = reshape(weights,nu*nv*nw,1);
```



- error estimate:

  - Matlab

  | n(h-refinement)   |    1    |    2    |    3    |    4    |    5    |    6    |
  | ----------------- | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: |
  | $L^2$-norm        | 6.47e-3 | 9.61e-4 | 9.60e-5 | 9.21e-6 | 1.02e-6 | 1.23e-7 |
  | convergence order |         | 2.7519  | 3.3233  | 3.3814  | 3.1702  | 3.0523  |
  | $H^1$-seminorm    | 4.44e-2 | 9.12e-3 | 2.33e-3 | 5.48e-4 | 1.33e-4 | 3.31e-5 |
  | convergence order |         | 2.2854  | 1.9646  | 2.0914  | 2.0392  | 2.0120  |
  
  - Octave
  
  | n(h-refinement)   |    1    |    2    |    3    |    4    |    5    |    6    |
  | ----------------- | :-----: | :-----: | :-----: | :-----: | :-----: | :-----: |
  | $L^2$-norm        | 6.99e-3 | 9.33e-4 | 9.47e-5 | 9.06e-6 | 1.00e-6 | 1.21e-7 |
  | convergence order |         | 2.0961  | 3.3001  | 3.3861  | 3.1755  | 3.0544  |
  | $H^1$-seminorm    | 4.71e-2 | 8.58e-3 | 2.24e-3 | 5.22e-4 | 1.27e-4 | 3.14e-5 |
  | convergence order |         | 2.4549  | 1.9396  | 2.0995  | 2.0436  | 2.0134  |



- Future improvement:
  - non-homogeneous boundary conditions;
  - $ru$ and $\mathbf{q}\cdot\nabla u$ in the model equation;
  - complex domain



- Analysis:
  - BSplines/NURBS is not interpolated at the interior knots.
  - not easy to deal with the non-homogeneous boundary conditions.(u=g on the $\partial\Omega$)
