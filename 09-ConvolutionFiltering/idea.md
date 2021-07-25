# Adaptivity + Convolution Filtering

1. Obtain $u_h$ on an adaptive mesh $T_h$ (__how and when to stop?__)
2. Obtain the background uniform mesh $T_{h_0}$ and Gauss quadrature points $\{\xi_n\}$ on $(-1,1)$. The Gauss quad points on each eleement $K_l$ are $\{\frac{h}{2}\xi_h+x_{mid}^l\}$
3. Figure out the relation of $T_h$ and $T_{h_0}$. 
    * Need to know what are the elements $K \in T_h$ that are contained in each $ K_0 \in T_{h_0}$
    * Next, which element $K \in T_h$ does each Gauss quad point locate?
4. Evaluate $u_h$ on all the Gauss quad points  $\{(u_h^l)_n\}$ on each element $K_l$. 
5. **Precompute** 
$$A_{n,s,r,j} = \psi(\frac{\xi_n}{2} - \frac{z_s}{2} - j - r)$$
where $\xi_n$ and $z_s$ are the GQ pts in step 2. We have $1\leq n,s \leq N_{GQ}$, $-2k\leq j \leq 2k $ and $-k\leq r\leq k$

5. Convolution at point $ x = \frac{h}{2} \xi_n + x_{mid}^m$ on element $K_m$: 
$$
 \begin{align*}
 K_h*u_h &= \frac{1}{h}\sum_{r}\sum_{l}C_r\int_{x_l}^{x_r}\psi(\frac{x-y}{h} - r )\ u_h^l(y) dy\\
 & = \frac{1}{2}\sum_{r}\sum_{l}C_r\int_{-1}^{1}\psi(\frac{x-\frac{h}{2} z - x_{mid}^l}{h} - r )\ u_h(\frac{h}{2} z + x_{mid}^l)  dz\\
 & = \frac{1}{2}\sum_{r}\sum_{l}C_r\int_{-1}^{1}\psi(\frac{\xi_n}{2} - \frac{z}{2} - (l-m) - r) \ u_h(\frac{h}{2}z+x_{mid}^l)dz\\
 & = \text{...Gauss quadrature of z ... and }j = l-m \\
 & = \frac{1}{2}\sum_{r}\sum_{j}\sum_s C_rA_{n,,s,r,j}\ (u_h^l)_s \ w_s\\
 & = \sum_{l}\sum_{s} M_{n,s,j}(u_h^l)_s w_s
 \end{align*}
 $$
 where $$M_{n,s,j} = \frac{1}{2}\sum_{r}C_r A_{n,s,r,j} $$
 Note $-2k \leq j \leq 2k$. This means $$ m-2k\leq l\leq m+2k $$
 In other words, to get $K_h*u_h \ (x)$ on element $K_m$, we need to use info from the its $2k$ neighbors on each direction

 

## Steps we done
1. figure out the relation of two meshes 
2. ~~Obtain uh at GQ points of coarse mesh~~
3. ~~check GQ points are correct~~
4. precompute convolution matrix
5. do convolution 
2. ~~Check simple smooth case if convolution~~


## Steps to do
1. L2 projection of uh from finer mesh to coarse


 