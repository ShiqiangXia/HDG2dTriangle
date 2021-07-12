# Adaptivity + Convolution Filtering

1. Obtain $u_h$ on an adaptive mesh $T_h$
2. Determin a background uniform mesh $T_{h_0}$ and Gauss quadrature points $\{\xi_n\}$
3. Figure out the relation of $T_h$ and $T_{h_0}$. 
    * Need to know what are the elements $K \in T_h$ that are contained in each $ K_0 \in T_{h_0}$
    * Next, which element $K \in T_h$ does each Gauss quad point localte?
4. Evaluate $u_h$ on all the Gauss quad points $\{(u_h^l)_n\}$
5. Precompute 
$$A_{n,j,s,r} = \psi(\frac{\xi_n}{2} - \frac{z_s}{2} - j - r)$$
where $1\leq n,s \leq N_{GQ}$, $-2k\leq j \leq 2k $ and $-k\leq r\leq r$

5. Convolution at $ x = \frac{h}{2} \xi_n + x_{mid}^m$ : 
$$
 \begin{align*}
 K_h*u_h &= \frac{1}{h}\sum_{r}\sum_{l}\int_{x_l}^{x_r}\psi(\frac{x-y}{h} - r )\ u_h^l(y) dy\\
 & = \frac{1}{2}\sum_{r}\sum_{l}\int_{-1}^{1}\psi(\frac{x-\frac{h}{2} z - x_{mid}^l}{h} - r )\ u_h(\frac{h}{2} z + x_{mid}^l)  dz\\
 & = \frac{1}{2}\sum_{r}\sum_{l}\int_{-1}^{1}\psi(\frac{\xi_n}{2} - \frac{z}{2} - (l-m) - r) \ u_h(\frac{h}{2}z+x_{mid}^l)dz\\
 & = \frac{1}{2}\sum_{r}\sum_{l}\sum_s A_{n,(l-m),s,r}(u_h^l)_s w_s
 \end{align*}
 $$

 Take $z$ be the Gauss quad pts, then $\frac{h}{2}z+x_{mid}^l$ are the Gauss quads on element $l$. Then $u_h(\frac{h}{2}z+x_{mid}^l)$ are the values computed in stpe 4

 