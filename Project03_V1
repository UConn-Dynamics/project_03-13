### A Pluto.jl notebook ###    
      2  # v0.20.21                                                           
      3          
      4 +#> risky_file_source = "https://github.com/UConn-Dynamics/projec
        +t_03-13/blob/main/project_03.jl?raw=true"                              
      5 +                                                                
      6  using Markdown                                                         
      7  using InteractiveUtils                                               
      8  
     ...
       12  ![Sliding compound pendulum with a support block connected to
           a spring and rotating compound pendulum](https://raw.githubuse
           rcontent.com/cooperrc/me5180-project_02/refs/heads/main/spring
           _compound-2_bodies.png)
       13  
       14  In this project, a rigid bar is connected to a sliding block a
           long a
       13 -horizontal tracks. The sliding block is connected to a spring 
          -that stretches and compresses. The rigid bar $L = 0.4~m$ acts 
          -as a compound pendulum.                                       
       15 +horizontal tracks. The sliding block is connected to a spring 
          +that stretches and compresses. The rigid bar $L = 0.4~m$ acts 
          +as a compound pendulum.                                       
       16  
       17  1. $x_1-y_1-$ describes block 1 position and orientation, $\th
           eta_1$
       18  2. $x_2-y_2-$ describes the rigid bar position and orientation
           , $\theta_2$
       19  
       18 -The applied forces are,                                       
       20 +The applied forces are,                                       
       21  
       22  1. Spring attached to block 1, $F = -k x_1$ where $k = 10~N/m$
       23  2. gravity acting on block 1 and the rigid bar, $F_1 = -m_1g\h
           at{j}$ and $F_2 = -m_2 g\hat{j}$ where $m_1 = 0.1~kg$ and $m_2
            = 0.3~kg$
       22 -                                                              
       23 -In this project, you need to                                  
       24  
       25 +In this project, you need to                                  
       26 +                                                              
       27  1. determine constraint equations $C(\mathbf{q},~t)$
       28  2. Create an augmented solution method for the dynamic motion
           of these two moving parts
       29  3. visualize the motion of the system as the two parts complet
           e at least one oscillation
       30  4. calculate and show (graph or vectors) the constraint forces
            acting on the 2-body system
       31  """
       32  
       33 +# ╔═╡ 1a2b3c4d-06c0-11f1-a2b0-79e68ed152eb                    
       34 +md"""## 1. Constraint Equations $C(\mathbf{q},\,t)$           
       35 +                                                              
       36 +The **generalized coordinate vector** is:                     
       37 +                                                              
       38 +$$\mathbf{q} = [x_1,\ y_1,\ \theta_1,\ x_2,\ y_2,\ \theta_2]^T
          +$$                                                            
       39 +                                                              
       40 +Body 1 is the **sliding block** (3 DOF), body 2 is the **rigid
          + bar** (3 DOF), giving 6 total coordinates.                   
       41 +Four holonomic constraints reduce the system to **2 independen
          +t DOF** ($x_1$ and $\theta_2$).                               
       42 +                                                              
       43 +| # | $C_i(\mathbf{q})$ | Physical meaning |                  
       44 +|---|---|---|                                                 
       45 +| $C_1$ | $y_1 = 0$ | Block constrained to horizontal track | 
       46 +| $C_2$ | $\theta_1 = 0$ | Block does not rotate |            
       47 +| $C_3$ | $x_1 - x_2 + \tfrac{L}{2}\sin\theta_2 = 0$ | Pin joi
          +nt — x direction |                                            
       48 +| $C_4$ | $y_2 + \tfrac{L}{2}\cos\theta_2 = 0$ | Pin joint — y
          + direction (using $y_1=0$) |                                  
       49 +                                                              
       50 +$C_3$ and $C_4$ enforce that the **top of the bar** (body 2 lo
          +cal offset $[0,\,L/2]$) coincides with the **block center** (p
          +in point).                                                    
       51 +                                                              
       52 +The **constraint Jacobian** $\mathbf{C_q} = \partial\mathbf{C}
          +/\partial\mathbf{q}$:                                         
       53 +                                                              
       54 +$$\mathbf{C_q} = \begin{bmatrix}                              
       55 +0 & 1 & 0 & 0 & 0 & 0 \\                                      
       56 +0 & 0 & 1 & 0 & 0 & 0 \\                                      
       57 +1 & 0 & 0 & -1 & 0 & \tfrac{L}{2}\cos\theta_2 \\              
       58 +0 & 0 & 0 & 0 & 1 & -\tfrac{L}{2}\sin\theta_2                 
       59 +\end{bmatrix}$$                                               
       60 +                                                              
       61 +Differentiating constraints twice gives the **acceleration RHS
          +** $\boldsymbol{\gamma} = -\dot{\mathbf{C}}_q\dot{\mathbf{q}}$
          +:                                                             
       62 +                                                              
       63 +$$\boldsymbol{\gamma} = \begin{bmatrix} 0 \\ 0 \\ \tfrac{L}{2}
          +\sin\theta_2\,\dot\theta_2^2 \\ \tfrac{L}{2}\cos\theta_2\,\dot
          +\theta_2^2 \end{bmatrix}$$                                    
       64 +"""                                                           
       65 +                                                              
       66 +# ╔═╡ 2b3c4d5e-06c0-11f1-a2b0-79e68ed152eb                    
       67 +md"""## 2. Augmented Solution Method                          
       68 +                                                              
       69 +The **augmented Lagrangian equations of motion** with Lagrange
          + multipliers $\boldsymbol{\lambda}$:                          
       70 +                                                              
       71 +$$\begin{bmatrix}\mathbf{M} & \mathbf{C_q}^T \\ \mathbf{C_q} &
          + \mathbf{0}\end{bmatrix} \begin{bmatrix}\ddot{\mathbf{q}} \\ \
          +boldsymbol{\lambda}\end{bmatrix} = \begin{bmatrix}\mathbf{Q} \
          +\ \boldsymbol{\gamma}_{\mathrm{stab}}\end{bmatrix}$$          
       72 +                                                              
       73 +where:                                                        
       74 +                                                              
       75 +- $\mathbf{M} = \mathrm{diag}(m_1,\ m_1,\ I_1,\ m_2,\ m_2,\ I_
          +2)$ — mass/inertia matrix                                     
       76 +- $\mathbf{Q} = [-kx_1,\ -m_1g,\ 0,\ 0,\ -m_2g,\ 0]^T$ — gener
          +alized applied forces                                         
       77 +- $I_2 = m_2 L^2/12$ — bar moment of inertia about its center 
          +of mass                                                       
       78 +                                                              
       79 +**Baumgarte stabilization** prevents constraint drift in the n
          +umerical integration:                                         
       80 +                                                              
       81 +$$\boldsymbol{\gamma}_{\mathrm{stab}} = \boldsymbol{\gamma} - 
          +2\alpha\,\mathbf{C_q}\dot{\mathbf{q}} - \beta^2\,\mathbf{C}$$ 
       82 +                                                              
       83 +with $\alpha = \beta = 5$.  The augmented $(10\times10)$ linea
          +r system is solved at each ODE step for $[\ddot{\mathbf{q}};\ 
          +\boldsymbol{\lambda}]$.                                       
       84 +                                                              
       85 +**Physical interpretation of $\boldsymbol{\lambda}$:**        
       86 +                                                              
       87 +| Multiplier | Force |                                        
       88 +|---|---|                                                     
       89 +| $\lambda_1$ | Track normal force on block |                 
       90 +| $\lambda_2$ | Moment preventing block rotation |            
       91 +| $\lambda_3$ | Pin joint force — x direction |               
       92 +| $\lambda_4$ | Pin joint force — y direction |               
       93 +"""                                                           
       94 +                                                              
       95 +# ╔═╡ 3c4d5e6f-06c0-11f1-a2b0-79e68ed152eb                    
       96 +begin                                                         
       97 +  using LinearAlgebra                                         
       98 +  using DifferentialEquations                                 
       99 +  using Plots                                                 
      100 +end                                                           
      101 +                                                              
      102 +# ╔═╡ 4d5e6f70-06c0-11f1-a2b0-79e68ed152eb                    
      103 +begin                                                         
      104 +  # Physical parameters                                       
      105 +  L_bar = 0.4     # bar length [m]                            
      106 +  k_spr = 10.0    # spring stiffness [N/m]                    
      107 +  m1    = 0.1     # block mass [kg]                           
      108 +  m2    = 0.3     # bar mass [kg]                             
      109 +  g_acc = 9.81    # gravitational acceleration [m/s²]         
      110 +  I1    = 1e-4    # block MOI [kg⋅m²] (small; θ₁ is constraine
          +d)                                                            
      111 +  I2    = m2 * L_bar^2 / 12   # bar MOI about center [kg⋅m²]  
      112 +                                                              
      113 +  M_sys = Diagonal([m1, m1, I1, m2, m2, I2])                  
      114 +                                                              
      115 +  md"""**Parameters loaded:**                                 
      116 +  - $L = $(L_bar)$ m, $k = $(k_spr)$ N/m                      
      117 +  - $m_1 = $(m1)$ kg, $m_2 = $(m2)$ kg                        
      118 +  - $I_2 = $(round(I2, digits=5))$ kg⋅m²"""                   
      119 +end                                                           
      120 +                                                              
      121 +# ╔═╡ 5e6f7071-06c0-11f1-a2b0-79e68ed152eb                    
      122 +begin                                                         
      123 +  # Constraint vector C(q)                                    
      124 +  function C_vec(q, L)                                        
      125 +    x1, y1, θ1, x2, y2, θ2 = q                                
      126 +    return [y1,                                               
      127 +            θ1,                                               
      128 +            x1 - x2 + (L/2)*sin(θ2),                          
      129 +            y2 + (L/2)*cos(θ2)]                               
      130 +  end                                                         
      131 +                                                              
      132 +  # Constraint Jacobian ∂C/∂q                                 
      133 +  function Cq_jac(q, L)                                       
      134 +    θ2 = q[6]                                                 
      135 +    return [0  1  0   0  0   0;                               
      136 +            0  0  1   0  0   0;                               
      137 +            1  0  0  -1  0   (L/2)*cos(θ2);                   
      138 +            0  0  0   0  1  -(L/2)*sin(θ2)]                   
      139 +  end                                                         
      140 +                                                              
      141 +  # Acceleration RHS  γ = -Ċq⋅q̇                               
      142 +  function γ_rhs(q, qdot, L)                                  
      143 +    θ2    = q[6]                                              
      144 +    θ2dot = qdot[6]                                           
      145 +    return [0.0,                                              
      146 +            0.0,                                              
      147 +            (L/2)*sin(θ2)*θ2dot^2,                            
      148 +            (L/2)*cos(θ2)*θ2dot^2]                            
      149 +  end                                                         
      150 +                                                              
      151 +  "Constraint functions defined ✓"                            
      152 +end                                                           
      153 +                                                              
      154 +# ╔═╡ 6f708182-06c0-11f1-a2b0-79e68ed152eb                    
      155 +begin                                                         
      156 +  # Solve augmented system for q̈ and λ at a given state       
      157 +  function solve_augmented(q, qdot, M, L, k, m1, m2, g; α=5.0,
          + β=5.0)                                                       
      158 +    Cq = Cq_jac(q, L)                                         
      159 +    Q  = [-k*q[1], -m1*g, 0.0, 0.0, -m2*g, 0.0]               
      160 +    γ  = γ_rhs(q, qdot, L)                                    
      161 +    C  = C_vec(q, L)                                          
      162 +    Cd = Cq * qdot                                            
      163 +                                                              
      164 +    # Baumgarte stabilization                                 
      165 +    γ_s = γ .- 2α .* Cd .- β^2 .* C                           
      166 +                                                              
      167 +    # Build and solve 10×10 augmented system                  
      168 +    n, nc = 6, 4                                              
      169 +    A = [Matrix(M) Cq'; Cq zeros(nc, nc)]                     
      170 +    b = [Q; γ_s]                                              
      171 +    x = A \ b                                                 
      172 +    return x[1:n], x[n+1:end]   # q̈, λ                        
      173 +  end                                                         
      174 +                                                              
      175 +  "Augmented solver defined ✓"                                
      176 +end                                                           
      177 +                                                              
      178 +# ╔═╡ 70819293-06c0-11f1-a2b0-79e68ed152eb                    
      179 +begin                                                         
      180 +  function dynamics!(dz, z, p, t)                             
      181 +    q    = z[1:6]                                             
      182 +    qdot = z[7:12]                                            
      183 +    M, L, k, m1, m2, g = p                                    
      184 +    qddot, _ = solve_augmented(q, qdot, M, L, k, m1, m2, g)   
      185 +    dz[1:6]  = qdot                                           
      186 +    dz[7:12] = qddot                                          
      187 +  end                                                         
      188 +                                                              
      189 +  "ODE function defined ✓"                                    
      190 +end                                                           
      191 +                                                              
      192 +# ╔═╡ 8192a3b4-06c0-11f1-a2b0-79e68ed152eb                    
      193 +begin                                                         
      194 +  # Initial conditions: bar tilted 30° from vertical, system a
          +t rest                                                        
      195 +  θ2_0 = π/6                                                  
      196 +  q0   = [0.0,                                                
      197 +          0.0,                                                
      198 +          0.0,                                                
      199 +          (L_bar/2)*sin(θ2_0),                                
      200 +         -(L_bar/2)*cos(θ2_0),                                
      201 +          θ2_0]                                               
      202 +  z0    = [q0; zeros(6)]                                      
      203 +  tspan = (0.0, 8.0)                                          
      204 +  p_ode = (M_sys, L_bar, k_spr, m1, m2, g_acc)                
      205 +                                                              
      206 +  prob = ODEProblem(dynamics!, z0, tspan, p_ode)              
      207 +  sol  = solve(prob, Rodas4(), reltol=1e-8, abstol=1e-8)      
      208 +                                                              
      209 +  md"**Simulation complete:** $(length(sol.t)) steps over $(ts
          +pan[2]) s"                                                    
      210 +end                                                           
      211 +                                                              
      212 +# ╔═╡ 92a3b4c5-06c0-11f1-a2b0-79e68ed152eb                    
      213 +begin                                                         
      214 +  t_arr  = sol.t                                              
      215 +  x1_arr = [u[1] for u in sol.u]                              
      216 +  θ2_arr = [u[6] for u in sol.u]                              
      217 +                                                              
      218 +  p_x = plot(t_arr, x1_arr,                                   
      219 +             label="x₁ — block position", lw=2, color=:royalbl
          +ue,                                                           
      220 +             xlabel="Time [s]", ylabel="Position [m]",        
      221 +             title="3. Motion of the 2-Body System")          
      222 +  p_θ = plot(t_arr, rad2deg.(θ2_arr),                         
      223 +             label="θ₂ — bar angle", lw=2, color=:crimson,    
      224 +             xlabel="Time [s]", ylabel="Angle [°]")           
      225 +                                                              
      226 +  plot(p_x, p_θ, layout=(2,1), size=(700, 480), legend=:toprig
          +ht)                                                           
      227 +end                                                           
      228 +                                                              
      229 +# ╔═╡ a3b4c5d6-06c0-11f1-a2b0-79e68ed152eb                    
      230 +begin                                                         
      231 +  step_a = max(1, length(sol.t) ÷ 200)                        
      232 +  anim   = @animate for i in 1:step_a:length(sol.t)           
      233 +    u = sol.u[i]                                              
      234 +    x1, _, _, x2, y2, θ2 = u[1:6]                             
      235 +                                                              
      236 +    # Bar endpoints in global frame                           
      237 +    x_top = x2 - (L_bar/2)*sin(θ2)                            
      238 +    y_top = y2 + (L_bar/2)*cos(θ2)                            
      239 +    x_bot = x2 + (L_bar/2)*sin(θ2)                            
      240 +    y_bot = y2 - (L_bar/2)*cos(θ2)                            
      241 +                                                              
      242 +    # Block outline                                           
      243 +    bw = 0.04                                                 
      244 +    bx = [x1-bw, x1+bw, x1+bw, x1-bw, x1-bw]                  
      245 +    by = [-bw,   -bw,    bw,     bw,    -bw]                  
      246 +                                                              
      247 +    plot(xlim=(-0.35, 0.35), ylim=(-0.45, 0.12),              
      248 +         aspect_ratio=:equal, xlabel="x [m]", ylabel="y [m]", 
      249 +         title="t = $(round(sol.t[i], digits=2)) s",          
      250 +         size=(520, 420), legend=:topright)                   
      251 +    hline!([0.0], lw=1, color=:black, linestyle=:dash, label="
          +Track")                                                       
      252 +    plot!(bx, by, seriestype=:shape, color=:gray,  label="Bloc
          +k")                                                           
      253 +    plot!([x_top, x_bot], [y_top, y_bot], lw=5, color=:steelbl
          +ue, label="Bar")                                              
      254 +    scatter!([x_top], [y_top], ms=7, color=:gray, label="Pin")
      255 +    scatter!([x2],    [y2],    ms=5, color=:red,  label="CoM₂"
          +)                                                             
      256 +  end                                                         
      257 +  gif(anim, fps=30)                                           
      258 +end                                                           
      259 +                                                              
      260 +# ╔═╡ b4c5d6e7-06c0-11f1-a2b0-79e68ed152eb                    
      261 +begin                                                         
      262 +  # Recover Lagrange multipliers (constraint forces) at every 
          +time step                                                     
      263 +  λ_all = [solve_augmented(sol.u[i][1:6], sol.u[i][7:12],     
      264 +                            M_sys, L_bar, k_spr, m1, m2, g_acc
          +)[2]                                                          
      265 +           for i in eachindex(sol.t)]                         
      266 +                                                              
      267 +  λ1_arr = [l[1] for l in λ_all]   # track normal force [N]   
      268 +  λ2_arr = [l[2] for l in λ_all]   # moment preventing block r
          +otation [N⋅m]                                                 
      269 +  λ3_arr = [l[3] for l in λ_all]   # pin force — x [N]        
      270 +  λ4_arr = [l[4] for l in λ_all]   # pin force — y [N]        
      271 +                                                              
      272 +  p_λ1 = plot(t_arr, λ1_arr, lw=2, color=:seagreen,           
      273 +              label="λ₁ — track normal [N]",                  
      274 +              title="4. Constraint Forces on 2-Body System",  
      275 +              ylabel="Force [N]")                             
      276 +  p_λ3 = plot(t_arr, λ3_arr, lw=2, color=:darkorange,         
      277 +              label="λ₃ — pin force x [N]", ylabel="Force [N]"
          +)                                                             
      278 +  p_λ4 = plot(t_arr, λ4_arr, lw=2, color=:mediumpurple,       
      279 +              label="λ₄ — pin force y [N]",                   
      280 +              xlabel="Time [s]", ylabel="Force [N]")          
      281 +                                                              
      282 +  plot(p_λ1, p_λ3, p_λ4, layout=(3,1), size=(700, 600), legend
          +=:topright)                                                   
      283 +end                                                           
      284 +                                                              
      285  # ╔═╡ 0d9be664-d7c5-4084-add2-25e5418742d6
      286  
      287  
      288  # ╔═╡ 00000000-0000-0000-0000-000000000001
      289  PLUTO_PROJECT_TOML_CONTENTS = """
      290  [deps]
      291 +DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
      292 +LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"        
      293 +Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"                
      294  """
      295  
      296  # ╔═╡ 00000000-0000-0000-0000-000000000002
     ...
      306  
      307  # ╔═╡ Cell order:
      308  # ╟─f17103ea-06bf-11f1-a2b0-79e68ed152eb
      309 +# ╟─1a2b3c4d-06c0-11f1-a2b0-79e68ed152eb                      
      310 +# ╟─2b3c4d5e-06c0-11f1-a2b0-79e68ed152eb                      
      311 +# ╠═3c4d5e6f-06c0-11f1-a2b0-79e68ed152eb                      
      312 +# ╠═4d5e6f70-06c0-11f1-a2b0-79e68ed152eb                      
      313 +# ╠═5e6f7071-06c0-11f1-a2b0-79e68ed152eb                      
      314 +# ╠═6f708182-06c0-11f1-a2b0-79e68ed152eb                      
      315 +# ╠═70819293-06c0-11f1-a2b0-79e68ed152eb                      
      316 +# ╠═8192a3b4-06c0-11f1-a2b0-79e68ed152eb                      
      317 +# ╠═92a3b4c5-06c0-11f1-a2b0-79e68ed152eb                      
      318 +# ╠═a3b4c5d6-06c0-11f1-a2b0-79e68ed152eb                      
      319 +# ╠═b4c5d6e7-06c0-11f1-a2b0-79e68ed152eb                      
      320  # ╠═0d9be664-d7c5-4084-add2-25e5418742d6
      321  # ╟─00000000-0000-0000-0000-000000000001
      322  # ╟─00000000-0000-0000-0000-000000000002
