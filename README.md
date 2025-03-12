\begin{frame}[fragile]
    \frametitle{Example 1: 二维泊松方程}
    \begin{minted}{julia}
        # FE space
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V = TestFESpace(model, reffe,
        dirichlet_tags="fixed")
        
        g(x) = 0.0
        U = TrialFESpace(V, g)
        
        f(x) = 1.0
        a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
        l(v) = ∫( v*f )dΩ
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 1: 二维泊松方程}
    \begin{minted}{julia}
        op = AffineFEOperator(a, l, U, V)
        
        uh = solve(op)
        
        writevtk(Ω, "results", cellfields=["uh"=>uh])
    \end{minted}
\end{frame}

\begin{frame}{The Nitsche Method}
    We introduce a method for treating general boundary conditions in the finite element method generalizing an approach, due to Nitsche, for approximating Dirichlet boundary conditions. \\
    Find $u_h \in V_h$ such that
    \begin{equation*}
        \mathcal{A}_h(u_h,v)=\mathcal{l}_h(v).
    \end{equation*}
    for all $v \in V_h$. \\
\end{frame}

\begin{frame}{The Nitsche Method}
    \begin{equation*}
        \mathcal{A}_h(u_h,v) = \int_{\Omega}{\nabla u \cdot \nabla v}\mathrm{d}\Omega + \int_{\Gamma_D}{\left( \lambda uv-v(\mathbf{n}\cdot \nabla u) -u(\mathbf{n}\cdot \nabla v)\right) }\mathrm{d}s
    \end{equation*}
    \begin{equation*}
        \mathcal{l}_h(v)=\int_{\Omega}{vf}\mathrm{d}\Omega  + \int_{\Gamma_D}{\left( \lambda vg_D-\dfrac{\partial v}{\partial \mathbf{n}}g_D \right) }\mathrm{d}s
    \end{equation*}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 2: Nitsche法求解}
    \begin{minted}{julia}
        # 读取网格
        Ω = Triangulation(model)
        dΩ = Measure(Ω, 2)
        
        Γd = BoundaryTriangulation(model,tags="fixed")
        dΓd = Measure(Γd, 2)
        n_Γd = get_normal_vector(Γd)
        
        # FE space
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V = TestFESpace(model, reffe, conformity=:H1)
        
        U = TrialFESpace(V)
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 2: Nitsche法求解}
    \begin{minted}{julia}
        γd = 10.0
        h = 1.0 / 400
        
        ud(x) = 0.0
        f(x) = 1.0
        
        a(u,v) =∫( ∇(v)⋅∇(u) )dΩ +
        ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u )dΓd
        
        l(v) =  ∫( v*f )dΩ  +
        ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd
    \end{minted}
\end{frame}

\begin{frame}{Example 2: Nitsche法求解}
    分析结果对比：
    \begin{figure}
        \centering %居中
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/002.png}
            \caption{一般方法}
        \end{minipage}%
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/003.png}
            \caption{The Nitsche Method}
        \end{minipage}
    \end{figure}
\end{frame}

\section{Part II: CutFEM有限元}

\begin{frame}{CutFEM有限元}
    Like in any other \textcolor{red}{\emph{embedded boundary method}}, we build the computational mesh by introducing an \textbf{artificial
        domain} $\Omega_{\mathrm{art}}$ such that it has a simple geometry that is easy to mesh using Cartesian grids and it includes the \textbf{physical domain} $\Omega \subset \Omega_{\mathrm{art}}$.
    \begin{figure}
        \centering
        \includegraphics[height=0.5\textheight]{./img/004.png}
        \caption{Embedded boundary setup.}
        \label{f004}
    \end{figure}
\end{frame}

\begin{frame}{CutFEM有限元}
    Let us construct a partition of $\Omega_{\mathrm{art}}$ into cells, represented by $T^h_{\mathrm{art}}$ , with characteristic cell size $h$. Cells in $T^h_{\mathrm{art}}$ can be classified as
    follows: a cell $K \in T^h_{\mathrm{art}}$ such that $K \subset \Omega$ is an \textcolor{red}{\textbf{internal cell}}; if $K\cap \Omega =\varnothing$, $K$ is an \textcolor{green}{\textbf{external cell}}; otherwise,
    $K$ is a \textcolor{blue}{\textbf{cut cell}} (see Fig. \ref{f004}). Furthermore, we define the set of \textbf{active cells} as $T^\mathrm{act}_h = T^\mathrm{in}_h  \cup T^\mathrm{cut}_h $ and its union $\Omega_{\mathrm{act}}$.
    \begin{figure}
        \centering
        \includegraphics[height=0.4\textheight]{./img/005.png}
        \caption{Finite Element spaces.}
        \label{f005}
    \end{figure}
\end{frame}

\begin{frame}{CutFEM有限元}
    We define the FE-wise operators:
    
    \begin{gather*}%不会产生编号
        \mathcal{A}_K(u,v)\doteq\int_{K\cap \Omega }\nabla u \cdot \nabla v \mathrm{d}V \\
        +\int_{\partial K\cap \Gamma_D}\left(h^{-1}\lambda uv -v(\mathbf{n}\cdot \nabla u)-u(\mathbf{n}\cdot \nabla v) \right) \mathrm{d}S\\
        +\int_{\partial K\cap \Gamma_\mathrm{ghost}}h^{-1}\mu [\![\partial_nu]\!] [\![\partial_nv]\!] \mathrm{d}S
    \end{gather*}
    
    \begin{equation*}
        \mathcal{l}_K(u,v)\doteq \int_{K\cap \Omega }vf \mathrm{d}V + \int_{\partial K\cap \Gamma_D } \left(h^{-1}\lambda vg_D -(\mathbf{n} \cdot \nabla v)g_D)\right) \mathrm{d}S
    \end{equation*}
    
\end{frame}

\begin{frame}{CutFEM有限元}
    Given a face $F \in \mathcal{F}^{\mathrm{ghost}}_h$ and the two cells $K$ and $K'$ sharing this face, we define the jump operator
    \[
    [\![\partial_nu]\!] \doteq \mathbf{n}_K \cdot \nabla u|_K + \mathbf{n}_{K'} \cdot \nabla u|_{K'},
    \]
    $h_F$ is some average of $h_K$ and $h_{K'}$. \\
    We define the GP stabilisation term:
    \begin{equation*}
        s_h(u,v) = \sum_{F \in \mathcal{F}^{\mathrm{ghost}}_h}{ \left( \lambda_G h_F [\![\partial_nu]\!], [\![\partial_nv]\!] \right)_F}\;.
    \end{equation*}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 3: CutFEM求解}
    \begin{minted}{julia}
        # 1. Build background mesh
        nn = 40
        partition = (nn,nn)
        pmin = 1.2*Point(-1,-1)
        pmax = 1.2*Point(1,1)
        bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
        
        # 2. Build CSG geometry
        R = 1.0
        geo = disk(R, name="csg")
        
        # 3. Cut the background model
        cutgeo = cut(bgmodel,geo)
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 3: CutFEM求解}
    \begin{minted}{julia}
        # 生成计算域
        Ω = Triangulation(cutgeo,PHYSICAL,"csg")
        Ω_act = Triangulation(cutgeo,ACTIVE,"csg")
    \end{minted}
    \begin{figure}
        \centering %居中
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/006.png}
            \caption{物理域$\Omega$}
        \end{minipage}%
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/007.png}
            \caption{求解域$\Omega_{\mathrm{act}}$}
        \end{minipage}
    \end{figure}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 3: CutFEM求解}
    \begin{minted}{julia}
        Γd = EmbeddedBoundary(cutgeo,"csg")
        Γg = GhostSkeleton(cutgeo,"csg")
    \end{minted}
    \begin{figure}
        \centering %居中
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/011.png}
            \caption{背景网格$\Omega_{\mathrm{bg}}$}
        \end{minipage}%
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/010.png}
            \caption{几何}
        \end{minipage}
    \end{figure}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 3: CutFEM求解}
    \begin{minted}{julia}
        # 方程弱形式
        a(u,v) = ∫( ∇(v)⋅∇(u) ) dΩ +
        ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) dΓd +
        ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) dΓg
        
        l(v) = ∫( v*f ) dΩ +
        ∫( (γd/h)*v*g - (n_Γd⋅∇(v))*g ) dΓd
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 4: 一个更复杂的例子}
    \begin{minted}{julia}
        R = 0.5
        geo1 = cylinder(R,v=VectorValue(1,0,0))
        geo2 = cylinder(R,v=VectorValue(0,1,0))
        geo3 = cylinder(R,v=VectorValue(0,0,1))
        geo4 = union(union(geo1,geo2),geo3,name="source")
        
        geo5 = sphere(1)
        geo6 = cube(L=1.5)
        geo7 = intersect(geo6,geo5)
        geo8 = setdiff(geo7,geo4,name="csg")
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 4: 一个更复杂的例子}
    \begin{figure}
        \centering %居中
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/Openscad.png}
            \caption{几何定义}
        \end{minipage}%
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/013.png}
            \caption{分析结果}
        \end{minipage}
    \end{figure}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 4: 一个更复杂的例子}
    \begin{figure}
        \centering %居中
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/012.png}
            \caption{背景网格}
        \end{minipage}%
        \begin{minipage}[t]{0.5\linewidth}
            \centering
            \includegraphics[height=0.5\textheight]{./img/014.png}
            \caption{边界条件}
        \end{minipage}
    \end{figure}
\end{frame}

\section{Part III: 线弹性结构有限元}

\begin{frame}{线弹性有限元}
    Let $\Omega$ be a domain in $\mathbb{R}^d$, $d$ = 2 or 3, with boundary $\partial \Omega = \Gamma_D \cup \Gamma_N$, $\Gamma_D \cap \Gamma_N = \varnothing$, and exterior unit normal $\mathbf{n}$. We consider the following problems: \\
    Find the displacement $\mathbf{u} : \Omega \to \mathbb{R}^d$ such that
    \begin{align*}
        -\nabla \cdot \boldsymbol{\sigma} (\mathbf{u}) &= \mathbf{f} \quad \mathrm{in} \enspace \Omega \\
        \boldsymbol{\sigma} (\mathbf{u}) \cdot \mathbf{n} &= \mathbf{T} \quad \mathrm{on} \enspace \Gamma_N \\
        \mathbf{u} &= \mathbf{g} \quad \mathrm{on} \enspace \Gamma_D
    \end{align*}
    
    where the strain and stress tensors are defined by
    \begin{align*}
        \boldsymbol{\varepsilon} (\mathbf{u}) &\doteq \dfrac{1}{2} \left( \nabla \cdot \mathbf{u} + (\nabla \cdot \mathbf{u} )^{\mathrm{T}} \right) \\
        \sigma (\mathbf{u}) &\doteq 2\mu \boldsymbol{\varepsilon} (\mathbf{u}) + \lambda \mathrm{tr}(\boldsymbol{\varepsilon} (\mathbf{u}))
    \end{align*}
    with Lamé parameters $\lambda$ and $\mu$.
    
\end{frame}

\begin{frame}{线弹性有限元}
    Let $\mathbf{V}_g = \left\{ \mathbf{u} \in \mathbf{H}^1(\Omega) : \mathbf{u} = \mathbf{g} \; \mathrm{on} \; \Gamma_D\right\}$, and define the bilinear form
    \begin{equation*}
        a(\mathbf{u},\mathbf{v})=2\mu \left( \boldsymbol{\varepsilon}(\mathbf{u}), \boldsymbol{\varepsilon}(\mathbf{v})\right)_\Omega + \lambda \left( \mathrm{tr}(\boldsymbol{\varepsilon}(\mathbf{u})), \mathrm{tr}(\boldsymbol{\varepsilon}(\mathbf{v}))\right)_\Omega
    \end{equation*}
    Find such that
    \begin{equation*}
        a(\mathbf{u}, \mathbf{v})=l(\mathbf{v}), \quad \forall \mathbf{v} \in \mathbf{V}_0
    \end{equation*}
    where the linear form on right hand side is defined by
    \begin{equation*}
        l(\mathbf{v})=\left( \mathbf{f} , \mathbf{v}\right)_\Omega + \left( \mathbf{T} , \mathbf{v}\right)_{\Gamma_N}
    \end{equation*}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 6: 悬臂梁}
    \begin{minted}{julia}
        # 读取网格文件
        model = GmshDiscreteModel("./src/beam.msh")
    \end{minted}
    \begin{figure}
        \centering %居中
        \includegraphics[width=1.0\textwidth]{./img/015.png}
        \caption{网格划分}
    \end{figure}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 6: 悬臂梁}
    \begin{minted}{julia}
        # 参考单元
        order = 1
        reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
        
        V = TestFESpace(model,reffe; conformity=:H1,
        dirichlet_tags=["fixed"],
        dirichlet_masks=[(true,true,true)])
        
        g(x) = VectorValue(0.0,0.0,0.0)
        U = TrialFESpace(V, g)
        
        degree = 2*order
        Ω = Triangulation(model)
        dΩ = Measure(Ω,degree)
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 6: 悬臂梁}
    \begin{minted}{julia}
        # 材料参数及本构关系
        const E = 210.0
        const ν = 0.3
        const λ = (E*ν)/((1+ν)*(1-2*ν))
        const μ = E/(2*(1+ν))
        
        σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
    \end{minted}
\end{frame}

\begin{frame}[fragile]
    \frametitle{Example 6: 悬臂梁}
    \begin{minted}{julia}
        # 弱形式
        a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )dΩ
        
        f(x) = VectorValue(0.0,-9.8,0.0)
        l(v) = ∫(v ⋅ f)dΩ
        
        # 求解及输出结果
        op = AffineFEOperator(a,l,U,V)
        uh = solve(op)
        
    \end{minted}
\end{frame}

\begin{frame}{线弹性有限元}
    Define the stabilized Nitsche form
    \begin{equation*}
        \mathcal{A}_h(\mathbf{v},\mathbf{w})=\mathcal{a}_h(\mathbf{v},\mathbf{w})-\left(\boldsymbol{\sigma}(\mathbf{v})\cdot\mathbf{n}, \mathbf{w} \right)_{\Gamma_D} - \left(\mathbf{v}, \boldsymbol{\sigma}(\mathbf{w})\cdot\mathbf{n}\right)_{\Gamma_D} + \beta h^{-1}b_h(\mathbf{v},\mathbf{w})
    \end{equation*}
    where $\beta > 0$ is a parameter and
    \begin{equation*}
        b_h(\mathbf{v}, \mathbf{w})=2\mu\left(  \mathbf{v}, \mathbf{w} \right)_{\Gamma_D} + \lambda\left( \mathbf{v} \cdot \mathbf{n}, \mathbf{w} \cdot \mathbf{n} \right)_{\Gamma_D}
    \end{equation*}
    
\end{frame}

\begin{frame}{线弹性有限元}
    Find $\mathbf{u}_h \in \mathbf{V}_h$ such that
    \begin{equation*}
        \mathcal{A}_h(\mathbf{u}_h, \mathbf{v})=\mathcal{L}_h(\mathbf{v}), \quad \forall\mathbf{v}_h \in \mathbf{V}_h
    \end{equation*}
    where the right hand side is given by
    \begin{equation*}
        \mathcal{L}_h(\mathbf{v})=\left( \mathbf{f} , \mathbf{v}\right)_\Omega + \left( \mathbf{T} , \mathbf{v}\right)_{\Gamma_N} - \left( \mathbf{g}, \boldsymbol{\sigma}(\mathbf{v}) \cdot \mathbf{n} \right)_{\Gamma_D} + \beta h^{-1}b_h\left(\mathbf{g}, \mathbf{v}  \right)_{\Gamma_D}.
    \end{equation*}
\end{frame}

\begin{frame}
    \frametitle{Example 6: 悬臂梁}
    \begin{minted}{julia}
        内容...
    \end{minted}
\end{frame}

\section{Part IV: 杂项}
\section{Part V: 项目规划}
