### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 51564fca-b572-11ec-1f86-2b443d5a1085
begin
	import Pkg;
	Pkg.activate();

	using DataFrames;

	using RecipesBase;
	using Plots; gr();
	#using Plots: plot, plot!, @recipe, @animate, @gif, scatter;
	using CairoMakie: Point, Figure, Axis, scatter!, streamplot, streamplot!, vlines!, (..);
		
	using DifferentialEquations: ODEProblem, solve;

	abstract type AbstractSIR end

		
	using PlutoUI;
	TableOfContents()
end

# ╔═╡ 030aecea-ed3d-403b-97e7-4d7916382556
md"""
# Лабораторная работа №4
# Математические модели распространения инфекционных заболеваний
Подготовил: Боровский Илья
"""

# ╔═╡ 4b079bbe-51e4-4714-aa48-36a214e0f20c
md"""
## Задание 1. SIR-модель

$\begin{cases}
	\frac{ⅆ S}{ⅆ t} = - ir \; S \; I \\
	\frac{ⅆ I}{ⅆ t} = ir \; S \; I - rr \; I \\
	\frac{ⅆ R}{ⅆ t} = rr \; I
\end{cases}$

$\begin{cases}
	S(0) = N - I₀ \\
	I(0) = I₀ \\
	R(0) = 0
\end{cases}$

где $ir = const > 0$ -- infection rate, $rr = const > 0$ -- recovery rate
"""

# ╔═╡ 6e458c4c-bf1a-4b16-9d47-ace152ce2abd
md"""
### Задание 1.1 (Пороговый эффект)
"""

# ╔═╡ 64012f30-afa5-4142-9a3a-1d8819b3e23c
begin 
	function sir!(du, u, p, t)
		du[1] = - p[1] * u[1] * u[2]
		du[2] = p[1] * u[1] * u[2] - p[2] * u[2]
		du[3] = p[2] * u[2]
	end
	
	struct SIR <: AbstractSIR
		ir::Float64
		rr::Float64
		ρ::Float64

		N::Int

		S₀::Int
		I₀::Int
		R₀::Float64

		solution

		SIR(ir::Float64, rr::Float64; N::Int=10^4, I₀::Int=10^2, tmax::Int=10^3) = new(
			ir, rr, rr / ir,
			N, 
			N - I₀,
			I₀, 
			ir / rr * (N - I₀),
			ODEProblem(
				sir!,
				[N - I₀, I₀, 0],
				(0, tmax),
				(ir, rr)
			) |> solve
		)
	end

	@recipe function f(sir::SIR)
		title := "SIR-model, $(sir.R₀)"
		label := ["S(t)" "I(t)" "R(t)"]
		sir.solution
	end

	plot3d(sir::SIR) = plot(
		sir.solution,
		vars=(1,2,3),
		label=:none,
		xlabel="S",
		ylabel="I",
		zlabel="R",
		title="S-I-R"
	)

	function streamplot_sir(sir::SIR)
		fig = Figure()
		ax = Axis(fig[1, 1])
		streamplot!(
			(S, I) -> Point(
				- sir.ir * S * I,
				sir.ir * S * I - sir.rr * I 
			), 0..(sir.N), 0..(sir.N)
		)
		vlines!(
			ax,
			sir.ρ
		)
		scatter!(
			[sir.S₀],
			[sir.I₀]
		)
		fig
	end
end;

# ╔═╡ a90bcd6e-bfaa-4cd4-ab76-c0c0619f40b1
md"""
ir = $(@bind ir_sir Slider(0.01:0.01:1; default=0.1, show_value=true))

rr = $(@bind rr_sir Slider(0.01:0.01:1; default=0.2, show_value=true))
"""

# ╔═╡ a4bbcf81-460a-45e9-8ecd-a74537ed5e36
sir = SIR(ir_sir, rr_sir; N=100, I₀ = 10, tmax=10);

# ╔═╡ 7ef9886d-1321-41c7-a89d-bbcfa5ac47ec
sir |> plot

# ╔═╡ e52c8a17-40c6-4d1c-8140-f1ec68e607f0
streamplot_sir(sir)

# ╔═╡ 281f448f-37f8-4ce1-89a1-1d7ecd39f72b
md"""
### Задание 1.2 (Учет вакцинации в модели)
"""

# ╔═╡ 1175e548-7b7c-487a-b616-c76b24846b4a
begin 
	struct SIRV <: AbstractSIR
		ir::Float64
		rr::Float64
		ρ::Float64
		p::Float64

		N::Int

		S₀::Int
		I₀::Int
		R₀::Float64

		solution

		function SIRV(ir::Float64, rr::Float64; N::Int=10^4, I₀::Int=10^2, tmax::Int=10^3, ϵ::Float64=0.01)
			p = 1 - 1 / (ir / rr * (N - I₀)) + ϵ
			new(
				ir, rr, rr / ir, p,
				N, 
				N - I₀,
				I₀, 
				ir / rr * (N - I₀),
				ODEProblem(
					sir!,
					[(1-p)*N - I₀, I₀, p*N],
					(0, tmax),
					(ir, rr)
				) |> solve
			)
		end
	end

	@recipe function f(sirv::SIRV)
		title := "SIR-model, $(sirv.R₀), p=$(sirv.p)"
		label := ["S(t)" "I(t)" "R(t)"]
		sirv.solution
	end
end

# ╔═╡ 27843c46-68e3-42b0-9787-5e2e9bd37d7f
SIRV(ir_sir, rr_sir; tmax=1) |> plot

# ╔═╡ 40803d60-3d88-4bd3-92b2-300d2b0863b1
md"""
### Задание 1.3 (Коллективный иммунитет)

$\begin{cases}
	\frac{ⅆ S}{ⅆ t} = - ir \; S \; I = S' \\
	\frac{ⅆ I}{ⅆ t} = ir \; S \; I - rr \; I = - S' - R'\\
	\frac{ⅆ R}{ⅆ t} = rr \; I = R'
\end{cases} \;\;\; \begin{cases}
	S(0) = N - I₀ \\
	I(0) = I₀ \\
	R(0) = 0
\end{cases}$

$\begin{cases}
	- ir \; S \; I = S' \\
	rr \; I = R' \\ 
	N = S + I + R
\end{cases} ⇒ \begin{cases}
	- ir \; S \; I = S' \\
	rr \; I = R' \\ 
	S = N - I - R
\end{cases} ⇒$

$⇒ \begin{cases}
	- ir \; (N - I - R) \; I = S' \\
	I = \frac{R'}{rr}
\end{cases} ⇒ - ir \; (N - I₀ - R₀) \; \frac{R'}{rr} = S' ⇒$

$⇒ R₀ = \frac{ir}{rr}(N - I₀) = -\frac{S'}{R'}$

Имеем: $R_{H3N2} = 18 ⇒ p^⋆_{H3N2} = \frac{17}{18}, R_{measles} = 144 ⇒ p^⋆_{measles} = \frac{143}{144}$
"""

# ╔═╡ f812af9f-31d5-43c7-a269-53c786b22c1b
md"""
## Задание 2. SEIR-модель

$\begin{cases}
	\frac{ⅆ S}{ⅆ t} = - ir \; S \; I \\
	\frac{ⅆ E}{ⅆ t} = ir \; S \; I - er \; E \\
	\frac{ⅆ I}{ⅆ t} = er \; E - rr \; I \\
	\frac{ⅆ R}{ⅆ t} = rr \; I
\end{cases}$

$\begin{cases}
	S(0) = N - I₀ \\
	E(0) = 0 \\
	I(0) = I₀ \\
	R(0) = 0
\end{cases}$

где $ir = const > 0$ -- infection rate, $rr = const > 0$ -- recovery rate, $er = const > 0$ -- exposed rate

"""

# ╔═╡ 5ca102fe-d6a1-42f3-8677-ddf9678a0258
md"""
### Задание 2.1 (Пороговый эффект)
"""

# ╔═╡ 650e613e-cbc2-4f6c-9011-319f406c1030
md"""
ir = $(@bind ir_seir Slider(0.01:0.01:1; default=0.1, show_value=true))

er = $(@bind er_seir Slider(0.01:0.01:1; default=0.6, show_value=true))

rr = $(@bind rr_seir Slider(0.01:0.01:1; default=0.2, show_value=true))
"""

# ╔═╡ 7b7b93cc-6fc4-4beb-8129-5728ce12bb84
begin 
	function seir!(du, u, p, t)
		du[1] = - p[1] * u[1] * u[3] # S
		du[2] = p[1] * u[1] * u[3] - p[2] * u[2] # E
		du[3] = p[2] * u[2] - p[3] * u[3] # I
		du[4] = p[3] * u[3] # R
	end
	
	struct SEIR <: AbstractSIR
		ir::Float64
		er::Float64
		rr::Float64
		ρ::Float64

		N::Int

		S₀::Int
		I₀::Int
		R₀::Float64

		solution

		SEIR(ir::Float64, er::Float64, rr::Float64; N::Int=10^4, I₀::Int=10^2, tmax::Int=10^3) = new(
			ir, rr, er, rr / ir,
			N, 
			N - I₀,
			I₀, 
			ir / rr * (N - I₀),
			ODEProblem(
				seir!,
				[N - I₀, 0, I₀, 0],
				(0, tmax),
				(ir, er, rr)
			) |> solve
		)
	end

	@recipe function f(seir::SEIR)
		title := "SEIR-model, $(seir.R₀)"
		label := ["S(t)" "E(t)" "I(t)" "R(t)"]
		seir.solution
	end

	function streamplot_seir(seir::SEIR)
		fig = Figure()
		ax = Axis(fig[1, 1])
		streamplot!(
			(S, I) -> Point(
				- seir.ir * S * I,
				#=seir.er * E=# - seir.rr * I 
			), 0..(seir.N), 0..(seir.N)#, 0..(seir.N)
		)
		vlines!(
			ax,
			seir.ρ
		)
		scatter!(
			[seir.S₀],
			[seir.I₀]
		)
		fig
	end
end;

# ╔═╡ 89188c94-d889-4c70-831c-d1f7bf7ecfe1
seir = SEIR(ir_seir, er_seir, rr_seir; N=100, I₀=1, tmax=30);

# ╔═╡ 3fbe9e75-eea4-4e88-a881-a1b3e7c583b5
seir |> plot

# ╔═╡ c77a3e78-1ecc-4bfe-b943-b08027bc695e
seir |> streamplot_seir

# ╔═╡ 4b0b8e9e-425b-4548-8bb7-d65458fc3ec4
md"""
### Задание 2.2 (Качественный анализ динамики распространения заболевания)
"""

# ╔═╡ d17fd9e0-79c4-4d26-b5a7-b401934019cf
md"""
## Задание 3. SIR-модель с дополнительным учетом рождаемости и смертности

$\begin{cases}
	\frac{ⅆ S}{ⅆ t} = - ir \; S \; I + br \; N  - dr \; S\\
	\frac{ⅆ I}{ⅆ t} = ir \; S \; I - rr \; I - dr \; I\\
	\frac{ⅆ R}{ⅆ t} = rr \; I - dr \; R
\end{cases}$

$\begin{cases}
	S(0) = N - I₀ \\
	I(0) = I₀ \\
	R(0) = 0
\end{cases}$

где $ir = const > 0$ -- infection rate, $rr = const > 0$ -- recovery rate, $br = const > 0$ -- birth rate, $dr = const > 0$ -- death rate

"""

# ╔═╡ 05aabd28-bc38-4dd2-9bbe-afde5656d746
md"""
### Задание 3.1 (Пороговый эффект)
"""

# ╔═╡ cfa0188e-6e1e-4992-bd05-07f3da27b175
begin 
	function msir!(du, u, p, t)
		du[1] = - p[1] * u[1] * u[2] + p[3] * u[4] - p[4] * u[1]
		du[2] = p[1] * u[1] * u[2] - p[2] * u[2] - p[4] * u[2]
		du[3] = p[2] * u[2] - p[4] * u[3]
		du[4] = (p[3] - p[4]) * u[4]
	end
	
	struct MSIR <: AbstractSIR
		ir::Float64
		rr::Float64
		br::Float64
		dr::Float64
		
		ρ::Float64

		N::Int

		S₀::Int
		I₀::Int
		R₀::Float64

		solution

		MSIR(ir::Float64, rr::Float64, br::Float64, dr::Float64; N::Int=10^4, I₀::Int=10^2, tmax::Int=10^3) = new(
			ir, rr, br, dr, (rr + dr) / ir,
			N, 
			N - I₀,
			I₀, 
			ir / (rr + dr) * (N - I₀),
			ODEProblem(
				msir!,
				[N - I₀, I₀, 0, N],
				(0, tmax),
				(ir, rr, br, dr)
			) |> solve
		)
	end

	@recipe function f(msir::MSIR)
		title := "MSIR-model, $(msir.R₀)"
		label := ["S(t)" "I(t)" "R(t)" "N(t)"]
		msir.solution
	end

	plot3d(msir::MSIR) = plot(
		msir.solution,
		vars=(1,2,3),
		label=:none,
		xlabel="S",
		ylabel="I",
		zlabel="R",
		title="S-I-R"
	)
end;

# ╔═╡ 08b15145-209f-4b0e-a275-2cc242ba62f8
sir |> plot3d

# ╔═╡ 4f72bac7-9232-4f5c-9f37-ddb1e32df5fe
md"""
ir = $(@bind ir_msir Slider(0.01:0.01:1; default=0.01, show_value=true))

rr = $(@bind rr_msir Slider(0.01:0.01:1; default=0.2, show_value=true))

br = $(@bind br_msir Slider(0.01:0.01:1; default=0.2, show_value=true))

dr = $(@bind dr_msir Slider(0.01:0.01:1; default=0.2, show_value=true))
"""

# ╔═╡ d9da28a0-4b82-41e2-8bbf-b0584b7e57d0
msir = MSIR(ir_msir, rr_msir, br_msir, dr_msir; N=100, I₀=1, tmax=30);

# ╔═╡ 982dc642-de8a-4ff4-95b4-d12fea6af0ec
msir |> plot

# ╔═╡ ba351b4d-287b-462c-9fee-05a37cd27b88
msir |> plot3d

# ╔═╡ c37fe7c4-d9ef-40fc-a374-3c9bd93051d5
md"""
### Задание 3.2 (Качественный анализ динамики распространения заболевания)

 - $br < dr$ -- с падением рождаемости падает и скорость заражения.

 - $br = dr$ -- графики выходят на плато.

 - $br > dr$ -- с ростом рождаемости растёт и скорость заражения.
"""

# ╔═╡ Cell order:
# ╟─51564fca-b572-11ec-1f86-2b443d5a1085
# ╟─030aecea-ed3d-403b-97e7-4d7916382556
# ╟─4b079bbe-51e4-4714-aa48-36a214e0f20c
# ╟─6e458c4c-bf1a-4b16-9d47-ace152ce2abd
# ╟─64012f30-afa5-4142-9a3a-1d8819b3e23c
# ╟─a90bcd6e-bfaa-4cd4-ab76-c0c0619f40b1
# ╠═a4bbcf81-460a-45e9-8ecd-a74537ed5e36
# ╠═7ef9886d-1321-41c7-a89d-bbcfa5ac47ec
# ╠═08b15145-209f-4b0e-a275-2cc242ba62f8
# ╠═e52c8a17-40c6-4d1c-8140-f1ec68e607f0
# ╟─281f448f-37f8-4ce1-89a1-1d7ecd39f72b
# ╟─1175e548-7b7c-487a-b616-c76b24846b4a
# ╠═27843c46-68e3-42b0-9787-5e2e9bd37d7f
# ╟─40803d60-3d88-4bd3-92b2-300d2b0863b1
# ╟─f812af9f-31d5-43c7-a269-53c786b22c1b
# ╟─5ca102fe-d6a1-42f3-8677-ddf9678a0258
# ╟─650e613e-cbc2-4f6c-9011-319f406c1030
# ╟─7b7b93cc-6fc4-4beb-8129-5728ce12bb84
# ╠═89188c94-d889-4c70-831c-d1f7bf7ecfe1
# ╠═3fbe9e75-eea4-4e88-a881-a1b3e7c583b5
# ╠═c77a3e78-1ecc-4bfe-b943-b08027bc695e
# ╟─4b0b8e9e-425b-4548-8bb7-d65458fc3ec4
# ╟─d17fd9e0-79c4-4d26-b5a7-b401934019cf
# ╟─05aabd28-bc38-4dd2-9bbe-afde5656d746
# ╟─cfa0188e-6e1e-4992-bd05-07f3da27b175
# ╟─4f72bac7-9232-4f5c-9f37-ddb1e32df5fe
# ╠═d9da28a0-4b82-41e2-8bbf-b0584b7e57d0
# ╠═982dc642-de8a-4ff4-95b4-d12fea6af0ec
# ╠═ba351b4d-287b-462c-9fee-05a37cd27b88
# ╟─c37fe7c4-d9ef-40fc-a374-3c9bd93051d5
