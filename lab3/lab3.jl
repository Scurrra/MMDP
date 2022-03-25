### A Pluto.jl notebook ###
# v0.18.4

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

# ╔═╡ d4131f40-a5f0-11ec-3145-277b974115e6
begin
	import Pkg;
	Pkg.activate();

	using DataFrames;

	using RecipesBase;
	using Plots
	#using Plots: plot, plot!, @recipe, @animate, @gif, scatter;
	using CairoMakie: Point, Figure, Axis, scatter!, streamplot, streamplot!, (..);
	#using ImplicitEquations

	using DifferentialEquations: ODEProblem, SecondOrderODEProblem, solve;
	using SymPy: @syms, Differential, dsolve, diff, integrate;
	#using ModelingToolkit #: @variables, @parameters, Differential, ODESystem
	
	abstract type Pendulum end;

	const G = 9.807
	
	using PlutoUI;
	TableOfContents()
end

# ╔═╡ 6ad86a53-633c-4226-aab9-62ac0455faee
md"""
# Лабораторная работа №3
# Нелинейные математические модели колебательных явлений
Подготовил: Боровский Илья


## Задание 1. Математический маятник под действием силы тяжести без учета сопротивления среды

Математическая модель представляет собой задачу Коши для нелинейного однородного дифференциального уравнения второго порядка следующего вида:

$\frac{ⅆ^2 α(t)}{ⅆ t^2} + \frac{g}{l} sin( α(t) ) = 0, α(0) = α₀, \frac{ⅆ α(0)}{ⅆ t} = ω₀$

В случае малых колебаний маятника математическая модель (1) упрощается за счет допущения, что $sin(α) ≈ α$ при $α << 1$ и сводится к модели гармонического осциллятора:

$\frac{ⅆ^2 α(t)}{ⅆ t^2} + \frac{g}{l} α(t) = 0, α(0) = α₀, \frac{ⅆ α(0)}{ⅆ t} = ω₀$

Динамическая система запишется в виде:

$\begin{cases}
	\frac{ⅆ α}{ⅆ t} = ω \\
	\frac{ⅆ ω}{ⅆ t} = - \frac{g}{l} sin(α)
\end{cases}$

### Задание 1.1 (Фазовый портрет)

Из записи динамической системы имеем: 

$\frac{ⅆ ω}{ⅆ α} = - \frac{g}{l ω} sin(α) \Rightarrow ω^2 = 2 \left( \frac{g}{l} cos(α) + C \right)$ -- уравнение фазовых траекторий.
"""

# ╔═╡ 9b142b22-d7f7-404f-9cda-f85b912df686
md"""
l = $(@bind L Slider(.1:.1:20, show_value=true, default=9.8))

C = $(@bind C Slider(0:.1:10, show_value=true))
"""

# ╔═╡ 80a8c108-df13-4e26-b46a-266e5bb654e0
streamplot(
	(α, ω) -> Point(
		ω, 
		- G / L * sin(α)
	), -2π..2π, -2π..2π
)

# ╔═╡ 53be9958-f14a-4821-949b-653b09bed530
begin
	Α = collect(-2π:(π/180):2π)
	plot(
		xlims=(-2π, 2π),
		ylims=(-2π, 2π),
		title="Фазовый портрет"
	)
	plot!(
		Α,
		sqrt.(abs.(2 * (G/L*cos.(Α) .+ C))),
		label=:none,
		color=:blue
	)
	plot!(
		Α,
		-sqrt.(abs.(2 * (G/L*cos.(Α) .+ C))),
		label=:none,
		color=:blue
	)
end

# ╔═╡ 707e80ba-6c79-4d52-ade5-7b23b5c9d342
md"""
$ω = 0, α = πn (n ∈ \mathbb(Z))$ -- особые точки системы. 

Если $C = - \frac{g}{l}$ вырождаются в точки с координатами $(2 \pi n, 0)$ -- это подмножество особых точек системы, в которых состояние системы устойчиво.

Если $|C| < \frac{g}{l}$ -- траектории являются замкнутыми кривыми -- колебательные движения.

Если $|C| > \frac{g}{l}$, то $\frac{g}{l} cos(α) + C > 0$ при любом значении α -- траектории являются волнистыми линиями -- вращательные движения.

Если $C = \frac{g}{l}$ -- сепаратрисы.
"""

# ╔═╡ ceafbdd3-ea9e-48c9-b80e-3cb7e29c4536
begin
	plot(
		xlims=(-2π, 2π),
		ylims=(-2π, 2π),
		title="Фазовый портрет"
	)
	
	plot!(
		Α,
		sqrt.(abs.(2 * (G/L*cos.(Α) .- G/L))),
		label="Устойчиво",
		color=:blue
	)
	plot!(
		Α,
		-sqrt.(abs.(2 * (G/L*cos.(Α) .- G/L))),
		label=:none,
		color=:blue
	)

	plot!(
		Α,
		sqrt.(abs.(2 * (G/L*cos.(Α) .- (0.5 * G/L)))),
		label="Колебания",
		color=:red
	)
	plot!(
		Α,
		-sqrt.(abs.(2 * (G/L*cos.(Α) .- (0.5 * G/L)))),
		label=:none,
		color=:red
	)

	plot!(
		Α,
		sqrt.(abs.(2 * (G/L*cos.(Α) .- (1.5 * G/L)))),
		label="Вращение",
		color=:yellow
	)
	plot!(
		Α,
		-sqrt.(abs.(2 * (G/L*cos.(Α) .- (1.5 * G/L)))),
		label=:none,
		color=:yellow
	)

	plot!(
		Α,
		sqrt.(abs.(2 * (G/L*cos.(Α) .+ G/L))),
		label="Сепаратриса",
		color=:green
	)
	plot!(
		Α,
		-sqrt.(abs.(2 * (G/L*cos.(Α) .+ G/L))),
		label=:none,
		color=:green
	)
end

# ╔═╡ a7b8d2b0-fcfe-4715-99b6-4b8bbd03e9bc
md"""
### Задание 1.2 (Динамическая визуализация)
"""

# ╔═╡ e88eee60-2477-4c99-9038-e0dbfe4ced2e
md"""
l = $(@bind L_vis Slider(.1:.1:20, show_value=true, default=9.8))

ω₀ = $(@bind ω₀ Slider(-π:.1π:π, show_value=true, default=0))

α₀ = $(@bind α₀ Slider(-π:.1π:π, show_value=true, default=0))
"""

# ╔═╡ 52cb2d17-a5be-46af-bb05-cd918df184c0
begin
	sol_vis = SecondOrderODEProblem(
		(dα, α, p, t) -> - G/L_vis*sin(α),
		ω₀,
		α₀,
		(0, 1000)
	) |> solve;
	
	@gif for i=1:length(sol_vis.t)
		plot(
			xlims=(-1.3L, 1.3L), 
			ylims=(-1.3L, 1.3L),
			title="Math Pendulum"
		)
		Plots.scatter!(
			[L_vis*sin(sol_vis[i][2])], 
			[-L_vis*cos(sol_vis[i][2])],
			label=:none
		)
		plot!(
			[0, L_vis*sin(sol_vis[i][2])], 
			[0, -L_vis*cos(sol_vis[i][2])],
			label=:none
		)
	end
end

# ╔═╡ f5786b75-496c-4821-be80-1e724d65dc9d
md"""
### Задание 1.3 (Период колебаний)

Запишем закон сохранения энергии для процесса колебаний математического маятника:

$\frac{m v^2}{2} = m g l (cos(α) - cos(α₀))$

Пусть длина дуги, по которой тело отклоняется от положения равновесия, вычисляется при малых углах α следующим образом: $s = l α ⇒ v = \frac{ⅆ s}{ⅆ t} = l \frac{ⅆ α}{ⅆ t}$

Подставим в закон сохранения энергии.

$\frac{\left(l \frac{ⅆ α}{ⅆ t} \right)^2}{2} = g l (cos(α) - cos(α₀))$

$\frac{l}{2} \left(\frac{ⅆ α}{ⅆ t} \right)^2 = g (cos(α) - cos(α₀))$

$ⅆ t = - √\left( \frac{l}{2 g} \right) \frac{ⅆ α}{√( cos(α) - cos(α₀) )}$

Для колебания маятника $T$ имеем:

$T = 4 √\left( \frac{l}{2 g} \right) \int \limits_{0}^{α} \frac{ⅆ α}{√( cos(α) - cos(α₀) )}$


"""

# ╔═╡ 32e566b2-9951-48d3-a9a0-fac330ea46ce
@syms t, g, l, φ, φ₀

# ╔═╡ a9a8daab-341c-416b-a898-9abf66c32d3a
solve(integrate(
	1 / sqrt( cos(φ) - cos(φ₀) ), (φ, 0, φ₀)
) * 4 * sqrt( l / 2 / g ) ~ 2, l)

# ╔═╡ 260045bb-5f5f-4049-911c-0561d31699a9
md"""
## Задание 2. Модель двухвидового взаимодействия “хищник-жертва”

Система уравнений Лотки-Вольтерра задает математическую модель взаимодействия популяции жертв с численностью $N = N(t) ≥ 0$ и популяции хищников с численностью $M = M(t) ≥ 0$:

$\begin{cases}
	\frac{ⅆ N}{ⅆ t} = α N - c N M \\
	\frac{ⅆ M}{ⅆ t} = - β M + d N M
\end{cases} \; , N(t₀) = N₀, M(t₀) = M₀$

где
$α = const > 0$ -- коэффициент прироста жертв при отсутствии хищников в условиях неограниченности ресурса для питания;
$c = const > 0$ -- коэффициент пропорциональности при взаимодействии жертв с хищниками, соответствующий уменьшению жертв;
$β = const > 0$ -- коэффициент смертности хищников при отсутствии жертв;
$d = const > 0$ -- коэффициент пропорциональности при взаимодействии жертв с хищниками, соответствующий увеличению жертв.
"""

# ╔═╡ 3be58d11-b11c-4e1d-9516-78066f63451b
md"""
### Задание 2.1 (Колебательное поведение)
"""

# ╔═╡ fe09e8dc-e05f-426b-a0bb-e9de46a69965
md"""
α = $(@bind α Slider(0:0.1:1, show_value=true, default=0.2))

β = $(@bind β Slider(0:0.1:1, show_value=true, default=0.1))

c = $(@bind c Slider(0:0.1:1, show_value=true, default=0.3))

d = $(@bind d Slider(0:0.1:1, show_value=true, default=0.4))

N₀ = $(@bind N₀ NumberField(0:100, default=10))

M₀ = $(@bind M₀ NumberField(0:100, default=10))
"""

# ╔═╡ bc1633bc-7806-4933-ba78-1307b91a715d
begin
	function pp!(du, u, p, t)
		du[1] = α * u[1] - c * u[1] * u[2]
		du[2] = -β * u[2] + d * u[1] * u[2]
	end
		
	struct PredatorPrey
		α
		β
		c
		d

		N₀
		M₀

		solution

		PredatorPrey(α=α, β=β, c=c, d=d, N₀=N₀, M₀=M₀; tmax=1000) = new(
			α, β, c, d, N₀, M₀,
			ODEProblem(
				pp!,
				[N₀, M₀],
				(0, tmax)
			) |> solve
		)
	end

	@recipe function f(pp::PredatorPrey)
		layout := 2
		widths := [0.5, 0.5]
		
		@series begin
			subplot := 1
			xlabel := "Predators"
			ylabel := "Preys"
			label := :none
			vars := (1,2)
			pp.solution
		end

		@series begin
			subplot := 2
			xlabel := "Predators"
			ylabel := "Preys"
			label := :none
			pp.solution
		end
	end

	function streamplot_pp(pp::PredatorPrey)
		fig = Figure()
		ax = Axis(fig[1, 1])
		streamplot!(
			(N, M) -> Point(
				pp.α * N - pp.c * N * M,
				- pp.β * M + pp.d * N * M 
			), 0..10, 0..10
		)
		scatter!([0, pp.β/pp.d], [0, pp.α/pp.c])
		fig
	end

	pp = PredatorPrey()
end;

# ╔═╡ 4bff40df-cac3-48ac-84bd-84460b59575d
pp |> plot

# ╔═╡ 811ce2fc-4a5f-4b9f-92bf-1b728d0f6e63
pp |> streamplot_pp

# ╔═╡ 85a92feb-fc0c-4357-94b6-19fe9f99a0e1
md"""
### Задание 2.2 (Неявная зависимость между N(t) и M(t))
"""

# ╔═╡ 96412f1d-3768-45c5-a071-11f8372b7a26
md"""
### Задание 2.3 (Структурная неустойчивость модели “хищник-жертва”)
"""

# ╔═╡ 2c56071e-aea8-473c-a1de-fb6a91fe1358
begin
	fig = Figure()
	ax = Axis(fig[1, 1])
	streamplot!(
		(N, M) -> Point(
			pp.α * N - pp.c * N * M + ℯ^(-N^10) * M,
			- pp.β * M + pp.d * N * M + ℯ^(-M^10) * N
		), -10..10, -10..10
	)
	scatter!([0, pp.β/pp.d], [0, pp.α/pp.c])
	fig
end

# ╔═╡ Cell order:
# ╟─d4131f40-a5f0-11ec-3145-277b974115e6
# ╟─6ad86a53-633c-4226-aab9-62ac0455faee
# ╟─9b142b22-d7f7-404f-9cda-f85b912df686
# ╟─80a8c108-df13-4e26-b46a-266e5bb654e0
# ╟─53be9958-f14a-4821-949b-653b09bed530
# ╟─707e80ba-6c79-4d52-ade5-7b23b5c9d342
# ╟─ceafbdd3-ea9e-48c9-b80e-3cb7e29c4536
# ╟─a7b8d2b0-fcfe-4715-99b6-4b8bbd03e9bc
# ╟─e88eee60-2477-4c99-9038-e0dbfe4ced2e
# ╟─52cb2d17-a5be-46af-bb05-cd918df184c0
# ╟─f5786b75-496c-4821-be80-1e724d65dc9d
# ╠═32e566b2-9951-48d3-a9a0-fac330ea46ce
# ╠═a9a8daab-341c-416b-a898-9abf66c32d3a
# ╟─260045bb-5f5f-4049-911c-0561d31699a9
# ╟─3be58d11-b11c-4e1d-9516-78066f63451b
# ╟─fe09e8dc-e05f-426b-a0bb-e9de46a69965
# ╟─bc1633bc-7806-4933-ba78-1307b91a715d
# ╟─4bff40df-cac3-48ac-84bd-84460b59575d
# ╠═811ce2fc-4a5f-4b9f-92bf-1b728d0f6e63
# ╟─85a92feb-fc0c-4357-94b6-19fe9f99a0e1
# ╟─96412f1d-3768-45c5-a071-11f8372b7a26
# ╟─2c56071e-aea8-473c-a1de-fb6a91fe1358
