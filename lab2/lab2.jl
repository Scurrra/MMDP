### A Pluto.jl notebook ###
# v0.18.2

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

# ╔═╡ 24859bea-a037-11ec-1529-ede51f9bb656
begin
	import Pkg;
	Pkg.activate();

	using DataFrames;

	using RecipesBase;
	using Plots: plot, plot!, @recipe;

	using DifferentialEquations: SecondOrderODEProblem, solve;
	
	abstract type Pendulum end;
	
	using PlutoUI;
	TableOfContents()
end

# ╔═╡ 45d13863-edb9-4eb1-879d-8d4221deca50
md"""
# Лабораторная работа №2
# Линейные математические модели колебательных явлений
Подготовил: Боровский Илья
"""

# ╔═╡ b0cfb272-ebe8-444a-b9e2-7ed8d3979003
md"""
## Задание 1. Гармонические колебания

$\frac{∂^2 r(t)}{∂ t^2} + ω^2 r(t) = 0, r(0) = r₀, \frac{∂ r(0)}{∂ t} = v₀$

где $ω = \sqrt{k / m} = const > 0$ -- собственная частота колебаний.

Решим данную модель. Характеристическое уравнение имеет вид:

$λ^2 + ω^2 = 0 \Rightarrow λ = ±i√ω$ -- центр.

### Задание 1.1 (Амплитуда колебаний)

Решение представимо в виде: $r(t) = A cos(ω t + ϕ₀)$

Пусть $L$ -- растояние до стенки ⇒ $A ≤ L$. Нужно найти условия на величины начальных данных $r₀$ и $v₀$.

$\begin{cases}
	r₀ = A cos(ϕ₀) \\
	v₀ = -A ω sin(ϕ₀)
\end{cases}$
где $ϕ₀ ≤ arcsin \left( \frac{A}{l} \right)$, $l$ -- длина нити.

В случае $r₀ ≡ 0$ имеем
$\begin{cases}
	r₀ = A \\
	v₀ = ± A ω
\end{cases}$
где знак $v₀$ задаёт направление скорости.

### Задание 1.2 (Динамическая визуализация с учетом соударения тела о стенку)
"""

# ╔═╡ 22ef3e0d-c2d7-4806-8d80-793f9ce960b2
md"""
ω = $(@bind ω Slider(0.1:0.1:1, show_value=true, default=1))

r₀ = $(@bind r₀ Slider(-1:0.1:1, show_value=true, default=0))

v₀ = $(@bind v₀ Slider(-1:0.1:1, show_value=true, default=0))

L = $(@bind L Slider(0.1:0.1:1, show_value=true, default=1))
"""

# ╔═╡ 7abf4e62-af22-4906-9e23-4f18dd7f67ab
begin
	struct NRPendulum <: Pendulum
		ω
		
		r₀
		v₀

		L

		solution

		function NRPendulum(ω::Number, r₀::Number, v₀::Number, L::Number; tspan=(0, 10))

			sol = SecondOrderODEProblem(
				(dr, r, p, t) -> - ω^2 * r,
				v₀,
				r₀,
				tspan
			) |> solve
			
			new(ω, r₀, v₀, L,
				sol
			)

		end
	end

	@recipe function f(pendulum::Pendulum)
		ylims := extrema(u -> u[2], pendulum.solution.u) .+ (-0.1-pendulum.L, 0.1)
		
		@series begin
			label := "Pendulum"
			vars := 2
			pendulum.solution
		end
		
		@series begin
			seriestype := :hline
			label := "Wall"
			linewidth := 5
			[-pendulum.L]
		end
	end
end

# ╔═╡ 03ecedd8-5700-4f61-b3e8-72c042c73c11
sol1 = NRPendulum(ω, r₀, v₀, L)

# ╔═╡ 5222ff14-5197-4af7-801d-bf652fa8a88a
sol1 |> plot

# ╔═╡ 167ae2e2-0e79-4bb6-a12d-f11f0a5ff0a9


# ╔═╡ Cell order:
# ╟─24859bea-a037-11ec-1529-ede51f9bb656
# ╟─45d13863-edb9-4eb1-879d-8d4221deca50
# ╟─b0cfb272-ebe8-444a-b9e2-7ed8d3979003
# ╟─22ef3e0d-c2d7-4806-8d80-793f9ce960b2
# ╠═7abf4e62-af22-4906-9e23-4f18dd7f67ab
# ╠═03ecedd8-5700-4f61-b3e8-72c042c73c11
# ╠═5222ff14-5197-4af7-801d-bf652fa8a88a
# ╠═167ae2e2-0e79-4bb6-a12d-f11f0a5ff0a9
