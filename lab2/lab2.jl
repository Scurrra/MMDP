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
	using CairoMakie: Point, streamplot, (..);

	using DifferentialEquations: SecondOrderODEProblem, solve;
	using SymPy: @syms, Differential, solve, diff, dsolve

	using LinearAlgebra
	
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

$\frac{ⅆ^2 r(t)}{ⅆ t^2} + ω^2 r(t) = 0, r(0) = r₀, \frac{ⅆ r(0)}{ⅆ t} = v₀$

где $ω = \sqrt{k / m} = const > 0$ -- собственная частота колебаний.

Решим данную модель. Характеристическое уравнение имеет вид:

$λ^2 + ω^2 = 0 \Rightarrow λ = ±i√ω$ -- центр.

Динамическая система имеет вид:

$\begin{cases}
	\frac{ⅆ r}{ⅆ t} = v \\
	\frac{ⅆ v}{ⅆ t} = -ω^2 t
\end{cases}$

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
ω = $(@bind Ω Slider(0.1:0.1:1, show_value=true, default=1))

r₀ = $(@bind R₀ Slider(-1:0.1:1, show_value=true, default=0))

v₀ = $(@bind V₀ Slider(-1:0.1:1, show_value=true, default=0))

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
sol1 = NRPendulum(Ω, R₀, V₀, L);

# ╔═╡ 5222ff14-5197-4af7-801d-bf652fa8a88a
sol1 |> plot

# ╔═╡ 167ae2e2-0e79-4bb6-a12d-f11f0a5ff0a9
md"""
## Задание 2. Колебания с учетом сопротивления среды

$m \frac{ⅆ^2 r(t)}{ⅆ t^2} = -k r(t) - μ \frac{ⅆ r(t)}{ⅆ t}, r(0) = r₀, \frac{ⅆ r(0)}{ⅆ t} = v₀$

где $μ = const > 0$

Динамическая система имеет вид:

$\begin{cases}
	\frac{ⅆ r}{ⅆ t} = v \\
	\frac{ⅆ v}{ⅆ t} = - \frac{k}{m} r - \frac{μ}{m} v
\end{cases}$

### Задание 2.1 (Качественный анализ по фазовому портрету)

Характерестическое уравнение:

$m λ^2 + μ λ + k = 0$

$λ = \frac{-μ ± √(μ^2 - 4 k m)}{2 m}$

$Re λ_{1,2} < 0 ⇒ \; положение \; равновесия r = 0, v = 0 \; всегда \; является \; устойчивым$
"""

# ╔═╡ 3755c8f9-466a-4cad-be08-385503dbcd6b
md"""
m = $(@bind M Slider(0:1:10, default=1, show_value=true))

k = $(@bind K Slider(0:0.1:1, default=0.1, show_value=true))
"""

# ╔═╡ aed932a2-e563-475d-84a8-2f6d0d35c32e
md"""
#### 1. $μ^2 < 4 k m$ (устойчивый фокус)

$r(t) = ℯ^{\frac{- μ t}{2 m}} \left( A cos \left( \frac{√(4 k m - μ^2) t}{2 m} \right) + B cos \left( \frac{√(4 k m - μ^2) t}{2 m} \right) \right)$

**Затухающие колебания** с экспоненциально убывающей амплитудой. Период колебаний не зависит от начальных данных.
"""

# ╔═╡ 6cd36f0b-f2b3-4885-8a4a-16b917e8aa7e
begin
	μ₁ = √(4 * K * M) - 1
	streamplot(
		(r, v) -> Point(
			v, 
			- K/M * r - μ₁/M * v
		), -10..10, -10..10
	)
end

# ╔═╡ c46c02cd-5d1a-4f8b-aadd-50e223dde894
md"""
#### 2. $μ^2 = 4 k m$ (устойчивый узел)

$r(t) = ℯ^{\frac{- μ t}{2 m}} (С₁ + С₂ t)$

**Критическое затухание.** Решение не является колеблющимся. Оно является формальным разделением колебательного и апериодического поведения. 
"""

# ╔═╡ 8088e6e3-6a90-46d8-9a92-ec143ce55685
begin
	μ₂ = √(4 * K * M)
	streamplot(
		(r, v) -> Point(
			v, 
			- K/M * r - μ₂/M * v
		), -10..10, -10..10
	)
end

# ╔═╡ 7b39fe50-fc4a-4b65-b339-d3f6c039c172
md"""
#### 3. $μ^2 > 4 k m$ (устойчивый узел, так как $ω < μ$)

$r(t) = ℯ^{\frac{-μ t}{2 m}} \left( C₁ ℯ^{\frac{t √(μ^2 - 4 k m)}{2 m}} + C₂ ℯ^{\frac{-t √(μ^2 - 4 k m)}{2 m}} \right)$

**Апериодическое затухание.** Соответствует большому сопротивлению среды.
"""

# ╔═╡ 86e09784-7678-43ac-885e-6d3eedca7f3e
begin
	μ₃ = √(4 * K * M) + 1	
	streamplot(
		(r, v) -> Point(
			v, 
			- K/M * r - μ₃/M * v
		), -10..10, -10..10
	)
end

# ╔═╡ f0eb22f4-b81b-46f3-8157-e99b4c6b5018
md"""
### Задание 2.2 (Качественный анализ по теореме Ляпунова)
"""

# ╔═╡ 12fc6293-a695-4c9c-8682-8892dfa081cb
@syms λ, μ, k, m

# ╔═╡ c2d02d98-3236-4f44-be6a-78524cc44e6c
A = [0 1; -k/m -μ/m]

# ╔═╡ 971c1060-6519-4ca2-a873-6c5fccca3caa
solve(det(A - λ * I) ~ 0, λ)

# ╔═╡ 54a9c715-1e37-4f02-8b61-1378fa48d54b
md"""
Оба имеют отрицательную действительную часть, следовательно решения ассимптотически устойчивы.
"""

# ╔═╡ 168e14d8-933f-445e-a9e6-905611b6e4d7
md"""
### Задание 2.3 (Анализ частных случаев)
"""

# ╔═╡ 529ec0e3-2a4b-422d-949e-ab3a6938ce9c
systype(μ, k, m) = μ^2 >= k * m ? "Устойчивый узел" : "Устойчивый фокус";

# ╔═╡ cea0ab76-7faf-48c3-b6e3-421c76a63082
md"""
В воздухе: $(systype(1.82*10^-5, 10, 10))

В воде: $(systype(1.002*10^-3, 10, 10))

В глицерине: $(systype(1.49, 10, 10))
"""

# ╔═╡ 118ae518-5e84-4884-af6b-14cd398aca03
md"""
## Задание 3. Вынужденные колебания

$m \frac{ⅆ^2 r(t)}{ⅆ t^2} = -k r(t) + F₀ sin(t ω), r(0) = r₀, \frac{ⅆ r(0)}{ⅆ t} = v₀$

где $μ = const > 0, ω₂ = const > 0$

Динамическая система имеет вид:

$\begin{cases}
	\frac{ⅆ r}{ⅆ t} = v \\
	\frac{ⅆ v}{ⅆ t} = - \frac{k}{m} r + \frac{F₀ sin(t ω)}{m}
\end{cases}$

### Задание 3.1 (Частное решение математической модели)
#### Резонансный случай

$r₁^*(t) = t (C₁ sin(ω t) + C₂ cos(ω t))$ где $ω = √\left( \frac{k}{m} \right)$ 

$\dot{r₁^*} = (C₁ sin(ω t) + C₂ cos(ω t)) + ω t (C₁ cos(ω t) - C₂ sin(ω t)) = (C₁ - C₂ ω t) sin(ω t) + (C₁ + C₂ ω t) cos(ω t)$

$\ddot{r₁^*} = 2 ω (C₁ cos(ω t) - C₂ sin(ω t)) - ω^2 t (C₁ sin(ω t) + C₂ cos(ω t)) = (2 ω C₁ - ω^2 t C₂) cos(ω t) - (2 ω C₂ + ω^2 t C₁) sin(ω t)$

$m ((2 ω C₁ - ω^2 t C₂) cos(ω t) - (2 ω C₂ + ω^2 t C₁) sin(ω t)) = -k t (C₁ sin(ω t) + C₂ cos(ω t)) + F₀ sin(t ω)$

При $cos(ω t)$: $m (2 ω C₁ - ω^2 t C₂) = -k t C₂ ⇒ 2 m ω C₁ - k t C₂ = -k t C₂ ⇒ 2 m ω C₁ = 0 ⇒ C₁ = 0$

При $sin(ω t)$: $- m (2 ω C₂ + ω^2 t C₁) = -k t C₁ + F₀ ⇒ - 2 m ω C₂ - k t C₁ = -k t C₁ + F₀ ⇒ - 2 m ω C₂ = F₀ ⇒ C₂ = - \frac{F₀}{2 m ω}$

$r₁^* = - \frac{F₀ t cos(ω t)}{2 m ω}$

#### Безрезонансный случай

$r₂^*(t) = C sin(ω₂ t)$ где $ω₂ ≠ ω$

$- m ω₂^2 C sin(ω₂ t) = -k C sin(ω₂ t) + F₀ sin(ω₂ t)$

$- m ω₂^2 C = -k C + F₀ ⇒ C (k - m ω₂^2) = F₀$

$C = \frac{F₀}{k - m ω₂^2} = \frac{F₀}{m(ω^2 - ω₂^2}$

$r₂^*(t) = \frac{F₀}{m(ω^2 - ω₂^2)} sin(ω₂ t)$
"""

# ╔═╡ db5f5fae-9ee8-4846-8e8a-0185421d1166
md"""
### Задание 3.2 (Графики вынужденных колебаний)
"""

# ╔═╡ 0bd590eb-9225-4b6b-bbcf-601e0ab074c3
md"""
k = $(@bind k̇ Slider(1:.1:10, show_value=true))

m = $(@bind ṁ Slider(1:.1:10, show_value=true))

F₀ = $(@bind F₀̇ Slider(1:.1:10, show_value=true))

r₀ = $(@bind r₀̇ Slider(1:.1:10, show_value=true))

v₀ = $(@bind v₀̇ Slider(1:.1:10, show_value=true))

ω₂ = $(@bind ω₂̇ Slider(1:.1:10, show_value=true))
"""

# ╔═╡ cf150f5e-c6cd-43b6-9409-12ffd2dc565b
ω = √(k̇ / ṁ)

# ╔═╡ 4a208463-d0b5-49ce-9cb6-bc6016eef865
md"""
#### Резонансный случай
"""

# ╔═╡ 7f5bad32-53dc-4ed8-b7ef-f9a76a00b85a
begin
	t = 0:.1:100
	res = r₀̇*cos.(t * ω) .+ v₀̇*sin.(t * ω) / ω .- F₀̇ * t .* cos.(t * ω) / 2 / ṁ / ω

	plot(
		t, res,
		label=:none
	)
end

# ╔═╡ 1b736870-409e-4efe-b724-5977352378ae
md"""
#### Нерезонансный случай
"""

# ╔═╡ be39fcb8-5592-4a7d-8d29-5dc1f28d22a5
begin
	nres = r₀̇*cos.(t * ω) .+ v₀̇*sin.(t * ω) / ω .+ F₀̇ * sin.(t * ω₂̇) / (k̇ - ṁ * ω₂̇)

	plot(
		t, nres,
		label=:none
	)
end

# ╔═╡ b5ac39a0-b3f0-4c94-933c-6f707df8989e
md"""
## Задание 4. Система линейных химических реакций

### Задание 4.1 (Математическая модель)

#### Концептуальная постановка задачи
 - Объектом исследование является изменение концентрации веществ X и Y.
 - Вещество X притекает извне с постоянной скоростью k₁. Превращение вещества X в вещество Y проходит со скоростью, пропорциональной концентрации X с коэффициентом пропорциональности k₂. Вещество Y также выводится из сферы реакции со скоростью, пропорциональной концентрации Y с коэффициентом пропорциональности k₃.
 - Концентрации веществ в начальный момент времени t = 0: X₀ и Y₀.
 - Требуется определить концентрации веществ X и Y как функции от времени: X(t) и Y(t).

#### Математическая модель
Нам известно, что вещество X притекает извне с постоянной скоростью k₁.
Скорость изменения концентрации Y повышается пропорционально концентрации X (с коэффициентом k₂) и одновременно понижается пропорционально концентрации Y (с коэффициентом k₃). Таким образом, можем записать систему уравнений следующим образом:

$\begin{cases}
	\frac{ⅆ X(t)}{ⅆ t} = k₁ - k₂ X \\
	\frac{ⅆ Y(t)}{ⅆ t} = k₂ X - k₃ Y
\end{cases}$

Учтем начальные условия (концентрацию вещества в момент времени t = 0):
$\begin{cases}
	X(0) = X₀ \\
	Y(0) = Y₀
\end{cases}$

Таким образом, математическую модель исходной задачи можно
записать следующим образом:
$\begin{cases}
	\frac{ⅆ X(t)}{ⅆ t} = k₁ - k₂ X \\
	\frac{ⅆ Y(t)}{ⅆ t} = k₂ X - k₃ Y
\end{cases}, \begin{cases}
	X(0) = X₀ \\
	Y(0) = Y₀
\end{cases}$
"""

# ╔═╡ 6c18fcfd-9f3e-4ef7-a273-1592a3ca0230
md"""
### Задание 4.2 (Качественный анализ устойчивости)

k₁ = $(@bind k₁ Slider(0:.1:10, default=1, show_value=true))

k₂ = $(@bind k₂ Slider(0:.1:10, default=1, show_value=true))

k₃ = $(@bind k₃ Slider(0:.1:10, default=1, show_value=true))
"""

# ╔═╡ 7dcf43d8-f9ae-4db1-bb95-e6e12fe305f6
streamplot(
	(x, y) -> Point(
		k₁ - k₂*x,
		k₂*x - k₃*y
	), -3..3, -3..3
)

# ╔═╡ Cell order:
# ╟─24859bea-a037-11ec-1529-ede51f9bb656
# ╟─45d13863-edb9-4eb1-879d-8d4221deca50
# ╟─b0cfb272-ebe8-444a-b9e2-7ed8d3979003
# ╟─22ef3e0d-c2d7-4806-8d80-793f9ce960b2
# ╟─7abf4e62-af22-4906-9e23-4f18dd7f67ab
# ╠═03ecedd8-5700-4f61-b3e8-72c042c73c11
# ╠═5222ff14-5197-4af7-801d-bf652fa8a88a
# ╟─167ae2e2-0e79-4bb6-a12d-f11f0a5ff0a9
# ╟─3755c8f9-466a-4cad-be08-385503dbcd6b
# ╟─aed932a2-e563-475d-84a8-2f6d0d35c32e
# ╟─6cd36f0b-f2b3-4885-8a4a-16b917e8aa7e
# ╟─c46c02cd-5d1a-4f8b-aadd-50e223dde894
# ╟─8088e6e3-6a90-46d8-9a92-ec143ce55685
# ╟─7b39fe50-fc4a-4b65-b339-d3f6c039c172
# ╟─86e09784-7678-43ac-885e-6d3eedca7f3e
# ╟─f0eb22f4-b81b-46f3-8157-e99b4c6b5018
# ╠═12fc6293-a695-4c9c-8682-8892dfa081cb
# ╠═c2d02d98-3236-4f44-be6a-78524cc44e6c
# ╠═971c1060-6519-4ca2-a873-6c5fccca3caa
# ╟─54a9c715-1e37-4f02-8b61-1378fa48d54b
# ╟─168e14d8-933f-445e-a9e6-905611b6e4d7
# ╠═529ec0e3-2a4b-422d-949e-ab3a6938ce9c
# ╟─cea0ab76-7faf-48c3-b6e3-421c76a63082
# ╟─118ae518-5e84-4884-af6b-14cd398aca03
# ╟─db5f5fae-9ee8-4846-8e8a-0185421d1166
# ╟─0bd590eb-9225-4b6b-bbcf-601e0ab074c3
# ╠═cf150f5e-c6cd-43b6-9409-12ffd2dc565b
# ╟─4a208463-d0b5-49ce-9cb6-bc6016eef865
# ╟─7f5bad32-53dc-4ed8-b7ef-f9a76a00b85a
# ╟─1b736870-409e-4efe-b724-5977352378ae
# ╟─be39fcb8-5592-4a7d-8d29-5dc1f28d22a5
# ╟─b5ac39a0-b3f0-4c94-933c-6f707df8989e
# ╟─6c18fcfd-9f3e-4ef7-a273-1592a3ca0230
# ╟─7dcf43d8-f9ae-4db1-bb95-e6e12fe305f6
