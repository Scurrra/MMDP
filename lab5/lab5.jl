### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 73b63744-c182-11ec-0094-6991e9545112
begin
	import Pkg;
	Pkg.activate();

	using DataFrames;

	using RecipesBase;
	using Plots; gr();
	#using Plots: plot, plot!, @recipe, @animate, @gif, scatter;
	using CairoMakie: Point, Figure, Axis, scatter!, streamplot, streamplot!, vlines!, (..);
	using LinearAlgebra
		
	using DifferentialEquations: ODEProblem, DDEProblem, MethodOfSteps, Tsit5, solve;

	using PlutoUI;
	TableOfContents()
end

# ╔═╡ 791a757d-5e14-4dbc-86de-40f11f92ad2a
md"""
# Лабораторная работа №5
# Математические модели с запаздыванием
Подготовил: Боровский Илья
"""

# ╔═╡ 75c83886-e32c-47fa-a26b-3dd6ba9266d9
md"""
## Задание 1. Метод последовательного интегрирования (метод шагов)
"""

# ╔═╡ 7b40731e-9551-4d9b-ab60-0ea2ad2abdd1
begin
	import Base.in, Base.∈, Base.show

	abstract type AbstractInterval end

	struct Interval{BEGIN, END} <: AbstractInterval
		Interval{BEGIN, END}() where {BEGIN, END} = BEGIN < END ? new{BEGIN, END}() : throw("Invalid interval boundary")
	end
	Base.show(io::IO, int::Interval{BEGIN, END}) where {BEGIN, END} = print(io, "($BEGIN, $END)")
	Base.in(item, int::Interval{BEGIN, END}) where {BEGIN, END} = item > BEGIN && item < END
	
	struct RightSemiInterval{BEGIN, END} <: AbstractInterval
		RightSemiInterval{BEGIN, END}() where {BEGIN, END} = BEGIN <= END ? new{BEGIN, END}() : throw("Invalid interval boundary")
	end
	Base.show(io::IO, int::RightSemiInterval{BEGIN, END}) where {BEGIN, END} = print(io, "($BEGIN, $END]")
	Base.in(item, int::RightSemiInterval{BEGIN, END}) where {BEGIN, END} = item > BEGIN && item <= END

	struct LeftSemiInterval{BEGIN, END} <: AbstractInterval
		LeftSemiInterval{BEGIN, END}() where {BEGIN, END} = BEGIN <= END ? new{BEGIN, END}() : throw("Invalid interval boundary")
	end
	Base.show(io::IO, int::LeftSemiInterval{BEGIN, END}) where {BEGIN, END} = print(io, "[$BEGIN, $END)")
	Base.in(item, int::LeftSemiInterval{BEGIN, END}) where {BEGIN, END} = item >= BEGIN && item < END
end;

# ╔═╡ 7c228e10-f2fb-4948-b35f-2bceb442bf02
begin	
	struct DDE{T, T₀, Tₙ}
		solution
		ϕ::Function

		function DDE{T, T₀, Tₙ}(f::Function, ϕ::Function) where {T, T₀, Tₙ}
			solution = Dict{AbstractInterval, Any}(
				RightSemiInterval{T₀-T, T₀}() => ϕ
			)

			k::Int = 0
			while k*T < Tₙ
				solution[RightSemiInterval{T₀+k*T, T₀+(k+1)*T}()] = ODEProblem(
					(x, p, t) -> f(
						t, 
						x, 
						solution[RightSemiInterval{T₀+(k-1)*T, T₀+k*T}()]
					),
					solution[RightSemiInterval{T₀+(k-1)*T, T₀+k*T}()](T₀+k*T),		
					(T₀+k*T, T₀+(k+1)*T)
				) |> solve
				
				k += 1
			end

			delete!(solution, RightSemiInterval{T₀-T, T₀}())
			
			new{T, T₀, Tₙ}(solution, ϕ)
		end
	end

	(dde::DDE{T, T₀, Tₙ})(t) where {T, T₀, Tₙ} = 
		t in RightSemiInterval{T₀, Tₙ}() ? 
			(filter(sol -> t ∈ sol.first, dde.solution) |> values |> collect)[1](t) : 
			(t in RightSemiInterval{T₀-T, T₀}() ?
				dde.ϕ(t) : throw("`t` is out of boundary"))

	@recipe function f(dde::DDE{T, T₀, Tₙ}; color::Symbol=:blue) where {T, T₀, Tₙ}
		@series begin
			xlims := (T₀-T, T₀)
			color := color
			label := :none
			dde.ϕ
		end
		
		xlims := (T₀-T, Tₙ)
		title := "DDE Solution"
		for sol ∈ dde.solution |> values
			@series begin
				color := color
				label := :none
				sol
			end
		end
	end
end

# ╔═╡ 25a58eb7-aae9-45c1-a15c-e4146558447f
md"""
## Задание 2. Линейное дифференциальное уравнение с запаздывающим аргументом

$x'(t) = -\frac{π x(t-T)}{2T}$

### Задание 2.1 (Аналитическое решение)

$x(t) = ℯ^{λ t} ⇒ x'(t) = λ ℯ^{λ t} = -\frac{π ℯ^{λ (t-T)}}{2T} ⇒ λ = -\frac{π ℯ^{-λ T}}{2T}$

Положим $λ = α + i β$:

$α + i β = -\frac{π ℯ^{-(α + i β) T}}{2T} = -\frac{π ℯ^{-α T}}{2T} \left[ cos(-β t) + i sin(-β t) \right] \Rightarrow$

$\begin{cases}
	α = -\frac{π ℯ^{-α T}}{2T} cos(-β t) \\
	β = -\frac{π ℯ^{-α T}}{2T} sin(-β t)
\end{cases}$

Если $β$ -- решение системы, то $-β$ -- тоже решение системы ⇔ $α^2 + β^2 = \frac{π^2 ℯ^{-2αT}}{4 T^2} ⇒ β = ± √ \left( \frac{π^2 ℯ^{-2αT}}{4 T^2} - α^2 \right)$

Положим $α = 0 ⇒ β = ± \frac{π}{2 T}, x(t) = C Re \left( ℯ^{i β t} \right) = C cos(β t) = C cos \left( \frac{π t}{2 T} \right)$

### Задание 2.2 (Численное решение)

T = $(@bind T₂₂ Slider(.1:.1:20; show_value=true, default=6))
"""

# ╔═╡ 70bbda40-29b9-419d-b0c7-853f71ec7c9d
sol = DDE{T₂₂, 0., 10.}((t, x, ϕ) -> -π*ϕ(t)/2/T₂₂, (t) -> cos(π*t/2/T₂₂));

# ╔═╡ ee310025-efb9-460b-8fd1-3d1700d439a3
sol(-0.5)

# ╔═╡ 4b0f0997-3f59-415d-a123-9a7331a979f1
plot(
	sol;
	color=:green
)

# ╔═╡ bf92b570-d46a-45fe-82f4-29d11344ab12
plot!(
	solve(
		DDEProblem(
			(x, ϕ, p, t) -> -π*ϕ(p, t)/2/T₂₂,
			(p, t) -> cos(π*t/2/T₂₂),
			(0, 10)
		),
		MethodOfSteps(Tsit5())
	);
	label="DifferentialEquations.jl",
	color=:red,
	xlims = (-T₂₂, 10)
)

# ╔═╡ d4757f2d-f36a-4001-91f4-175b879efa91
md"""
## Задание 3. Начальная задача для уравнения Хатчинсона

### Математическая модель

$\begin{cases}
	N`(t) = N(t)*(1 - N(t-T)), & t>0 \\
	N(t) = ϕ(t), & -T \le t \le 0
\end{cases}$

где $ϕ(t) = 1+\frac{t}{T}, t ∈ [-T, 0], T=const>0$

### Задание 3.1 (Численное решение)

T = $(@bind T₃₁ Slider(0:(π/16):π; show_value=true, default=(π/2)))
"""

# ╔═╡ 2c3c9e88-f17a-4499-ab54-9cd21f75088a
Hutchinson = DDE{T₃₁, 0., 10π}(
	(t, x, ϕ) -> x*(1-ϕ(t)), 
	(t) -> 1 + t/T₃₁
);

# ╔═╡ e97d0019-e581-46da-aff6-9ef9abd3c405
plot(
	Hutchinson;
	color=:green
)

# ╔═╡ 2a138947-a217-4f22-9a3d-0ac6d1979992
plot!(
	solve(
		DDEProblem(
			(x, ϕ, p, t) -> x*(1-ϕ(p, t)),
			(p, t) -> 1 + t/T₃₁,
			(0, 10π)
		),
		MethodOfSteps(Tsit5())
	);
	label="DifferentialEquations.jl",
	color=:red,
	xlims = (-T₃₁, 10π)
)

# ╔═╡ f2e9448a-9056-409e-84f8-eff637cc0a02
md"""
### Задание 4. Модель регуляции концентрации клеток крови

$\begin{cases}
	с`(t) = \frac{λ \; c(t-T)}{1 + c(t-T)^m} - g \; c(t), & t>0 \\
	с(t) = ϕ(t), & -T \le t \le 0
\end{cases}$

### Задание 4.1 (Численное решение)

T = $(@bind T₄₁ Slider(1:1:20; show_value=true, default=6))

λ = $(@bind λ Slider(.1:.1:1; show_value=true, default=.2))

m = $(@bind m Slider(1:1:20; show_value=true, default=10))

g = $(@bind g Slider(.1:.1:1; show_value=true, default=.1))
"""

# ╔═╡ 2af66f38-50b2-4534-9242-e931f6ced415
blood = DDE{T₄₁, 0, 600}(
	(t, x, ϕ) -> λ * ϕ(t) / (1+ϕ(t)^m) - g*x, 
	(t) -> .1
);

# ╔═╡ 4bd8ef51-5e5f-4d96-b64a-11da02d56a5f
plot(
	blood;
	color=:green
)

# ╔═╡ 5ce80336-a2ef-49fa-b263-3ce540f4a961
plot!(
	solve(
		DDEProblem(
			(x, ϕ, p, t) -> λ * ϕ(p, t) / (1+ϕ(p, t)^m) - g*x,
			(p, t) -> .1,
			(0, 600)
		),
		MethodOfSteps(Tsit5())
	);
	label="DifferentialEquations.jl",
	color=:red,
	xlims = (-T₄₁, 600)
)

# ╔═╡ cc68a013-6d37-44dc-87e5-1fb95df8bdbd
md"""
### Задание 4.2 (Фазовый портрет)
"""

# ╔═╡ 8975b924-05df-4029-b765-bd0f6709afab
function streamplot_delayed(dde::DDE{T, T₀, Tₙ}; params...) where {T, T₀, Tₙ}
	t = ((T/100):(T/100):T |> collect) .+ T₀
	
	plots = [plot(
		dde.(t),
		dde.ϕ.(t .- T);
		title="t ∈ (T₀, T₀+T]",
		label=:none,
		params...
	)]

	
	k = 1
	while k*T < Tₙ-T
		push!(
			plots,
			plot(
				dde.(t .+ k*T),
				dde.(t .+ (k-1)*T);
				title="t ∈ (T₀+$k*T, T₀+$(k+1)*T]",
				label=:none,
				params...
			)
		)
		k+=1
	end

	return plots
end

# ╔═╡ 991f8ee9-db98-4b7a-9adb-6c0ab6c097b0
blood_streamplots = streamplot_delayed(
	DDE{2, 0, 600}(
		(t, x, ϕ) -> 2 * ϕ(t) / (1+ϕ(t)^20) - 1*x, 
		(t) -> .9
	); color=:green,
	xlims=(0, 1.5),
	ylims=(0, 1.5)
);

# ╔═╡ 48df330d-39ca-4992-8eb5-e155f4789e3f
md"""
k = $(@bind k Slider(0:(length(blood_streamplots)-1); show_value=true))
"""

# ╔═╡ e9c8a21b-fbfb-4182-8792-820d7875ac36
blood_streamplots[k+1]

# ╔═╡ Cell order:
# ╟─73b63744-c182-11ec-0094-6991e9545112
# ╟─791a757d-5e14-4dbc-86de-40f11f92ad2a
# ╟─75c83886-e32c-47fa-a26b-3dd6ba9266d9
# ╠═7b40731e-9551-4d9b-ab60-0ea2ad2abdd1
# ╠═7c228e10-f2fb-4948-b35f-2bceb442bf02
# ╟─25a58eb7-aae9-45c1-a15c-e4146558447f
# ╠═70bbda40-29b9-419d-b0c7-853f71ec7c9d
# ╠═ee310025-efb9-460b-8fd1-3d1700d439a3
# ╠═4b0f0997-3f59-415d-a123-9a7331a979f1
# ╠═bf92b570-d46a-45fe-82f4-29d11344ab12
# ╟─d4757f2d-f36a-4001-91f4-175b879efa91
# ╠═2c3c9e88-f17a-4499-ab54-9cd21f75088a
# ╠═e97d0019-e581-46da-aff6-9ef9abd3c405
# ╠═2a138947-a217-4f22-9a3d-0ac6d1979992
# ╟─f2e9448a-9056-409e-84f8-eff637cc0a02
# ╠═2af66f38-50b2-4534-9242-e931f6ced415
# ╠═4bd8ef51-5e5f-4d96-b64a-11da02d56a5f
# ╠═5ce80336-a2ef-49fa-b263-3ce540f4a961
# ╟─cc68a013-6d37-44dc-87e5-1fb95df8bdbd
# ╠═8975b924-05df-4029-b765-bd0f6709afab
# ╠═991f8ee9-db98-4b7a-9adb-6c0ab6c097b0
# ╟─48df330d-39ca-4992-8eb5-e155f4789e3f
# ╠═e9c8a21b-fbfb-4182-8792-820d7875ac36
