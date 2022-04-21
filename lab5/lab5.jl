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
			
			new{T, T₀, Tₙ}(solution)
		end
	end

	(dde::DDE{T, T₀, Tₙ})(t) where {T, T₀, Tₙ} = t in RightSemiInterval{T₀-T, Tₙ}() ? (filter(sol -> t ∈ sol.first, dde.solution) |> values |> collect)[1](t) : throw("`t` is out of boundary")

	@recipe function f(dde::DDE{T, T₀, Tₙ}; color::Symbol=:blue) where {T, T₀, Tₙ}
		xlims := (T₀, Tₙ)
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

### Задание 2.1 (Аналитическое решение)


### Задание 2.2 (Численное решение)

T = $(@bind T₂₂ Slider(.1:.1:10; show_value=true))
"""

# ╔═╡ 70bbda40-29b9-419d-b0c7-853f71ec7c9d
sol = DDE{T₂₂, 0., 10.}((t, x, ϕ) -> -π*ϕ(t)/2/T₂₂, (t) -> cos(π*t/2/T₂₂));

# ╔═╡ ee310025-efb9-460b-8fd1-3d1700d439a3
sol(0.5)

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
	label="DifferentialEquations.jl"
)

# ╔═╡ 069da2ba-1b1a-43bb-ab59-6cb6c6b754c6


# ╔═╡ Cell order:
# ╠═73b63744-c182-11ec-0094-6991e9545112
# ╟─791a757d-5e14-4dbc-86de-40f11f92ad2a
# ╟─75c83886-e32c-47fa-a26b-3dd6ba9266d9
# ╠═7b40731e-9551-4d9b-ab60-0ea2ad2abdd1
# ╠═7c228e10-f2fb-4948-b35f-2bceb442bf02
# ╟─25a58eb7-aae9-45c1-a15c-e4146558447f
# ╠═70bbda40-29b9-419d-b0c7-853f71ec7c9d
# ╠═ee310025-efb9-460b-8fd1-3d1700d439a3
# ╠═4b0f0997-3f59-415d-a123-9a7331a979f1
# ╠═bf92b570-d46a-45fe-82f4-29d11344ab12
# ╠═069da2ba-1b1a-43bb-ab59-6cb6c6b754c6
