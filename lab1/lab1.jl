### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ c5717b80-88bf-11ec-18bd-258b0006abe5
begin
	import Pkg
	Pkg.activate()

	using RecipesBase
	using Plots: plot, @recipe
	using CairoMakie: Point, streamplot, (..)
	
	using DifferentialEquations: ODEProblem, solve
	using SymPy: @syms, Differential, dsolve, diff
	using LsqFit
	
	using CSV, DataFrames, OffsetArrays
	
	import Base: show, getindex

	using PlutoUI
	TableOfContents()
end

# ╔═╡ a27fb1c6-9a0e-47e7-a264-09a820aa9598
@syms t, α(), β(), N(), t₀, N₀

# ╔═╡ 61b99f88-623a-4fdf-beb5-03f0ec1bcc9f
∂ₜ = Differential(t)

# ╔═╡ bc3ba429-f8c3-415c-b557-096a9b2ae27d
eq₁ = ∂ₜ(N) ~ (α(t) - β(t)) * N(t)

# ╔═╡ 7780cacc-c747-4f5a-b80c-6c3b5cdd7d38
dsolve(eq₁)

# ╔═╡ f5d39fcd-0558-4221-a711-8d30e0da737c
dsolve(eq₁, ics=Dict(N(t₀) => N₀))

# ╔═╡ 387c8fa8-65f8-48e3-b402-bd9ec0526e4d
"""
	eq₁cauchy(; α₀=α, β₀=β, N₀=N₀, t₀=t₀)

Linear Maltus equation Cauchy problem solution
"""
function eq₁cauchy(; α₀=α, β₀=β, N₀=N₀, t₀=t₀)
	dsolve(
		∂ₜ(N) ~ (α₀(t) - β₀(t)) * N(t), 
		ics=Dict(N(t₀) => N₀)
	)
end

# ╔═╡ 69495793-fbbd-4abb-820d-b0d94025138d
eq₁cauchy(N₀ = 10, t₀=1)

# ╔═╡ b9bc3371-1fdd-456f-839c-0a4d490f7d4f
t0, N0 = 1, 10

# ╔═╡ f98a4b52-bb8c-4a9b-9081-1924882013d4
streamplot(
	(x, y) -> Point(
		# (α₀ - β₀) * x,
		(1 - 2) * x,
		0
	), 0..10, -0.1..0.1
)

# ╔═╡ 4ee66613-3f1e-4d97-a793-b4c98678ce06
@syms Nₚ, k

# ╔═╡ 3c9fe6f1-a98b-4361-9e78-34eed8e67cd8
eq₂ = ∂ₜ(N) ~ N(t) * k * (1 - N(t) / Nₚ)

# ╔═╡ dab7c447-aeea-4fc5-a0e9-932648379f6f
dsolve(eq₂)

# ╔═╡ 66b1ccb8-ae23-4dda-aa59-eca861abaf86
dsolve(eq₂, ics=Dict(N(t₀) => N₀))

# ╔═╡ b1cfa4e2-6346-4420-9e46-fc3f24bd6099
diff(
	dsolve(eq₂, ics=Dict(N(t₀) => N₀)),
	t
)

# ╔═╡ 4f421093-c724-43c9-a149-a1b04538ac5d
"""
	eq₂cauchy(; k=k, Nₚ=Nₚ, N₀=N₀, t₀=t₀)

Verhulst equation Cauchy problem solution
"""
function eq₂cauchy(; k=k, Nₚ=Nₚ, N₀=N₀, t₀=t₀)
	dsolve(
		∂ₜ(N) ~ N(t) * k * (1 - N(t) / Nₚ), 
		ics=Dict(N(t₀) => N₀)
	)
end

# ╔═╡ 56c29f63-0b6e-4f0f-80c2-239d22756068
eq₂cauchy(t₀=1, N₀=10^9)

# ╔═╡ ece00c65-098b-4862-a55e-fb5ed1301466
streamplot(
	(x, y) -> Point(
		#x * kₚ * (1 - x / Np),
		x * 0.1 * (1 - x / 12*10^9),
		0
	), 0..10, -0.1..0.1
)

# ╔═╡ 899a2de6-a92e-403e-85c9-e98dca4688d1
@syms α₀₀, β₀₀

# ╔═╡ a4c52dc1-51fe-420c-94d0-979c19c507df
eq₃ = ∂ₜ(N) ~ (α₀₀ * N(t) - β₀₀) * N(t)

# ╔═╡ b885f5ec-06e8-4278-a115-c80fc50942d5
dsolve(eq₃)

# ╔═╡ 436712ac-e656-47f7-b486-ae3d6789b2d6
dsolve(eq₃, ics=Dict(N(t₀) => N₀))

# ╔═╡ 061cfd80-4783-4399-8814-ec9dc3b00ae7
"""
	eq₃cauchy(; α₀=α, β₀=β, N₀=N₀, t₀=t₀)

Non-Linear Maltus equation Cauchy problem solution
"""
function eq₃cauchy(; α₀=α, β₀=β, N₀=N₀, t₀=t₀)
	dsolve(
		∂ₜ(N) ~ (α₀(t) * N(t) - β₀(t)) * N(t), 
		ics=Dict(N(t₀) => N₀)
	)
end

# ╔═╡ 5a39345e-b1dc-4cf5-89e3-0d408eb6e87c
eq₃cauchy(N₀ = 10, t₀=1)

# ╔═╡ 929d3561-f660-4f70-9ab9-a55dc1ecfcf2
begin
	struct CountryKnown
		name::String
		year::Int
		year_end::Int
	
		α₀
		β₀
		N₀

		solution
	
		CountryKnown(α₀, β₀, N₀; name::String, year::Int=2022, year_end::Int=2122) = new(
			name, year, year_end, α₀, β₀, N₀,
			solve(
				ODEProblem(
					(N, p, t) -> ((α₀/N₀) * N - β₀) * N,
					N₀,
					(year, year_end)
				)
			)
		)
	end
	
	#Base.show(io::IO, c::Country) = print(name, 
	#	eq₃cauchy(α₀=(t)->c.α₀/c.N₀, β₀=(t)->c.β₀, N₀=c.N₀, t₀=c.year)
	#)

	(c::CountryKnown)(year::Int) = c.solution(year) |> round |> Int;
end

# ╔═╡ 0499e28c-cf62-46be-a7ca-5b1d386f0564
Belarus = CountryKnown(0.011874, 0.013346, 9450233, name="Belarus", year=2017);

# ╔═╡ 586e57ce-af4b-4629-9221-c8745e2d720e
plot(
	Belarus.solution,
	label="Population",
	title=Belarus.name
)

# ╔═╡ c15f1c0f-e4cc-4d71-af6f-ca0b64353895
population = CSV.read("population.csv", DataFrame);

# ╔═╡ 3abda224-3fbd-420c-9fa7-eeffd33c7637
@bind country Select(names(population), default="Belarus")

# ╔═╡ a10590bb-1d7c-4488-82bf-6e2c2f74da7b
begin
	mutable struct CountryPopulation
		name::String
		year_start::Int
		year_end::Int
		data::OffsetArray

		LM
		V
		NM
		
		CountryPopulation(name, data; year_start::Int=1960, year_end::Int=2020) = new(
			name,
			year_start,
			year_end,
			OffsetArray(data, year_start:year_end),
			@NamedTuple{αβ₀::Float64}(curve_fit(
				(t, p) -> data[1] * exp.(p[1] * (t .- year_start)), 
				year_start:year_end |> collect, 
				data |> collect, 
				[.0]
			).param |> Tuple),
			# doesn`t fit
			@NamedTuple{k::Float64, Nₚ::Int}(curve_fit((t, p) -> data[1] * p[2] * exp.(p[1] * (t .- year_start)) ./ (data[1] - p[2]) ./ (data[1] * exp.(p[1] * (t .- year_start)) / (data[1] - p[2]) .- 1), 
				year_start:year_end |> collect, 
				data |> collect, 
				[1., data[1]]
			).param |> Tuple),
			# fit`s into a straight horisontal line
			@NamedTuple{α₀::Float64, β₀::Float64}(curve_fit(
			(t, p) -> p[2] ./ (p[1] .- exp.( p[2] * (t .- year_start .+ log( Complex(p[1] - p[2]/data[1]) ) / p[2]) )), 
				year_start:year_end |> collect, 
				data |> collect, 
				[1., 1.]
			).param |> Tuple)
		)
	end

	@inline Base.getindex(p::CountryPopulation, i::Integer) = p.data[i]

	@recipe function f(c::CountryPopulation)
		
		lm = OffsetArray(
			c.data[c.year_start] * exp.(c.LM.αβ₀ * (0:(c.year_end - c.year_start))),
			c.year_start:c.year_end
		)
		
		nm = OffsetArray(
			c.NM.β₀ ./ (c.NM.α₀ .- exp.( c.NM.β₀ * ((0:(c.year_end - c.year_start)) .+ log( Complex(c.NM.α₀ - c.NM.β₀/c.data[c.year_start]) ) / c.NM.β₀) )) .|> Float64,
			c.year_start:c.year_end
		)
		
		title := c.name
		label := ["Real Population" "Linear Maltus Model" "Nonlinear Maltus Model"]
		y := [c.data, lm, nm]
		
	end
end

# ╔═╡ 4e7d36f9-dd1e-4ab6-8055-f9b124e54e68
md"""
# Лабораторная работа №1
# Математические модели динамики численности популяции одного вида
Подготовил: Боровский Илья
"""

# ╔═╡ 14a0f2c9-7a3b-4806-b988-f5ac7b5e5eeb
md"""
## Задание 1. Модель Мальтуса
### Задание 1.1 (аналитическое решение)

$\frac{∂ N(t)}{∂ t} = (α(t) + β(t))N(t)$
$\frac{d N(t)}{N} = (α(t) + β(t)) d t$
$ln(N(t)) = \int (α(t) + β(t)) d t + C₁$
$N(t) = C₁ e^{\int (α(t) + β(t)) d t}$
"""

# ╔═╡ 94873f6e-e5a9-4e42-b693-30bef5901c75
md"""
### Задание 1.2 (график аналитического решения)
"""

# ╔═╡ acc363ae-780d-417f-adab-c1bac579038a
md"""
α₀ = $(@bind α₀ Slider(-5:0.1:5, default=1, show_value=true))
β₀ = $(@bind β₀ Slider(-5:0.1:5, default=2, show_value=true))
"""

# ╔═╡ 8a1ff32e-954e-4e61-9482-1ca7952044d9
sol₁ = solve(
	ODEProblem(
		(N, p, t) -> (α₀ - β₀) * N,
		N0,
		(t0, 1000)
	)
);

# ╔═╡ 863770d3-357d-489c-9833-8b14a4813e6a
plot(
	sol₁,
	title="Maltus Linear Equation",
	label="Solution"
)

# ╔═╡ 0c62ffde-5772-4ada-8e6b-9f4af75b1255
begin 
	N0₃ = Dict(
		"N₀" => N0,
		"N₀ < Nₖₚ" => β₀/α₀-0.0001,
		"N₀ ⩵ Nₖₚ" => β₀/α₀,
		"N₀ > Nₖₚ" => β₀/α₀+0.0001,
	)
	@bind key₃ Radio(["N₀", "N₀ < Nₖₚ", "N₀ ⩵ Nₖₚ", "N₀ > Nₖₚ"], default="N₀")
end

# ╔═╡ 2a353848-15e8-4bb6-8944-66aeacd059d3
sol₃ = solve(
	ODEProblem(
		(N, p, t) -> (α₀ * N - β₀) * N,
		N0₃[key₃],
		(t0, 1000)
	)
);

# ╔═╡ 3af61166-36c6-47bf-9c7f-e3a0a84195a9
plot(
	sol₃,
	title="Maltus Non-Linear Equation",
	label="Solution"
)

# ╔═╡ 325cc855-878e-4375-ae94-71c368771a41
streamplot(
	(x, y) -> Point(
		(α₀ * x - β₀) * x,
		# (1 * x - 2) * x,
		0
	), 0..10, -0.1..0.1
)

# ╔═╡ 5a6195f9-d27a-4d3d-9df5-f832b9321874
md"""
### Задание 1.3 (качественный анализ модели по решению)

При α₀ < β₀ наблюдается экспоненциальное убывание, а при α₀ > β₀ -- экспоненциальное возрастание. При α₀ ⩵ β₀ N(t₀) ≡ N₀ является положением равновесия.
"""

# ╔═╡ ac096f46-799d-4f1a-88cc-793f076bc425
md"""
### Задание 1.4 (качественный анализ модели по фазовому портрету)
"""

# ╔═╡ 61b9e905-da33-4802-a223-3e8aa42d5165
md"""
## Задание 2. Логистическая модель или модель Ферхюльста
### Задание 2.1 (аналитическое решение)

$\frac{∂ N(t)}{∂ t} = N(t) k \left( 1 - \frac{N(t)}{Nₚ} \right)$
"""

# ╔═╡ 48be88ff-282a-4ada-9f07-9e6a37588127
md"""
### Задание 2.2 (график аналитического решения)
"""

# ╔═╡ 2c52968b-3886-4b6b-9aa3-9419bdef8273
md"""
Nₚ = $(@bind Np Slider(1000000000:1000000000:50000000000, default=12000000000, show_value=true))

N₀ = $(@bind N0ₚ Slider(1000000000:1000000000:50000000000, default=12000000000, show_value=true))

k = $(@bind kₚ Slider(0:.1:10, default=1, show_value=true))
"""

# ╔═╡ db2c3767-5030-4756-9614-7f506b147709
sol₂ = solve(
	ODEProblem(
		(N, p, t) -> N * kₚ * (1 - N / Np),
		N0ₚ,
		(t0, 1000)
	)
);

# ╔═╡ 6d96f216-407f-40cb-88bd-6a37f06a9173
plot(
	sol₂,
	title="Verhulst Model",
	label="Solution"
)

# ╔═╡ 6a06dfcc-a3cf-4659-9c5b-507bc3813c47
md"""
### Задание 2.3 (качественный анализ модели по решению)

При Nₚ < N₀ наблюдается экспоненциальное убывание до Nₚ, а при Nₚ > N₀ наблюдается экспоненциальный рост до Nₚ, а далее -- график константный. N(t₀) ≡ N₀ при N₀ ⩵ Nₚ является точкой равновесия.
"""

# ╔═╡ 3f1005b7-2204-4f3e-97a9-cf81a6c881b3
md"""
### Задание 2.4 (качественный анализ модели по фазовому портрету)
"""

# ╔═╡ ad440370-79d3-4510-8630-5c4d77e0615f
md"""
## Задание 3. Нелинейный аналог модели Мальтуса
### Задание 3.1 (аналитическое решение)

$\frac{∂ N(t)}{∂ t} = (α(t)*N(t) + β(t))N(t)$
"""

# ╔═╡ d2851b4c-36e7-4331-aad7-1237abe3ce45
md"""
### Задание 3.2 (график аналитического решения)
"""

# ╔═╡ 69647171-e695-4a71-8990-0eab57fc436a
md"""
### Задание 3.3 (качественный анализ модели по решению)

При N₀ < Nₖₚ наблюдается экспоненциальное убывание, а при N₀ > Nₖₚ -- экспоненциальный рост. При N₀ ⩵ Nₖₚ численность населения не изменяется.
"""

# ╔═╡ 60195422-7e30-4133-99e0-f9bba7f1b203
md"""
### Задание 3.4 (качественный анализ модели по фазовому портрету)
"""

# ╔═╡ 44950393-ac60-479a-b589-f98bb0bb8cb7
md"""
## Задание 4. Предсказание численности населения по заданным параметрам модели
"""

# ╔═╡ 4dcf2096-7657-4ce9-af1b-94cc8fa23659
md"""
$(@bind year_Belarus Slider(Belarus.year:Belarus.year_end, show_value=true))
"""

# ╔═╡ 31ab9e1d-2c11-4f61-81c2-444ef82d48a9
md"""
Belarus population in $(year_Belarus): $(Belarus(year_Belarus))
"""

# ╔═╡ 8af93dde-6486-4639-8cc9-862f05876bd1
md"""
## Задание 5-6. Изменение численности населения страны и мира. Определение параметров модели по реальным данным
"""

# ╔═╡ 760fba49-e682-4fbf-8e04-54d92d1546aa
begin
	dsolve(		∂ₜ(N) ~ (α₀₀ * N(t) - β₀₀) * N(t)
		, ics=Dict(N(t₀) => N₀)
	)
end

# ╔═╡ 3e71a374-9df6-4db1-877e-baf371cc4fa9
country_population = CountryPopulation(country, population[!, country])

# ╔═╡ c2a61651-9101-4d62-9282-e472c7f1bdb9
plot(
	country_population
)

# ╔═╡ Cell order:
# ╠═c5717b80-88bf-11ec-18bd-258b0006abe5
# ╟─4e7d36f9-dd1e-4ab6-8055-f9b124e54e68
# ╟─14a0f2c9-7a3b-4806-b988-f5ac7b5e5eeb
# ╠═a27fb1c6-9a0e-47e7-a264-09a820aa9598
# ╠═61b99f88-623a-4fdf-beb5-03f0ec1bcc9f
# ╠═bc3ba429-f8c3-415c-b557-096a9b2ae27d
# ╠═7780cacc-c747-4f5a-b80c-6c3b5cdd7d38
# ╠═f5d39fcd-0558-4221-a711-8d30e0da737c
# ╟─387c8fa8-65f8-48e3-b402-bd9ec0526e4d
# ╠═69495793-fbbd-4abb-820d-b0d94025138d
# ╟─94873f6e-e5a9-4e42-b693-30bef5901c75
# ╟─acc363ae-780d-417f-adab-c1bac579038a
# ╠═b9bc3371-1fdd-456f-839c-0a4d490f7d4f
# ╠═8a1ff32e-954e-4e61-9482-1ca7952044d9
# ╟─863770d3-357d-489c-9833-8b14a4813e6a
# ╟─5a6195f9-d27a-4d3d-9df5-f832b9321874
# ╟─ac096f46-799d-4f1a-88cc-793f076bc425
# ╠═f98a4b52-bb8c-4a9b-9081-1924882013d4
# ╟─61b9e905-da33-4802-a223-3e8aa42d5165
# ╠═4ee66613-3f1e-4d97-a793-b4c98678ce06
# ╠═3c9fe6f1-a98b-4361-9e78-34eed8e67cd8
# ╠═dab7c447-aeea-4fc5-a0e9-932648379f6f
# ╠═66b1ccb8-ae23-4dda-aa59-eca861abaf86
# ╠═b1cfa4e2-6346-4420-9e46-fc3f24bd6099
# ╟─4f421093-c724-43c9-a149-a1b04538ac5d
# ╠═56c29f63-0b6e-4f0f-80c2-239d22756068
# ╟─48be88ff-282a-4ada-9f07-9e6a37588127
# ╟─2c52968b-3886-4b6b-9aa3-9419bdef8273
# ╠═db2c3767-5030-4756-9614-7f506b147709
# ╟─6d96f216-407f-40cb-88bd-6a37f06a9173
# ╟─6a06dfcc-a3cf-4659-9c5b-507bc3813c47
# ╟─3f1005b7-2204-4f3e-97a9-cf81a6c881b3
# ╠═ece00c65-098b-4862-a55e-fb5ed1301466
# ╟─ad440370-79d3-4510-8630-5c4d77e0615f
# ╠═899a2de6-a92e-403e-85c9-e98dca4688d1
# ╠═a4c52dc1-51fe-420c-94d0-979c19c507df
# ╠═b885f5ec-06e8-4278-a115-c80fc50942d5
# ╠═436712ac-e656-47f7-b486-ae3d6789b2d6
# ╟─061cfd80-4783-4399-8814-ec9dc3b00ae7
# ╠═5a39345e-b1dc-4cf5-89e3-0d408eb6e87c
# ╟─d2851b4c-36e7-4331-aad7-1237abe3ce45
# ╠═0c62ffde-5772-4ada-8e6b-9f4af75b1255
# ╠═2a353848-15e8-4bb6-8944-66aeacd059d3
# ╟─3af61166-36c6-47bf-9c7f-e3a0a84195a9
# ╟─69647171-e695-4a71-8990-0eab57fc436a
# ╟─60195422-7e30-4133-99e0-f9bba7f1b203
# ╠═325cc855-878e-4375-ae94-71c368771a41
# ╟─44950393-ac60-479a-b589-f98bb0bb8cb7
# ╠═929d3561-f660-4f70-9ab9-a55dc1ecfcf2
# ╠═0499e28c-cf62-46be-a7ca-5b1d386f0564
# ╟─4dcf2096-7657-4ce9-af1b-94cc8fa23659
# ╟─31ab9e1d-2c11-4f61-81c2-444ef82d48a9
# ╟─586e57ce-af4b-4629-9221-c8745e2d720e
# ╟─8af93dde-6486-4639-8cc9-862f05876bd1
# ╠═c15f1c0f-e4cc-4d71-af6f-ca0b64353895
# ╟─3abda224-3fbd-420c-9fa7-eeffd33c7637
# ╠═a10590bb-1d7c-4488-82bf-6e2c2f74da7b
# ╠═760fba49-e682-4fbf-8e04-54d92d1546aa
# ╠═3e71a374-9df6-4db1-877e-baf371cc4fa9
# ╠═c2a61651-9101-4d62-9282-e472c7f1bdb9
