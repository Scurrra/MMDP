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

	using DifferentialEquations: ODEProblem, solve
	using Plots: plot
	using CairoMakie: Point, streamplot, (..)
	using SymPy: @syms, Differential, dsolve, diff

	using PlutoUI
	TableOfContents()
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

# ╔═╡ 94873f6e-e5a9-4e42-b693-30bef5901c75
md"""
### Задание 1.2 (график аналитического решения)
"""

# ╔═╡ acc363ae-780d-417f-adab-c1bac579038a
md"""
α₀ = $(@bind α₀ Slider(-5:0.1:5, default=0, show_value=true))
β₀ = $(@bind β₀ Slider(-5:0.1:5, default=0, show_value=true))
"""

# ╔═╡ b9bc3371-1fdd-456f-839c-0a4d490f7d4f
t0, N0 = 1, 10

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

# ╔═╡ 5a6195f9-d27a-4d3d-9df5-f832b9321874
md"""
### Задание 1.3 (качественный анализ модели по решению)

При α₀ < β₀ наблюдается экспоненциальное убывание, а при α₀ > β₀ -- экспоненциальное возрастание. При α₀ ⩵ β₀ N(t₀) ≡ N₀ является положением равновесия.
"""

# ╔═╡ ac096f46-799d-4f1a-88cc-793f076bc425
md"""
### Задание 1.4 (качественный анализ модели по фазовому портрету)
"""

# ╔═╡ f98a4b52-bb8c-4a9b-9081-1924882013d4
streamplot(
	(x, y) -> Point(
		# (α₀ - β₀) * x,
		(1 - 2) * x,
		0
	), 0..10, -0.1..0.1
)

# ╔═╡ 61b9e905-da33-4802-a223-3e8aa42d5165
md"""
## Задание 2. Логистическая модель или модель Ферхюльста
### Задание 2.1 (аналитическое решение)

$\frac{∂ N(t)}{∂ t} = N(t) k \left( 1 - \frac{N(t)}{Nₚ} \right)$
"""

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

# ╔═╡ ece00c65-098b-4862-a55e-fb5ed1301466
streamplot(
	(x, y) -> Point(
		#x * kₚ * (1 - x / Np),
		x * 0.1 * (1 - x / 12*10^9),
		0
	), 0..10, -0.1..0.1
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
