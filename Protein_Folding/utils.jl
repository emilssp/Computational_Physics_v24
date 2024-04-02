### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ b69b5ff0-f0ef-11ee-246f-a19167d7f18d
module Aminoacids
export Aminoacid, isdistance1
#define aminoacids according to the standard one letter abbreviation and assign 
#them to a number
@enum Aminoacid begin
	A = 1
	R = 2
	N = 3
	D = 4
	C = 5
	E = 6
	Q = 7
	G = 8
	H = 9
	I = 10
	L = 11
	K = 12
	M = 13
	F = 14
	P = 15
	S = 16
	T = 17
	W = 18
	Y = 19
	V = 20
end	

function isdistance1(pos1, pos2)
    return norm(pos1 - pos2) <= 1
end

end

# ╔═╡ 9439b7f1-005e-4522-b7a3-f64684ad7303


# ╔═╡ be8c5275-d6ec-484a-8041-fa17088ed685
module Monomers
	import ..Aminoacids: Aminoacid
	export Monomer
	mutable struct Monomer
		id::Int64
		kind::Aminoacid #type of the monomer
		nearest_neighbors::Vector{Monomer}
		pos::Vector{Int64}
		available_moves::Vector{Vector{Int64}}
		function Monomer(id::Int64, kind::Aminoacid ,pos::Vector{Int64})
			new(id, kind, Vector{Monomer}(), pos)
		end
	end
end

# ╔═╡ Cell order:
# ╠═b69b5ff0-f0ef-11ee-246f-a19167d7f18d
# ╠═9439b7f1-005e-4522-b7a3-f64684ad7303
# ╠═be8c5275-d6ec-484a-8041-fa17088ed685
