# Independent enumeration of the extreme rays of the vertex cones by the double description
# method (a separate implementation from checkers/dd.c, and from cdd and lrs).
# The cone of a vertex J with minimal vectors k_i is { Y : Y[k_i] >= 0 } in the coordinates
# y_ab (a <= b); the rays found are compared with the RAYS of the certificate.
# Usage: julia dd.jl certificate.txt

const D = 21

fail(m) = (println(stderr, "FAIL: ", m); exit(1))

symc(k) = [Int64((a == b ? 1 : 2) * k[a] * k[b]) for a in 1:6 for b in a:6]

function prim(c::Vector{BigInt})
    g = reduce(gcd, c)
    g == 0 && fail("zero ray")
    c .÷ g
end

# a nonzero kernel vector of a 20 x 21 integer matrix of rank 20, by exact rational elimination
function kernel_vector(B::Matrix{Int64})
    M = Rational{BigInt}.(B)
    piv = Int[]; r = 1
    for c in 1:D
        r > 20 && break
        p = findfirst(i -> M[i, c] != 0, r:20)
        p === nothing && continue
        p += r - 1
        M[[r, p], :] = M[[p, r], :]
        M[r, :] ./= M[r, c]
        for i in 1:20
            if i != r && M[i, c] != 0
                M[i, :] .-= M[i, c] .* M[r, :]
            end
        end
        push!(piv, c); r += 1
    end
    length(piv) == 20 || fail("initial constraints are dependent")
    fr = first(setdiff(1:D, piv))
    x = zeros(Rational{BigInt}, D); x[fr] = 1
    for (i, c) in enumerate(piv); x[c] = -M[i, fr]; end
    L = lcm(denominator.(x))
    prim(numerator.(x .* L))
end

function rank_mod(rows)
    P = Int128(2)^61 - 1
    A = [mod.(Int128.(r), P) for r in rows]
    rk = 0; pivrows = Int[]
    idx = collect(1:length(A))
    for c in 1:D
        rk == length(A) && break
        p = findfirst(i -> A[i][c] != 0, rk+1:length(A))
        p === nothing && continue
        p += rk
        A[p], A[rk+1] = A[rk+1], A[p]; idx[p], idx[rk+1] = idx[rk+1], idx[p]
        rk += 1
        inv = powermod(A[rk][c], P - 2, P)
        for r in 1:length(A)
            if r != rk && A[r][c] != 0
                f = mod(A[r][c] * inv, P)
                A[r] = mod.(A[r] .- f .* A[rk], P)
            end
        end
        push!(pivrows, idx[rk])
    end
    rk, pivrows
end

function enumerate_rays(a::Vector{Vector{Int64}})
    m = length(a)
    rk, basis = rank_mod(a)
    rk == D || fail("constraints do not have rank 21")
    rest = setdiff(1:m, basis)
    rays = Vector{Vector{BigInt}}(); tight = UInt64[]
    for i in 1:D
        B = vcat([a[basis[j]]' for j in 1:D if j != i]...)
        x = kernel_vector(B)
        s = sum(a[basis[i]] .* x)
        s == 0 && fail("initial ray tight on its own constraint")
        s < 0 && (x = -x)
        t = UInt64(0)
        for j in 1:D; j != i && (t |= UInt64(1) << (basis[j] - 1)); end
        push!(rays, x); push!(tight, t)
    end
    for k in rest
        ak = a[k]; bit = UInt64(1) << (k - 1)
        sv = [sum(ak .* r) for r in rays]
        pos = findall(>(0), sv); neg = findall(<(0), sv)
        nr = Vector{Vector{BigInt}}(); nt = UInt64[]
        for i in eachindex(rays)
            sv[i] < 0 && continue
            push!(nr, rays[i]); push!(nt, sv[i] == 0 ? tight[i] | bit : tight[i])
        end
        for ip in pos, iq in neg
            S = tight[ip] & tight[iq]
            count_ones(S) < D - 2 && continue
            adj = true
            for r in eachindex(rays)
                (r == ip || r == iq) && continue
                if S & ~tight[r] == 0; adj = false; break; end
            end
            adj || continue
            push!(nr, prim(sv[ip] .* rays[iq] .- sv[iq] .* rays[ip])); push!(nt, S | bit)
        end
        rays, tight = nr, nt
    end
    rays
end

function main(path)
    lines = readlines(path)
    pos = Ref(1)
    nextline() = (l = lines[pos[]]; pos[] += 1; l)
    w = split(nextline()); w[1] == "VERTICES" || fail("format")
    nv = parse(Int, w[2]); total = 0
    for v in 0:nv-1
        split(nextline()) == ["VERTEX", string(v)] || fail("format vertex")
        nextline()
        w = split(nextline()); w[1] == "MIN" || fail("format MIN")
        m = parse(Int, w[2])
        a = [symc(parse.(Int, split(nextline()))) for _ in 1:m]
        w = split(nextline()); w[1] == "GENS" || fail("format GENS")
        for _ in 1:parse(Int, w[2]); nextline(); end
        w = split(nextline()); w[1] == "RAYS" || fail("format RAYS")
        cert = Set([parse.(BigInt, split(nextline())) for _ in 1:parse(Int, w[2])])
        w = split(nextline()); w[1] == "REPS" || fail("format REPS")
        for _ in 1:parse(Int, w[2]); nextline(); end
        nextline() == "END" || fail("format END")
        rays = enumerate_rays(a)
        found = Set(rays)
        ok = length(found) == length(rays) && found == cert
        println("vertex $v: $m constraints, $(length(rays)) extreme rays, certificate $(length(cert)), equal: $(ok ? "yes" : "no")")
        flush(stdout)
        ok || fail("the extreme rays differ from the certificate")
        total += length(rays)
    end
    println("dd.jl: all $nv ray sets agree with the certificate ($total rays)")
end

main(ARGS[1])
