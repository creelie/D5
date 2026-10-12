# Independent checker (Julia, exact rationals) for the certificate of the two-translate problem in R^5.
# Usage: julia check.jl certificate.txt

const QQ = Rational{BigInt}

fail(m) = (println(stderr, "FAIL: ", m); exit(1))

qform(J, k) = sum(J[i, j] * k[i] * k[j] for i in 1:6, j in 1:6)
symc(k) = [k[i] * k[j] * (i == j ? 1 : 2) for i in 1:6 for j in i:6]

function det_exact(A0)
    A = copy(A0); n = size(A, 1); d = one(QQ)
    for c in 1:n
        p = findfirst(r -> A[r, c] != 0, c:n)
        p === nothing && return zero(QQ)
        p += c - 1
        if p != c
            A[[p, c], :] = A[[c, p], :]; d = -d
        end
        d *= A[c, c]
        for r in c+1:n
            if A[r, c] != 0
                f = A[r, c] / A[c, c]
                A[r, c:n] .-= f .* A[c, c:n]
            end
        end
    end
    d
end

# Q[x] = sum_i d_i (x_i + sum_{j>i} m_ij x_j)^2
function ldl(A)
    n = size(A, 1); d = zeros(QQ, n); m = zeros(QQ, n, n)
    for i in 1:n
        v = A[i, i] - sum(QQ[d[k] * m[k, i]^2 for k in 1:i-1]; init = zero(QQ))
        v <= 0 && return nothing
        d[i] = v
        for j in i+1:n
            w = A[i, j] - sum(QQ[d[k] * m[k, i] * m[k, j] for k in 1:i-1]; init = zero(QQ))
            m[i, j] = w / v
        end
    end
    (d, m)
end

# all integer x with Q[x + c] <= B
function enumerate_short(d, m, c, B, cb)
    n = length(d); x = zeros(Int, n)
    function rec(i, R)
        if i == 0
            cb(x); return
        end
        u = c[i] + sum(QQ[m[i, j] * (x[j] + c[j]) for j in i+1:n]; init = zero(QQ))
        x0 = Int(floor(-u))
        v = x0 + 1
        while true
            val = d[i] * (v + u)^2
            val > R && break
            x[i] = v; rec(i - 1, R - val); v += 1
        end
        v = x0
        while true
            val = d[i] * (v + u)^2
            val > R && break
            x[i] = v; rec(i - 1, R - val); v -= 1
        end
        x[i] = 0
    end
    B >= 0 && rec(n, B)
end

function in_R(J)
    A = J[1:5, 1:5]
    f = ldl(A); f === nothing && fail("Q block not positive definite")
    d, m = f
    ok = true; mins = Vector{Vector{Int}}()
    enumerate_short(d, m, zeros(QQ, 5), QQ(4), x -> begin
        all(==(0), x) && return
        x[findfirst(!=(0), x)] < 0 && return
        k = [x; 0]; v = qform(J, k)
        if v < 4; ok = false; elseif v == 4; push!(mins, k); end
    end)
    r = J[1:5, 6]
    q = A \ r
    sigma = J[6, 6] - sum(r .* q)
    enumerate_short(d, m, q, QQ(4) - sigma, x -> begin
        k = [x; 1]; v = qform(J, k)
        if v < 4; ok = false; elseif v == 4; push!(mins, k); end
    end)
    ok, mins
end

const P61 = Int128(2)^61 - 1
function rank_mod(rows)
    isempty(rows) && return 0
    A = [mod.(Int128.(r), P61) for r in rows]
    n = length(A[1]); rk = 0
    for c in 1:n
        rk == length(A) && break
        p = findfirst(r -> A[r][c] != 0, rk+1:length(A))
        p === nothing && continue
        p += rk
        A[p], A[rk+1] = A[rk+1], A[p]
        rk += 1
        inv = powermod(A[rk][c], P61 - 2, P61)
        for r in 1:length(A)
            if r != rk && A[r][c] != 0
                f = mod(A[r][c] * inv, P61)
                A[r] = mod.(A[r] .- f .* A[rk], P61)
            end
        end
    end
    rk
end

# polynomials: coefficient vectors, low degree first
trim(p) = (q = copy(p); while !isempty(q) && q[end] == 0; pop!(q); end; q)
evalp(p, t) = foldr((a, s) -> s * t + a, p; init = zero(QQ))
deriv(p) = trim([p[i] * (i - 1) for i in 2:length(p)])
function mulp(a, b)
    (isempty(a) || isempty(b)) && return QQ[]
    c = zeros(QQ, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b); c[i+j-1] += a[i] * b[j]; end
    trim(c)
end
function divmodp(a, b)
    r = trim(a); q = zeros(QQ, max(length(r) - length(b) + 1, 0))
    while !isempty(r) && length(r) >= length(b)
        f = r[end] / b[end]; sh = length(r) - length(b)
        q[sh+1] = f
        for i in eachindex(b); r[sh+i] -= f * b[i]; end
        r = trim(r)
    end
    trim(q), r
end
function gcdp(a, b)
    a = trim(a); b = trim(b)
    while !isempty(b); _, r = divmodp(a, b); a, b = b, r; end
    isempty(a) ? a : a ./ a[end]
end
function squarefree(p)
    q, r = divmodp(p, gcdp(p, deriv(p))); isempty(r) || fail("sqf"); q
end
function sturm(A, B)
    s = [A]; b = trim(B)
    isempty(b) && return s
    push!(s, b)
    while true
        _, r = divmodp(s[end-1], s[end])
        isempty(r) && break
        push!(s, -r)
    end
    s
end
function changes(s, t)
    prev = 0; ch = 0
    for p in s
        v = t === nothing ? sign(p[end]) : sign(evalp(p, t))
        v == 0 && continue
        (prev != 0 && v != prev) && (ch += 1)
        prev = v
    end
    ch
end
taq(P, G, a, b) = length(P) <= 1 ? 0 : (s = sturm(P, mulp(deriv(P), G)); changes(s, a) - changes(s, b))
function interp(v)
    n = length(v); c = copy(v)
    for j in 1:n-1, i in n:-1:j+1
        c[i] = (c[i] - c[i-1]) / j
    end
    p = QQ[c[n]]
    for i in n-1:-1:1
        p = mulp(p, QQ[-(i - 1), 1]); isempty(p) && (p = QQ[0]); p[1] += c[i]; p = trim(p)
    end
    trim(p)
end
function drop_root(P, t)
    evalp(P, t) == 0 || return P
    q, r = divmodp(P, QQ[-t, 1]); q
end

act(T, Y) = transpose(T) * Y * T
function symm(y)
    Y = zeros(QQ, 6, 6); t = 1
    for i in 1:6, j in i:6; Y[i, j] = y[t]; Y[j, i] = y[t]; t += 1; end
    Y
end
# the integer 21x21 matrix of Y -> T^T Y T in the coordinates y_ij, i <= j
function symmap(T)
    Ti = Int.(T); L = zeros(Int, 21, 21); t = 1
    for i in 1:6, j in i:6
        S = zeros(Int, 6, 6); S[i, j] = 1; S[j, i] = 1
        Y = transpose(Ti) * S * Ti; u = 1
        for a in 1:6, b in a:6; L[u, t] = Y[a, b]; u += 1; end
        t += 1
    end
    L
end

function check_gamma(T)
    all(T[6, 1:5] .== 0) || fail("T not in Gamma (last row)")
    abs(T[6, 6]) == 1 || fail("T not in Gamma (eps)")
    abs(det_exact(QQ.(T[1:5, 1:5]))) == 1 || fail("T not in Gamma (det U)")
end

function parse_q(s)
    if occursin('/', s)
        a, b = split(s, '/'); return QQ(parse(BigInt, a), parse(BigInt, b))
    end
    QQ(parse(BigInt, s))
end

function main(path)
    toks = split(read(path, String))
    pos = Ref(1)
    nexttok() = (t = toks[pos[]]; pos[] += 1; String(t))
    nextint() = parse(Int, nexttok())
    nexttok() == "VERTICES" || fail("format")
    nv = nextint()
    V = []
    for v in 0:nv-1
        (nexttok() == "VERTEX" && nextint() == v) || fail("format vertex")
        nexttok()
        J = zeros(QQ, 6, 6)
        for i in 1:6, j in 1:6; J[i, j] = parse_q(nexttok()); end
        nexttok(); m = nextint(); mins = [[nextint() for _ in 1:6] for _ in 1:m]
        nexttok(); g = nextint(); gens = [permutedims(reshape([nextint() for _ in 1:36], 6, 6)) for _ in 1:g]
        nexttok(); nr = nextint(); rays = [[nextint() for _ in 1:21] for _ in 1:nr]
        nexttok(); nreps = nextint(); reps = []
        for _ in 1:nreps
            kind = nexttok(); idx = nextint(); sz = nextint()
            if kind == "B"
                ts = parse_q(nexttok()); to = nextint()
                T = permutedims(reshape([nextint() for _ in 1:36], 6, 6))
                push!(reps, (kind, idx, sz, ts, to, T))
            else
                push!(reps, (kind, idx, sz, nothing, -1, nothing))
            end
        end
        nexttok() == "END" || fail("format end")
        push!(V, (J, mins, gens, rays, reps))
    end
    nexttok(); ne = nextint()
    eqpts = [(nexttok(); (nexttok(), nexttok(), nexttok())) for _ in 1:ne]

    total_rays = 0; total_edges = 0; eq_found = 0; eq_mine = []
    for v in 0:nv-1
        J, lmins, gens, rays, reps = V[v+1]
        J == transpose(J) || fail("J not symmetric")
        ok, mins = in_R(J); ok || fail("vertex not in R")
        sort(lmins) == sort(mins) || fail("Min(J) differs at vertex $v")
        rank_mod([symc(k) for k in mins]) == 21 || fail("not a vertex")
        for T in gens
            check_gamma(T); act(QQ.(T), J) == J || fail("generator is not an automorphism")
        end
        posd = Dict{Vector{Int}, Int}()
        for (i, y) in enumerate(rays)
            haskey(posd, y) && fail("repeated ray")
            posd[y] = i
            gcd(y) == 1 || fail("ray not primitive")
            tight = Vector{Vector{Int}}()
            for k in mins
                a = symc(k); s = sum(a .* y)
                s < 0 && fail("ray outside the cone")
                s == 0 && push!(tight, a)
            end
            rank_mod(tight) == 20 || fail("ray not extreme")
        end
        par = collect(1:length(rays))
        findr(a) = (while par[a] != a; par[a] = par[par[a]]; a = par[a]; end; a)
        for T in gens
            L = symmap(T)
            for (i, y) in enumerate(rays)
                yy = L * y
                haskey(posd, yy) || fail("image of a ray is not a ray")
                a = findr(i); b = findr(posd[yy]); a != b && (par[a] = b)
            end
        end
        osize = Dict{Int, Int}()
        for i in eachindex(rays); r = findr(i); osize[r] = get(osize, r, 0) + 1; end
        seen = Set{Int}()
        for (kind, idx, sz, ts, to, T) in reps
            o = findr(idx + 1)
            o in seen && fail("two representatives in one orbit")
            push!(seen, o)
            osize[o] == sz || fail("orbit size")
        end
        length(seen) == length(osize) || fail("an orbit has no representative")
        total_rays += length(rays)
        dJ = det_exact(J); dQ = det_exact(J[1:5, 1:5])
        (dJ <= 0 && dQ < 512) && fail("vertex with det J <= 0 and det Q < 512")
        (dJ == 0 && dQ == 512) && push!(eq_mine, (string(v), "-1", "0"))
        for (kind, idx, sz, ts, to, T) in reps
            Y = symm(rays[idx+1]); Tp = nothing
            if kind == "U"
                i0 = findfirst(i -> Y[i, i] != 0, 1:5); s = Y[6, 6]
                if i0 === nothing
                    any(Y[1:5, :] .!= 0) && fail("unbounded: Q_Y = 0 but Y != s e e^T")
                    s < 0 && fail("unbounded: s < 0")
                else
                    Y[i0, i0] < 0 && fail("unbounded: negative diagonal")
                    w = Y[1:5, i0] ./ Y[i0, i0]
                    L = lcm(denominator.(w)); w = w .* L; w = w ./ gcd(numerator.(w))
                    lam = Y[i0, i0] / w[i0]^2
                    Y[1:5, 1:5] == lam .* (w * w') || fail("unbounded: Q_Y not of rank one")
                    mu = Y[i0, 6] / w[i0]
                    Y[1:5, 6] == mu .* w || fail("unbounded: r_Y not parallel to w")
                    n0 = floor(-mu / lam)
                    for n in (n0, n0 + 1)
                        lam * n^2 + 2 * mu * n + s < 0 && fail("unbounded: negative value")
                    end
                end
            else
                ts <= 0 && fail("t* <= 0")
                M = J .+ ts .* Y
                ok2, mm = in_R(M); ok2 || fail("end point not in R")
                any(qform(Y, k) < 0 for k in mm) || fail("t* not maximal")
                (0 <= to < nv) || fail("target out of range")
                check_gamma(T)
                act(QQ.(T), V[to+1][1]) == M || fail("end point not equivalent to the target vertex")
                Tp = ts
            end
            vJ = [det_exact(J .+ t .* Y) for t in 0:6]
            vQ = [det_exact((J .+ t .* Y)[1:5, 1:5]) for t in 0:5]
            pJ = interp(QQ.(vJ)); pQ = interp(QQ.(vQ))
            g = isempty(pQ) ? QQ[-512] : (gg = copy(pQ); gg[1] -= 512; trim(gg))
            if isempty(pJ)
                isempty(g) && fail("det Q constant 512 on an edge")
                mid = Tp === nothing ? QQ(1) : Tp / 2
                evalp(g, mid) > 0 || fail("zero edge: det Q < 512 inside")
                h = drop_root(squarefree(g), QQ(0))
                Tp !== nothing && (h = drop_root(h, Tp))
                taq(h, QQ[1], QQ(0), Tp) == 0 || fail("zero edge: det Q - 512 has a root inside")
                total_edges += 1
                continue
            end
            P = drop_root(squarefree(pJ), QQ(0))
            Tp !== nothing && (P = drop_root(P, Tp))
            N = taq(P, QQ[1], QQ(0), Tp)
            S = taq(P, g, QQ(0), Tp)
            h = isempty(g) ? P : gcdp(P, g)
            N0 = taq(h, QQ[1], QQ(0), Tp)
            iseven(N - N0 - S) || fail("parity")
            (N - N0 - S) ÷ 2 == 0 || fail("a point with det J = 0 and det Q < 512")
            if N0 > 0
                for e in eqpts
                    if e[1] == string(v) && e[2] == string(idx)
                        t = parse_q(e[3])
                        (evalp(pJ, t) == 0 && evalp(g, t) == 0) || fail("listed equality point is not one")
                        push!(eq_mine, e)
                    end
                end
                eq_found += N0
            end
            total_edges += 1
        end
        println("vertex $v: ok"); flush(stdout)
    end
    nvx = count(e -> e[2] == "-1", eq_mine)
    (length(eq_mine) == eq_found + nvx && length(eqpts) == length(eq_mine)) || fail("equality points do not match the certificate")
    println("vertices: $nv, rays: $total_rays, edge representatives: $total_edges, equality points: $(length(eq_mine))")
    println("Julia checker: all checks passed")
end

main(ARGS[1])
