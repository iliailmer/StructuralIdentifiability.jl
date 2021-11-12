using ModelingToolkit, StructuralIdentifiability

@parameters b0 g zeta0 k gam dlt N
@variables t S(t) A(t) I(t) R(t) Dd(t)  y1(t) [output=true] y2(t) [output=true]

D = Differential(t)

eqs = [
    D(S) ~  -b0 * g * S * I/N - zeta0 * g * S * A / N,
    D(A) ~ b0 * g * S * I/N + zeta0 * g * S * A / N - k * A,
    D(I) ~ k * A - (gam + dlt) * I,
    D(R) ~ gam * I,
    D(Dd) ~ dlt * I,
    y1 ~ I,
    y2 ~ R
]

ode = ODESystem(eqs, t, name=:SAIRD)

res = assess_identifiability(ode)
println(res)