@testset "IO-equation of predator-prey model" begin
    # predator-prey model
    var_names = [
        "x_0", "x_1",
        "a", "b", "c", "d"
    ]
    R, (x_0, x_1, a, b, c, d) = PolynomialRing(QQ, var_names)
    f = [a * x_0 - b * x_0 * x_1, - c * x_1 + d * x_0 * x_1]
    g = x_0

    ode = ODE([x_0, x_1], [a, b, c, d], f, g)
    find_ioequation(ode, true)

    (x_0, x_1, a, b, c, d, x_0_dot, x_1_dot, y_0, y_1, y_2) = gens(ode.poly_ring)
    ind_multi = ode.io_equation // (y_2*y_0 - y_1^2 - y_1*y_0^2*d + y_1*y_0*c + y_0^3*a*d - y_0^2*a*c)
    @test numerator(ind_multi) == to_base_ring(numerator(ind_multi))
end