function point_source_image(alpha_0, beta_x_src, beta_y_src, theta_min, theta_max, N, tol)

    theta_x_vals = range(theta_min, theta_max, length = N)
    theta_y_vals = range(theta_min, theta_max, length = N)

    theta_x_sol = Float64[]
    theta_y_sol = Float64[]

    for theta_x in theta_x_vals, theta_y in theta_y_vals
        r2 = theta_x^2 + theta_y^2

        if r2 == 0.0
            continue
        end

        beta_x = theta_x - (alpha_0^2 / r2) * theta_x
        beta_y = theta_y - (alpha_0^2 / r2) * theta_y

        dist = sqrt((beta_x - beta_x_src)^2 + (beta_y - beta_y_src)^2)

        if dist < tol
           push!(theta_x_sol, theta_x)
           push!(theta_y_sol, theta_y)
        end
    end

    return theta_x_sol, theta_y_sol
end

theta_x_pts, theta_y_pts = point_source_image(1.0, 0.0, 0.0, -2.0, 2.0, 5000, 1e-2)

using Plots
p = scatter(theta_x_pts, theta_y_pts, markersize = 1.0, xlabel = "theta_x [arcs]", ylabel = "theta_y [arcs]",
            title = "Lens Plane", size = (1000, 1000), dpi = 300)
savefig(p, "point_source_image.png")

