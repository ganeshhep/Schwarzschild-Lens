using Plots

function grlens_plots(alpha_0, beta_x_0, beta_y_0, sigma_x, sigma_y; N, theta_min, theta_max, tol)

    theta_vals = range(theta_min, theta_max, length = N)
    beta_vals = range(theta_min, theta_max, length = N)

    I_source = zeros(N, N)
    I_image = zeros(N, N)
 
    amp = 1/(pi * sigma_x * sigma_y)

    for i_src in 1:N, j_src in 1:N
        beta_x_src, beta_y_src = beta_vals[i_src], beta_vals[j_src]

        I_source[i_src, j_src] = amp * exp(- 0.5 * ((beta_x_src - beta_x_0)^2/sigma_x^2 + 
                                                    (beta_y_src - beta_y_0)^2/sigma_y^2))

        intensity = I_source[i_src, j_src]
 
        for i in 1:N, j in 1:N
            theta_x, theta_y = theta_vals[i], theta_vals[j]
       
            r2 = theta_x^2 + theta_y^2

            if r2 == 0
               continue
            end

            beta_x = theta_x - (alpha_0^2/r2) * theta_x
            beta_y = theta_y - (alpha_0^2/r2) * theta_y
        
            if sqrt((beta_x - beta_x_src)^2 + (beta_y - beta_y_src)^2) < tol
               I_image[i, j] = intensity
            end 
        end
    end

    return theta_vals, beta_vals, I_source, I_image
end

theta_vals, beta_vals, I_source, I_image = grlens_plots(0.5, 0.0, 0.0, 0.5, 0.5; N = 1500,
                                                        theta_min = - 5.0, theta_max = 5.0, tol = 1e-6)


p1 = heatmap(beta_vals, beta_vals, I_source', aspect_ratio = 1, xlabel = "βx", ylabel = "βy",
             title = "Source Plane", color=:inferno)

p2 = heatmap(theta_vals, theta_vals, I_image', aspect_ratio = 1, xlabel = "θx", ylabel = "θy",
             title = "Image Plane", color=:inferno)

plot(p1, p2, layout = (1,2), size = (1000, 400))

savefig("lensing_gaussian.png")
        
