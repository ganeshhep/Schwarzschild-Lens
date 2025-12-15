using Plots

function lensing_plots(D_s, D_d, D_ds, R_s; eta_tilde_min, eta_tilde_max, neta_tilde, nphi)

    # Lensing Scales

    alpha_0 = sqrt((2 * R_s * D_ds) / (D_d * D_s))
    eta_0   = alpha_0 * D_s
    xi_0    = alpha_0 * D_d

    # Uniformly Distributed Sources in a Annular Plane Disk

    phi = range(0, 2pi, length = nphi)
    
    eta_tilde_src = range(eta_tilde_min, eta_tilde_max, length = neta_tilde)

    p = scatter([], [], label = false, xlabel = "x [pc]", ylabel = "y [pc]", title = "Source Plane",
                size = (1000, 1000), dpi = 300)
    q = scatter([], [], label = false, xlabel = "x [pc]", ylabel = "y [pc]", title = "Lens Plane",
                size = (1000, 1000), dpi = 300)

    for eta_tilde in eta_tilde_src

        # Source Positions

        x_src = @. eta_tilde * eta_0 * cos(phi)
        y_src = @. eta_tilde * eta_0 * sin(phi)

        # Image Positions

        xi_plus  = 0.5 * (eta_tilde + sqrt(eta_tilde^2 + 4)) * xi_0
        xi_minus = 0.5 * (eta_tilde - sqrt(eta_tilde^2 + 4)) * xi_0

        x_img_1 = @. xi_plus * cos(phi)
        y_img_1 = @. xi_plus * sin(phi)

        x_img_2 = @. xi_minus * cos(phi)
        y_img_2 = @. xi_minus * sin(phi)

        # Source Plot in the Source Plane

        scatter!(p, x_src, y_src, markersize = 1.0, label = (eta_tilde == eta_tilde_min ?
                 "Source" : false))

        # Image Plot in the Lens Plane

        scatter!(q, x_img_1, y_img_1, markersize = 1.0, label = (eta_tilde == eta_tilde_min ?
                 "Source Image" : false))
        
        scatter!(q, x_img_2, y_img_2, markersize = 1.0, label = false)

    end
    
    # Save Figures

    savefig(p, "source_annular_plane_disk.png")
    savefig(q, "source_annular_plane_disk_image.png")

end

# Distances in Parsec

D_s = 2e3
D_d = 1e3
D_ds = 1e3
R_s = 5.4776e-12
 
lensing_plots(D_s, D_d, D_ds, R_s; eta_tilde_min = 1.5, eta_tilde_max = 3.0,
              neta_tilde = 500, nphi = 1500)
