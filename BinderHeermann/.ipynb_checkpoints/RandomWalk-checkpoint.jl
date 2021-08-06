module RandomWalk

using PyPlot
PyPlot.rc("figure", figsize = (4, 4));

function single_2DRW(ran_vec::Vector{Float64})
    
    """ Generate a single 2D random walk
    ran_vec = vector of random numbers
    
    Returns a tuple of x and y position vectors (x⃗,y⃗)
    """
    x = 0
    y = 0
    x_vec = zeros(Int64,length(ran_vec)+1)
    y_vec = zeros(Int64,length(ran_vec)+1)
    for (i,ran_f) in enumerate(ran_vec)
        ran_i::Int64 = ceil(4*ran_f)
        if ran_i == 1
            x += 1
        elseif ran_i == 2
            y += 1
        elseif ran_i == 3
            x -= 1
        else
            y -= 1
        end
        x_vec[i+1] = x
        y_vec[i+1] = y
    end
    return x_vec, y_vec
end;

function plot_single_walk(x_vec::Vector{Int64}, y_vec::Vector{Int64})
    """ Plot a single 2D random walk
    x_vec = x-coordinate vector
    y_vec = y-coordinate vector
    """
    # get arrows
    xx_vec = [xi + (xf - xi)/2 for (xi, xf) in zip(x_vec[1:1:end], x_vec[2:1:end])]
    yy_vec = [yi + (yf - yi)/2 for (yi, yf) in zip(y_vec[1:1:end], y_vec[2:1:end])]
    vx_vec = [xf - xi for (xi, xf) in zip(x_vec[1:1:end], x_vec[2:1:end])]
    vy_vec = [yf - yi for (yi, yf) in zip(y_vec[1:1:end], y_vec[2:1:end])]

    # plot
    plot(x_vec, y_vec, color="b", lw=3)
    quiver(xx_vec, yy_vec, vx_vec, vy_vec,
        color="b", scale_units="xy", scale=2, pivot="mid", 
        headwidth=10, headlength=10, headaxislength=5)
    scatter(first(x_vec), first(y_vec), s=60, marker="o", color="b")
    scatter(last(x_vec), last(y_vec), s=60, marker="s", color="b")

    # formatting
    xlo, xhi = minimum(x_vec), maximum(x_vec)
    ylo, yhi = minimum(y_vec), maximum(y_vec)
    xlabel("X")
    ylabel("Y")
    vlines(xlo:xhi, ylo, yhi, color="k", alpha=0.2)
    hlines(ylo:yhi, xlo, xhi, color="k", alpha=0.2)
    xticks([]), yticks([])
    axis("equal")
    axis("off");
end

function sample_2DRW(ran_vec::Vector{Float64}, nsamples::Int64, nsteps::Int64)
    """ Sample random walks in two dimensions
    nsamples = number of random walk samples
    nsteps = number of steps of each random walk
    
    Returns a named tuple of averages (x2=<x²>, y2=<y²>, x4=<x⁴>, y4=<y⁴>)
    """
    if length(ran_vec) != nsamples*nsteps
        throw("Wrong ran_vec dimensions! It should be nsamples*nsteps")
    end
    x2av = 0.0
    y2av = 0.0
    x4av = 0.0
    y4av = 0.0
    for i ∈ 1:nsamples
        xi = 0.0
        yi = 0.0
        for j ∈ 1:nsteps
            rn::Int64 = ceil(4*ran_vec[nsteps*(i-1)+j])
            if rn == 1
                xi += 1.0
            elseif rn == 2
                yi += 1.0
            elseif rn == 3
                xi -= 1.0
            else
                yi -= 1.0
            end
        end
        xsq = xi^2
        ysq = yi^2
        x2av += xsq
        y2av += ysq
        x4av += xsq*xsq
        y4av += ysq*ysq
    end
    x2av /= nsamples
    y2av /= nsamples
    x4av /= nsamples
    y4av /= nsamples
    av = (x2=x2av, y2=y2av, x4=x4av, y4=y4av)
    return av
end;

function plot_RNG_correlations(ran_vec::Vector{Float64})
    """ Plot correlations between random numbers
    ran_vec = vector of random numbers (0,1)
    """
    if isodd(length(ran_vec))
        throw("Wrong ran_vec dimensions! its length should be an even number")
    end
    len = Int(length(ran_vec)/2)
    x = 2*len*ran_vec[1:len]
    y = 2*len*ran_vec[len+1:2*len]
    scatter(x,y, s=4, alpha=0.5)
    xticks([])
    yticks([])
end;

end