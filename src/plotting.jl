export plot, plot_rectangle, plot_circle, plot_ellipse, plot_bounds, plot_graph, plot_tree, plot_path, plot_line_segments, save_plot

rectangle_corners(lo, hi) = ([lo[1],hi[1],hi[1],lo[1],lo[1]], [lo[2],lo[2],hi[2],hi[2],lo[2]])

function plot_rectangle(r; kwargs...)
    plt.fill(rectangle_corners(r[:,1], r[:,2])...,
             edgecolor="black", zorder=0; kwargs...)
end

function plot_polygon(pts; xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, kwargs...)
    XY = hcat(pts...)'
    plt.fill(clamp(XY[:,1], xmin, xmax),
             clamp(XY[:,2], ymin, ymax),
             edgecolor="black", zorder=0; kwargs...)
end

function plot_circle(c, r; xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, kwargs...)
    plt.fill(clamp(c[1] + r*cos(linspace(0, 2pi, 40)), xmin, xmax),
             clamp(c[2] + r*sin(linspace(0, 2pi, 40)), ymin, ymax),
             edgecolor="black", zorder=0; kwargs...)
end

function plot_ellipse(c, a, b, t; xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, kwargs...)
    XY = [a*cos(linspace(0, 2pi, 40)) b*sin(linspace(0, 2pi, 40))]*[cos(t) sin(t); -sin(t) cos(t)]
    plt.fill(clamp(c[1] + XY[:,1], xmin, xmax),
             clamp(c[2] + XY[:,2], ymin, ymax),
             edgecolor="black", zorder=0; kwargs...)
end

function plot_ellipse(c, Sigma; xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, kwargs...)
    XY = [cos(linspace(0, 2pi, 40)) sin(linspace(0, 2pi, 40))]*chol(Sigma)
    plt.fill(clamp(c[1] + XY[:,1], xmin, xmax),
             clamp(c[2] + XY[:,2], ymin, ymax),
             edgecolor="black", zorder=0; kwargs...)
end

function plot_bounds(lo = zeros(2), hi = ones(2))
    plt.plot(rectangle_corners(lo, hi)..., color="black", linewidth=1.0, linestyle="-")
    axis("equal")
end

function plot_graph(V::Matrix, F; kwargs...)  # learn how to just pass kwargs
    scatter(V[1,:], V[2,:], zorder=1; kwargs...)
    X = vcat([V[1,idx_list] for idx_list in findn(triu(F))]..., fill(nothing, 1, sum(triu(F))))[:]
    Y = vcat([V[2,idx_list] for idx_list in findn(triu(F))]..., fill(nothing, 1, sum(triu(F))))[:]
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end
plot_graph(V::Vector, F; kwargs...) = plot_graph(hcat(V...), F; kwargs...)

function plot_tree(V::Matrix, A; kwargs...)
    scatter(V[1,:], V[2,:], zorder=1; kwargs...)
    X = vcat(V[1,find(A)], V[1,A[find(A)]], fill(nothing, 1, countnz(A)))[:]
    Y = vcat(V[2,find(A)], V[2,A[find(A)]], fill(nothing, 1, countnz(A)))[:]
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end
plot_tree(V::Vector, A; kwargs...) = plot_tree(hcat(V...), A; kwargs...)

function plot_line_segments(P1::Vector, P2::Vector; kwargs...)
    X = [[v[1] for v in P1]'; [v[1] for v in P2]'; fill(nothing, 1, length(P1))][:]
    Y = [[v[2] for v in P1]'; [v[2] for v in P2]'; fill(nothing, 1, length(P1))][:]
    plt.plot(X, Y; kwargs...)
end

function plot_path(V::Matrix, idx_list = 1:size(V,2); kwargs...)
    plt.plot(V[1,idx_list]', V[2,idx_list]', linewidth=1.0, linestyle="-", zorder=2; kwargs...)
end
plot_path(V::Vector, idx_list = 1:length(V); kwargs...) = plot_path(hcat(V...), idx_list; kwargs...)

function save_plot(fname, title = "")
    axis("off")
    axes()[:get_xaxis]()[:set_visible](false)
    axes()[:get_yaxis]()[:set_visible](false)
    title != "" && title(title, fontsize=20)
    savefig(fname, bbox_inches="tight")
end

# using Compose

# function composable(C::Circle, opts = (fill("red"), fillopacity(1.)))
#     compose(context(), circle(C.c..., C.r), opts...)
# end

# function composable(P::Polygon, opts = (fill("red"), fillopacity(1.)))
#     compose(context(), polygon([(p[1],p[2]) for p in P.points]), opts...)
# end

# function composable(C::Compound2D, opts = (fill("red"), fillopacity(1.)))
#     compose(context(), [composable(P, opts) for P in C.parts]...)
# end

# function composable(L::Line, opts = (stroke("black"),))
#     compose(context(), line([(L.v[1], L.v[2]), (L.w[1], L.w[2])]), opts...)
# end

# function composable(SS::StateSpace, opts = (linewidth(2px), fill(nothing), stroke("black")); flipy=true)
#     w, h = SS.hi - SS.lo
#     compose(context(-SS.lo[1]/w, -SS.lo[2]/h, 1/w, 1/h),
#             compose(context(), rectangle(SS.lo[1], flipy ? SS.hi[2] : SS.lo[2], w, h), opts...))
# end

# flip(C::Context) = compose(context(mirror=Mirror()), C)

# function plot(P::MPProblem)
#     p = compose(composable(P.SS), composable(P.CC))
#     if P.status == :solved
#         println("IMPLEMENT ME!")
#     end
#     flip(p)
# end