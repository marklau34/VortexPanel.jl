# Minimal SVG line plots for the validation scripts, without plotting dependencies. Requires Printf.

const SVG_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8e44ad", "#52514e", "#c0392b", "#16a085"]

svg_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

"""Tick positions covering [lo, hi] at a 1-2-5 step giving about n ticks"""
function nice_ticks(lo::Real, hi::Real; n::Integer=6)
    span = hi - lo
    span > 0 || return [lo]
    raw = span / n
    mag = 10.0^floor(log10(raw))
    step = mag * (raw / mag <= 1 ? 1 : raw / mag <= 2 ? 2 : raw / mag <= 5 ? 5 : 10)
    return collect(ceil(lo / step)*step:step:hi+1e-9step)
end

tick_label(v) = abs(v) < 1e-12 ? "0" : (a = abs(v); a >= 1e4 || a < 1e-3 ? @sprintf("%.1e", v) : @sprintf("%.4g", v))

"""
    svgplot(file; series, title="", xlabel="", ylabel="", xlims=nothing, ylims=nothing, hline=nothing, width=660, height=430)

Write a line plot as SVG. `series` holds NamedTuples with `x`, `y` and `label`, and optionally `color`, `dash` (true for
a dashed line), `marker` (false for none) and `line` (false for markers only). `hline` draws a dashed reference line.
"""
function svgplot(file::AbstractString; series, title="", xlabel="", ylabel="", xlims=nothing, ylims=nothing, hline=nothing,
    width::Integer=660, height::Integer=430)
    L, R, T, B = 78, 22, 42, 58
    pw, ph = width - L - R, height - T - B
    xs = [x for s in series for (x, y) in zip(s.x, s.y) if isfinite(x) && isfinite(y)]
    ys = [y for s in series for (x, y) in zip(s.x, s.y) if isfinite(x) && isfinite(y)]
    isnothing(hline) || push!(ys, hline)
    xlo, xhi = isnothing(xlims) ? extrema(xs) : xlims
    ylo, yhi = isnothing(ylims) ? extrema(ys) : ylims
    xhi == xlo && (xlo -= 0.5; xhi += 0.5)
    yhi == ylo && (ylo -= 0.5; yhi += 0.5)
    if isnothing(ylims) # a little headroom
        pad = 0.05 * (yhi - ylo)
        ylo, yhi = ylo - pad, yhi + pad
    end
    px(x) = L + (x - xlo) / (xhi - xlo) * pw
    py(y) = T + (yhi - y) / (yhi - ylo) * ph
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height" font-family="Helvetica, Arial, sans-serif" font-size="12">""")
    println(io, """<rect width="$width" height="$height" fill="#fcfcfb"/>""")
    for t in nice_ticks(xlo, xhi)
        x = px(t)
        @printf(io, "<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e4e3df\"/>\n", x, T, x, T + ph)
        @printf(io, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" fill=\"#3a3936\">%s</text>\n", x, T + ph + 16, svg_escape(tick_label(t)))
    end
    for t in nice_ticks(ylo, yhi)
        y = py(t)
        @printf(io, "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e4e3df\"/>\n", L, y, L + pw, y)
        @printf(io, "<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" fill=\"#3a3936\">%s</text>\n", L - 6, y + 4, svg_escape(tick_label(t)))
    end
    @printf(io, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#52514e\"/>\n", L, T, pw, ph)
    if !isnothing(hline) && ylo <= hline <= yhi
        @printf(io, "<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#52514e\" stroke-dasharray=\"5,4\"/>\n", L, py(hline), L + pw, py(hline))
    end
    println(io, """<defs><clipPath id="area"><rect x="$L" y="$T" width="$pw" height="$ph"/></clipPath></defs><g clip-path="url(#area)">""")
    for (i, s) in enumerate(series)
        color = get(s, :color, SVG_COLORS[mod1(i, length(SVG_COLORS))])
        pts = [(px(x), py(y)) for (x, y) in zip(s.x, s.y) if isfinite(x) && isfinite(y)]
        isempty(pts) && continue
        if get(s, :line, true) && length(pts) > 1
            dash = get(s, :dash, false) ? " stroke-dasharray=\"7,4\"" : ""
            println(io, "<polyline fill=\"none\" stroke=\"$color\" stroke-width=\"2\"$dash points=\"", join([@sprintf("%.1f,%.1f", p...) for p in pts], " "), "\"/>")
        end
        if get(s, :marker, true)
            for p in pts
                @printf(io, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"3.2\" fill=\"%s\"/>\n", p..., color)
            end
        end
    end
    println(io, "</g>")
    # legend, top right
    labelled = [(i, s) for (i, s) in enumerate(series) if !isempty(string(get(s, :label, "")))]
    if !isempty(labelled)
        lw = 26 + 7 * maximum(length(string(s.label)) for (_, s) in labelled)
        lx, ly = L + pw - lw - 8, T + 8
        @printf(io, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"#ffffff\" fill-opacity=\"0.88\" stroke=\"#d5d4cf\"/>\n", lx, ly, lw, 18 * length(labelled) + 6)
        for (k, (i, s)) in enumerate(labelled)
            color = get(s, :color, SVG_COLORS[mod1(i, length(SVG_COLORS))])
            y = ly + 15 + 18 * (k - 1)
            dash = get(s, :dash, false) ? " stroke-dasharray=\"7,4\"" : ""
            get(s, :line, true) && @printf(io, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"2\"%s/>\n", lx + 6, y - 4, lx + 22, y - 4, color, dash)
            get(s, :marker, true) && @printf(io, "<circle cx=\"%d\" cy=\"%d\" r=\"3.2\" fill=\"%s\"/>\n", lx + 14, y - 4, color)
            @printf(io, "<text x=\"%d\" y=\"%d\" fill=\"#1f1e1c\">%s</text>\n", lx + 28, y, svg_escape(s.label))
        end
    end
    @printf(io, "<text x=\"%d\" y=\"%d\" font-size=\"14\" font-weight=\"bold\" fill=\"#1f1e1c\">%s</text>\n", L, 24, svg_escape(title))
    @printf(io, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" fill=\"#1f1e1c\">%s</text>\n", L + pw / 2, height - 14, svg_escape(xlabel))
    @printf(io, "<text x=\"16\" y=\"%.1f\" text-anchor=\"middle\" fill=\"#1f1e1c\" transform=\"rotate(-90 16 %.1f)\">%s</text>\n", T + ph / 2, T + ph / 2, svg_escape(ylabel))
    println(io, "</svg>")
    write(file, String(take!(io)))
    return file
end
