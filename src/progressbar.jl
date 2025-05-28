using Printf

# see https://github.com/JuliaLang/julia/blob/8f5b7ca12ad48c6d740e058312fc8cf2bbe67848/base/util.jl#L117-L140 for colours

const default_str = "\033[39m\033[0m"
const bold_str = "\033[1m"
const blue_str = "\033[34m"
const red_str = "\033[31m"
const green_str = "\033[32m"
const begin_str = "\r\033[K" # clear the line from the beginning, then reset the cursor to the beginning
const format_str = Printf.Format(*(begin_str, bold_str, green_str, "%i iterations", default_str, ", action/(8π²/%i) = ", bold_str, blue_str, "%.5f", default_str, " + ", bold_str, red_str, "%.5f", default_str, " = ", bold_str, "%.12f", default_str, ", log(convergence) = (%c)%.4e"))

function charactersign(n::Real)::Char
    if n >= 0
        return '+'
    else
        return '-'
    end
end


function printprogressbar(printoutput::IO, L::Lattice, action::Vector{Float64})
    print(printoutput, Printf.format(format_str, length(action)-1, L.N, L.electricAction/(8*pi^2/L.N), L.magneticAction/(8*pi^2/L.N), L.action/(8*pi^2/L.N), charactersign(action[end]-action[end-1]), log10(abs(action[end]-action[end-1]))))
    flush(printoutput)
end