using HDF5

function createdatafile(filename::AbstractString; directory::AbstractString="../data/")::HDF5.File
    f = h5open(*(directory,filename), "w", swmr=true)
    create_group(f, "lattice")
    create_group(f, "metadata")
    create_group(f, "analysis")
    create_group(f, "snapshots")
    return f
end

function savelattice!(L::Lattice, f::Union{HDF5.File,HDF5.Group}; calculateaction::Bool=true)
    if ! haskey(f, "lattice")
        create_group(f, "lattice")
    end
    g = f["lattice"]

    # delete current entries in the file if they exist
    for param in fieldnames(Lattice)
        if startswith(String(param), "_") # skip neighbours, checkerboards, etc
            continue
        end
        if haskey(g, String(param))
            delete_object(g, String(param))
        end
    end
    
    g["sites"] = L.sites
    g["N"] = L.N
    g["size"] = [L.size...]
    
    if calculateaction
        setactiondensity!(L)
        g["actionDensity"] = L.actionDensity
        g["action"] = wilsonaction!(L)
        g["electricAction"] = L.electricAction
        g["magneticAction"] = L.magneticAction
    else
        g["action"] = L.action
        g["electricAction"] = L.electricAction
        g["magneticAction"] = L.magneticAction
        g["actionDensity"] = L.actionDensity
    end

    g["twists"] = L.twists
end

function savesnapshot!(L::Lattice, f::HDF5.File; overwriteprevious::Bool=false)
    if ! haskey(f, "snapshots")
        create_group(f, "snapshots")
    end
    g = f["snapshots"]
    n = 1
    while haskey(g, "snapshot_$n")
        n += 1
    end
    if overwriteprevious && n > 1
        savelattice!(L, g["snapshot_$(n-1)"]; calculateaction=false)
    else
        create_group(g, "snapshot_$n")
        savelattice!(L, g["snapshot_$n"]; calculateaction=false)
    end
end 

function Lattice(f::Union{HDF5.File,HDF5.Group})
    g = f["lattice"]

    L = Lattice()
    for field in fieldnames(Lattice)
        if startswith(String(field), "_") # skip neighbours, checkerboards, etc
            continue
        end
        if haskey(g, String(field))
            if field == :size
                setfield!(L, field, Tuple(read(g[String(field)])))
            else
                setfield!(L, field, read(g[String(field)]))
            end
        end
    end

    setneighbours!(L)
    setcheckerboards!(L)

    return L
end

function getgroup!(fid::HDF5.File, groupname::AbstractString)::HDF5.Group
    if ! haskey(fid, groupname)
        create_group(fid, groupname)
        return fid[groupname]
    elseif isempty(fid[groupname])
        return fid[groupname]
    else
        n = 2
        while haskey(fid, "$(groupname)_run$n")
            n += 1
        end
        create_group(fid, "$(groupname)_run$n")
        return fid["$(groupname)_run$n"]
    end
end


function dumpmetadata!(fid::HDF5.File, metadata::Dict{String,<:Any})
    g = getgroup!(fid, "metadata")
    for (key, value) in metadata
        g[key] = value
    end
end

function dumpMCparams!(fid::HDF5.File, MCParams::MCParameters)
    g = getgroup!(fid, "MCParameters")
    for param in fieldnames(MCParameters)
        if isa(getfield(MCParams, param), Symbol)
            g[String(param)] = String(getfield(MCParams, param))
        else
            g[String(param)] = getfield(MCParams, param)
        end
    end
end