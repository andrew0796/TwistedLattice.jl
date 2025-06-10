using HDF5

""" Convert HDF5 files with SWMR possibly enabled to those without """

function copy_group!(dest::HDF5.Group, source::HDF5.Group)
    for obj_name in keys(source)
        obj = source[obj_name]
        if isa(obj, HDF5.Dataset)
            dest[obj_name] = read(obj)
        elseif isa(obj, HDF5.Group)
            create_group(dest, obj_name)
            new_dest = dest[obj_name]
            copy_group!(new_dest, obj)
        end
    end
end

function convert_file(filename::AbstractString, destination::AbstractString)
    f = h5open(filename, "r")
    f_new = h5open(destination, "w")

    for obj_name in keys(f)
        obj = f[obj_name]
        if isa(obj, HDF5.Group)
            create_group(f_new, obj_name)
            new_group = f_new[obj_name]
            copy_group!(new_group, obj)
        elseif isa(obj, HDF5.Dataset)
            f_new[obj_name] = read(obj)
        end
    end

    close(f)
    close(f_new)
end

function convert_file_with_suffix(filename::AbstractString, suffix::AbstractString="_non_SWMR", directory::AbstractString="")
    """ if directory is an relative path, the new files will be placed in a path relative to the original files """
    if isempty(suffix) && abspath(directory) == pwd()
        error("no suffix given and the directory given is the same as the original files")
    end

    abs_dir = ""
    f_name = ""
    f_dir, f_name = splitdir(abspath(filename))
    if !isabspath(directory)
        abs_dir = abspath(joinpath(f_dir, directory))
    else
        abs_dir = abspath(directory)
    end

    if !isdir(abspath(directory))
        mkdir(abspath(directory))
    end

    new_path = joinpath(abs_dir, replace(f_name, ".h5"=>suffix * ".h5"))

    convert_file(filename, new_path)
end


function run_convert()
    if length(ARGS) == 0
        println("convertfiles: copy HDF5 files potentially with SWMR enabled to those without, signature:")
        println("   convertfiles.jl [-s SUFFIX] [-d DIRECTORY] object1 object2 ... objectN")
        println("where (optional) SUFFIX is a suffix which is added to the filenames at the end (but before .h5)")
        println("and (optional) DIRECTORY is a directory to put the new files in. If DIRECTORY is relative, new files will be created in a directory relative to the original files.")
        println("If both SUFFIX and DIRECTORY are empty an error will be thrown.")
        println("The objects given are either HDF5 files or folders with HDF5 files.")
    end

    suffix = "_non_SWMR"
    directory = ""

    i = 1
    while i <= length(ARGS) && ARGS[i] in ["-s", "-d"]
        if ARGS[i] == "-s"
            if ARGS[i+1] == "-d" || isfile(ARGS[i+1]) || isdir(ARGS[i+1])
                suffix = ""
                i += 1
            else
                suffix = ARGS[i+1]
                i += 2
            end
        elseif ARGS[i] == "-d"
            if ARGS[i+1] == "-s" || isfile(ARGS[i+1])
                directory = ""
                i += 1
            else
                directory = ARGS[i+1]
                i += 2
            end
        end
    end

    while i <= length(ARGS)
        if isfile(ARGS[i]) && endswith(ARGS[i], ".h5")
            convert_file_with_suffix(ARGS[i], suffix, directory)
        elseif isdir(ARGS[i])
            for f in readdir(ARGS[i])
                if isfile(joinpath(ARGS[i], f)) && endswith(f, ".h5")
                    convert_file_with_suffix(joinpath(ARGS[i], f), suffix, directory)
                end
            end
        end
        i += 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_convert()
end