# Archived from the original repository. This experimental script is not part
# of the BilinearProjections.jl package API.

const datanodes = 256
const codenodes = 64
const batchsize = 100
const NEWTONITER = 10
const ERRFILE = "error.txt"

x = zeros(batchsize, datanodes)
z = zeros(batchsize, codenodes)
wp2 = zeros(codenodes, datanodes)
zp1 = zeros(batchsize, codenodes)
zp2 = zeros(batchsize, datanodes)
wref = zeros(codenodes)
zref = zeros(codenodes)

data = nothing
errptr = nothing
batchcount = 0
netiter = 0

function opendata(datafile)
    global data
    data = open(datafile, "r")
    if data === nothing
        println("Failed to open data file.")
        return false
    end
    return true
end

closedata() = close(data)

function readdata(count)
    for i in 1:datanodes
        value = readline(data)
        if value === nothing
            return false
        end
        x[count, i] = parse(Float64, value)
    end
    return true
end

function setup()
    global batchcount, netiter
    batchcount = 0
    netiter = 0
end

function encode(k)
    for i in 1:codenodes
        zp1[k, i] = tanh(sum(x[k, j] * wp2[i, j] for j in 1:datanodes))
    end
end

function updateweight()
    gram = zeros(codenodes, codenodes)
    gram2 = zeros(codenodes)
    for k in 1:batchsize
        gram += zp1[k, :] * zp1[k, :]'
        gram2 .+= zp1[k, :] * x[k, :]
    end
    for i in 1:codenodes
        nrm = dot(gram[i, :], gram2)
        if nrm != 0.0
            wn = wref[i] - nrm
            wn /= gram[i, i]
            wref[i] = wn
        end
    end
    for i in 1:codenodes, j in 1:datanodes
        wp2[i, j] += wref[i] * (x[batchcount, j] - zp2[batchcount, j])
    end
end

function initparam()
    fill!(wref, 0.0)
    fill!(wp2, 0.0)
end

function updateencoder()
    for i in 1:codenodes
        zref[i] = sum(zp1[j, i] for j in 1:batchsize) / batchsize
    end
end

function updateparam()
    for i in 1:codenodes, j in 1:datanodes
        wp2[i, j] -= zref[i] * x[batchcount, j]
    end
end

function calcerr(k)
    err = sum((x[k, i] - zp2[k, i])^2 for i in 1:datanodes)
    return sqrt(err / datanodes)
end

function decode(k)
    for i in 1:datanodes
        zp2[k, i] = dot(z[k, :], wp2[:, i])
    end
end

writeerror(err) = println(errptr, "%.15f", err)

function train(datafile)
    global errptr
    errptr = open(ERRFILE, "w")
    if errptr === nothing
        println("Failed to open error file.")
        return
    end
    setup()
    if !opendata(datafile)
        close(errptr)
        return
    end
    initparam()
    while true
        full_batch = true
        for k in 1:batchsize
            if !readdata(k)
                full_batch = false
                break
            end
        end
        full_batch || break
        for k in 1:batchsize
            encode(k)
            decode(k)
        end
        updateweight()
        updateencoder()
        updateparam()
        for k in 1:batchsize
            writeerror(calcerr(k))
        end
        global batchcount += 1
        batchcount % 1000 == 0 && println("Completed ", batchcount, " batches.")
    end
    closedata()
    close(errptr)
    println("Training completed.")
end
