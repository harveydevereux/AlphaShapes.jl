using LinearAlgebra, Random, AlphaShapes, Plots, Statistics, DelimitedFiles, ProgressMeter

Random.seed!(31415926535897)

r1 = 0.75
r2 = 2.0

# bake 

X = randn(2,100)

ind = [((norm(X[:,i])>r1) + (norm(X[:,i])<r2)) == 2 for i in axes(X,2)]

X = X[:,ind.==true]

time, mem = (@timed AlphaShape(X))[[:time,:bytes]]

N = 2 .^ collect(3:10)
reps = 10
result = zeros(N|>length,reps,2) 

@showprogress for (k,n) in enumerate(N)
    
    for r in 1:reps

        X = randn(2,n)

        ind = [((norm(X[:,i])>r1) + (norm(X[:,i])<r2)) == 2 for i in axes(X,2)]
        
        X = X[:,ind.==true]

        time, mem = (@timed AlphaShape(X,MaxSteps=1,silent=true))[[:time,:bytes]]

        result[k,r,:] = [time, mem]

    end

end

mu = mean(result,dims=2)[:,1,:]
sig = std(result,dims=2)[:,1,:]

p = scatter(N, mu[:,1],label="runtime")
p = yerror!(N, mu[:,1], yerror=sig[:,1])
savefig(p,"benchmark-time.png")

p = scatter(N, mu[:,2],label="allocations")
p = yerror!(N, mu[:,2], yerror=sig[:,2])
savefig(p,"benchmark-mem.png")

data = zeros(N|>length,1+reps,2)
data[:,2:end,:] = result
data[:,1,1] = N
data[:,1,2] = N

writedlm("benchmark-time",data[:,:,1])
writedlm("benchmark-mem",data[:,:,2])
