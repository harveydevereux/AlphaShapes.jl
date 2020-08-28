function DrawTri!(T)
        plot!([T[1,1],T[2,1]],[T[1,2],T[2,2]],label="")
        plot!([T[2,1],T[3,1]],[T[2,2],T[3,2]],label="")
        plot!([T[3,1],T[1,1]],[T[3,2],T[1,2]],label="")
end

function Ellipse(n;a=1,b=1,max=2π,θ=0,center=[0.0,0.0])
    x = [a*cos(t)*cos(θ) - b*sin(t)*sin(θ) for t in 0:(2π/n):2π] .+ center[1]
    y = [a*cos(t)*sin(θ) + b*sin(t)*cos(θ) for t in 0:(2π/n):2π] .+ center[2]
    return cat(x,y,dims=2)
end

function PackingFraction(X,ParticleArea)
    return (size(X,2)*ParticleArea) ./ Alphashapes.AlphaShapeVolume(AlphaShape(X,FindAlpha(X)))
end
