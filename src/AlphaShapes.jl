module AlphaShapes

    using Distributed

    @everywhere import BlackBoxOptim.bboptimize, BlackBoxOptim.best_candidate
    @everywhere import Distances.pairwise, Distances.Euclidean
    @everywhere import LinearAlgebra.det, LinearAlgebra.inv, LinearAlgebra.norm

    using MiniQhull

    export AlphaShape, AlphaShapeVolume

    const EQUAL_TOLERANCE = 6 #... decimal places, switching from Qhull to julia can lead to rounding errors

    #= e.g julia gave Qhull the point

    3-element Array{Float64,1}:
     1.2246467991473532e-16
     1.9999999999999998
     2.0000000000000004

    and this point in Qhulls triangulation became

    3-element Array{Float64,1}:
     -1.2246467991473532e-16
      2.0
      1.9999999999999998

    =#
    """
        CayleyMenger

    Compute the Cayley-Menger matrix from squared
    distances dij^2.

    https://mathworld.wolfram.com/Cayley-MengerDeterminant.html

    For a 2-d triangle the following example would symbolically
    represent the matrix for points x1, x2, and x3 where dij^2
    is the square Euclidean distance between i and j
    Example
```julia-repl
julia>:([0    1      1      1
         1    0      d12^2  d13^2
         1    d21^2  0      d23^2
         1    d32^2  d32^2    0])
:([0 1 1 1; 1 0 d12 ^ 2 d13 ^ 2; 1 d21 ^ 2 0 d23 ^ 2; 1 d32 ^ 2 d32 ^ 2 0])
```
    """
    function CayleyMenger(points::Array{Float64,2})::Array{Float64}
        d = pairwise(Euclidean(),points,dims=1)
        n = size(points,1)
        CM = ones(n+1,n+1)
        CM[2:end,2:end] = d.^2.0
        CM[end,end] = 0.0
        CM[1,1] = 0.0
        return CM
    end

    """
        SimplexVolume

    Calculate the volume of a simplex using the
    Cayley Menger matrix

    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexVolume(points::Array{Float64,2})::Float64
        try
            CM = CayleyMenger(points)
            n = size(CM,1)-2
            return sqrt(((-1.0)^(n+1))/((2.0^n)*factorial(n)^2.0)*det(CM))
        catch e
            return 0.0
        end
    end

    """
        SimplexCircumSphere

    Find the centre and radius of the circumsphere of a simplex
    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexCircumSphere(points::Array{Float64,2})::Tuple{Array{Float64,1},Float64}
        CM = CayleyMenger(points)
        cminv = inv(CM)
        R = sqrt(cminv[1,1]/(-2.0))
        # convert for barycentric to cartesian coordinates
        λ = cminv[1,2:end]
        c = (λ'*points)[1,:]
        return c,R
    end

    """
        SimplexCircumRadiusSquared

    Finds the squared cirum radius of a simplex
    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexCircumRadiusSquared(points::Array{Float64,2})::Float64
        CM = CayleyMenger(points)
        try
            cminv = inv(CM)
            return cminv[1,1]/(-2.0)
        catch e
            return Inf # Inf will be rejected by the alpha shape
        end
    end

    function VertexInSimplex(x::Array{Float64,1},T::Array{Float64,2})::Bool
        for i in 1:size(T,1)
            if x == T[i,:]
                return true
            end
        end
        return false
    end
    function VertexInTriangulation(x::Array{Float64,1},T::Array{Float64,3})::Bool
        for i in 1:size(T,1)
            if VertexInSimplex(x,T[i,:,:])
                return true
            end
        end
        return false
    end
    function AllPointsInAlphaShape(X::Array{Float64,2},A::Array{Float64,3},leak=1)::Bool
        x = round.(X,digits=EQUAL_TOLERANCE)
        a = round.(A,digits=EQUAL_TOLERANCE)
        n = 0
        for i in 1:size(X,1)
            if VertexInTriangulation(x[i,:],a) == false
                n += 1
                if n >= leak
                    return false
                end
            end
        end
        return true
    end

    """
        AlphaShapeVolume(A::Array{Float64,3})

    return the sum of volumes of all simplices in the alpha shapes A
    """
    function AlphaShapeVolume(A::Array{Float64,3})::Float64
        area = 0.0
        for t in 1:size(A,1)
            area += SimplexVolume(A[t,:,:])
        end
        return area
    end

    """
        GetDelaunayTriangulation(points::Array{Float64,2})::Array{Float64,3}

    Wrap MiniQhull.jl's delaunay to get a delaunay triangualation in any
    dimension
    """
    function GetDelaunayTriangulation(points::Array{Float64,2})::Array{Float64,3}
        tess = delaunay(permutedims(points,(2,1)))
        Triangles = zeros(size(tess,2),size(tess,1),size(tess,1)-1)
        for i in 1:size(tess,2)
            for j in 1:size(tess,1)
                Triangles[i,j,:] = points[tess[j,i],:]
            end
        end
        return Triangles
    end

    """
        FindAlpha(X::Array{Float64,2};
            search::Tuple{Float64,Float64}=(0.0, 10.0),
            MaxSteps::Int=100)::Float64

    Use BlackBocOptim to find the optimal alpha value (volume minimiser which cotains all
    points from the input).
    """
    function FindAlpha(X::Array{Float64,2};
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100,AllPoints::Bool=true,leak::Int64=0)::Float64
        objective = function(α)
            α = α[1]
            A = AlphaShape(X,α=α);
            t = AllPoints ? AllPointsInAlphaShape(X,A,leak) : true
            o =  AlphaShapeVolume(A)
            # minimise the area but if not all points are
            # included then set to extreme value
            if t
                return o
            else
                return Inf
            end
        end
        res = bboptimize(objective;
            SearchRange = search,
            NumDimensions = 1,
            MaxSteps=MaxSteps,
            Workers=workers());
        return best_candidate(res)[1]
    end

    """
        AlphaShape(X::Array{Float64,2};α::Union{Nothing,Float64}=nothing,
            search::Tuple{Float64,Float64}=(0.0, 10.0),
            MaxSteps::Int=100)::Array{Float64,3}

    Find the alpha shape corresponding to the 2D array of points
    X: [npoints,dimension], and the α value α.

    If α == nothing then a search over the range of values search is done
    to find the best value. The optimisation objective is the alpha shape
    area (minimise) subject to all points in X being included in the
    alpha shape triangulation.
    """
    function AlphaShape(X::Array{Float64,2};α::Union{Nothing,Float64}=nothing,
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100,AllPoints::Bool=true,leak::Int64=0)::Array{Float64,3}
        T = GetDelaunayTriangulation(X)
        if α == nothing
            println("Finding the optimum α value...\n")
            α = FindAlpha(X;search=search,MaxSteps=MaxSteps,AllPoints=AllPoints,leak=leak)
        end
        α2 = α^2.0
        A = [SimplexCircumRadiusSquared(T[i,:,:]) < α2 for i in 1:size(T,1)]
        return T[A.==1,:,:]
    end

    include("utils.jl")

end # module AlphaShapes
