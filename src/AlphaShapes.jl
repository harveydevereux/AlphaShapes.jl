module AlphaShapes
    import BlackBoxOptim.bboptimize, BlackBoxOptim.best_candidate
    import Distances.pairwise, Distances.Euclidean
    import LinearAlgebra.det, LinearAlgebra.inv

    using MiniQhull

    export AlphaShape, AlphaShapeVolume

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
    function CayleyMenger(points::AbstractArray{Float64,2})::AbstractArray{Float64}
        d = pairwise(Euclidean(),points,dims=2)
        n = size(points,2)
        CM = ones(n+1,n+1)
        view(CM, 2:n + 1,2:n + 1) .= d.^2.0
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
    function SimplexVolume(points::AbstractArray{Float64,2})::Float64

        CM = CayleyMenger(points)
        n = size(CM,2)-2
        return sqrt(((-1.0)^(n+1))/((2.0^n)*factorial(n)^2.0)*det(CM))
    end

    """
        SimplexCircumSphere

    Find the centre and radius of the circumsphere of a simplex
    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexCircumSphere(points::AbstractArray{Float64,2})::Tuple{AbstractArray{Float64,1},Float64}
        CM = CayleyMenger(points)
        cminv = inv(CM)
        R = sqrt(cminv[1,1]/(-2.0))
        # convert for barycentric to cartesian coordinates
        λ = cminv[2:end, 1]
        c = (points * λ)[:, 1]
        return c,R
    end

    """
        SimplexCircumSphere

    Find the centre and radius of the circumsphere of a simplex
    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexCircumRadiusSquared(points::AbstractArray{Float64,2})::Float64
        CM = CayleyMenger(points)
        cminv = inv(CM)
        return cminv[1,1]/(-2.0)
    end

    function VertexInTriangle(x::AbstractArray{Float64,1},T::AbstractArray{Float64,2})::Bool
        if x == T[:,1] || x == T[:,2] || x == T[:,3]
            return true
        else
            return false
        end
    end
    function VertexInTriangulation(x::AbstractArray{Float64,1},T::AbstractArray{Float64,3})::Bool
        for i in 1:size(T,3)
            if VertexInTriangle(x,view(T, :,:,i))
                return true
            end
        end
        return false
    end
    function AllPointsInAlphaShape(X::AbstractArray{Float64,2},A::AbstractArray{Float64,3})::Bool
        for i in 1:size(X,2)
            if VertexInTriangulation(view(X, :,i),A) == false
                return false
            end
        end
        return true
    end

    """
        AlphaShapeVolume(A::AbstractArray{Float64,3})

    return the sum of volumes of all simplices in the alpha shapes A
    """
    function AlphaShapeVolume(A::AbstractArray{Float64,3})::Float64
        area = 0.0
        for t in 1:size(A,3)
            area += SimplexVolume(view(A, :,:,t))
        end
        return area
    end

    """
        GetDelaunayTriangulation(points::AbstractArray{Float64,2})::AbstractArray{Float64,3}

    Wrap MiniQhull.jl's delaunay to get a delaunay triangualation in any
    dimension
    """
    function GetDelaunayTriangulation(points::AbstractArray{Float64,2})::AbstractArray{Float64,3}
        tess = delaunay(points)
        Triangles = zeros(size(tess,1)-1,size(tess,1),size(tess,2))
        for i in axes(tess, 1)
            for j in axes(tess, 2)
                view(Triangles, :, i, j) .= view(points, :, tess[i, j])
            end
        end
        return Triangles
    end

    """
        FindAlpha(X::AbstractArray{Float64,2};
            search::Tuple{Float64,Float64}=(0.0, 10.0),
            MaxSteps::Int=100)::Float64

    Use BlackBocOptim to find the optimal alpha value (volume minimiser which cotains all
    points from the input).
    """
    function FindAlpha(X::AbstractArray{Float64,2};
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100, silent::Bool = false)::Float64
        objective = function(α)
            α = α[1]
            A = AlphaShape(X,α=α);
            t = AllPointsInAlphaShape(X,A)
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
            TraceMode = silent ? :silent : :compact
        );
        return best_candidate(res)[1]
    end

    """
        AlphaShape(X::AbstractArray{Float64,2};α::Union{Nothing,Float64}=nothing,
            search::Tuple{Float64,Float64}=(0.0, 10.0),
            MaxSteps::Int=100)::AbstractArray{Float64,3}

    Find the alpha shape corresponding to the 2D AbstractArray of points
    X: [2,npoints], and the α value α.

    If α == nothing then a search over the range of values search is done
    to find the best value. The optimisation objective is the alpha shape
    area (minimise) subject to all points in X being included in the
    alpha shape triangulation.
    """
    function AlphaShape(X::AbstractArray{Float64,2};α::Union{Nothing,Float64}=nothing,
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100,
        silent::Bool=false)::AbstractArray{Float64,3}
        T = GetDelaunayTriangulation(X)
        if α == nothing
            if (!silent) println("Finding the optimum α value...\n") end
            α = FindAlpha(X;search=search,MaxSteps=MaxSteps,silent=silent)
        end
        α2 = α^2.0
        A = falses(size(T, 3))
        for i ∈ axes(T, 3)
            if SimplexCircumRadiusSquared(T[:, :, i]) < α2
                A[i] = true
            end
        end
        return T[:,:,A]
    end
end # module AlphaShapes
