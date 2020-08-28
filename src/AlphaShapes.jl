module AlphaShapes
    import BlackBoxOptim.bboptimize, BlackBoxOptim.best_candidate
    import Distances.pairwise, Distances.Euclidean
    import LinearAlgebra.det, LinearAlgebra.inv

    using MiniQhull

    export AlphaShape, AlphaShapeArea

    function CayleyMenger(points::Array{Float64,2})::Array{Float64}
        """
            CayleyMenger

        Compute the Cayley-Menger matrix from squared
        distances dij^2.

        e.g for a 2-d triangle

         0    1      1      1
         1    0      d12^2  d13^2
         1    d21^2  0      d23^2
         1    d32^2  d32^2    0
        """
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

        CM = CayleyMenger(points)
        n = size(CM,1)-2
        return sqrt(((-1.0)^(n+1))/((2.0^n)*factorial(n)^2.0)*det(CM))
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
        SimplexCircumSphere

    Find the centre and radius of the circumsphere of a simplex
    https://westy31.home.xs4all.nl/Circumsphere/ncircumsphere.htm
    """
    function SimplexCircumRadiusSquared(points::Array{Float64,2})::Float64
        CM = CayleyMenger(points)
        cminv = inv(CM)
        return cminv[1,1]/(-2.0)
    end

    function VertexInTriangle(x::Array{Float64,1},T::Array{Float64,2})::Bool
        if x == T[1,:] || x == T[2,:] || x == T[3,:]
            return true
        else
            return false
        end
    end
    function VertexInTriangulation(x::Array{Float64,1},T::Array{Float64,3})::Bool
        for i in 1:size(T,1)
            if VertexInTriangle(x,T[i,:,:])
                return true
            end
        end
        return false
    end
    function AllPointsInAlphaShape(X::Array{Float64,2},A::Array{Float64,3})::Bool
        for i in 1:size(X,1)
            if VertexInTriangulation(X[i,:],A) == false
                return false
            end
        end
        return true
    end

    function AlphaShapeVolume(A::Array{Float64,3})::Float64
        area = 0.0
        for t in 1:size(A,1)
            area += SimplexVolume(A[t,:,:])
        end
        return area
    end

    function GetTriangulation(points::Array{Float64,2})::Array{Float64,3}
        tess = delaunay(permutedims(points,(2,1)))
        Triangles = zeros(size(tess,2),size(tess,1),size(tess,1)-1)
        for i in 1:size(tess,2)
            for j in 1:size(tess,1)
                Triangles[i,j,:] = points[tess[j,i],:]
            end
        end
        return Triangles
    end

    function FindAlpha(X::Array{Float64,2};
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100)::Float64
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
        res = bboptimize(objective; SearchRange = search, NumDimensions = 1,MaxSteps=MaxSteps);
        return best_candidate(res)[1]
    end

    function AlphaShape(X::Array{Float64,2};α::Union{Nothing,Float64}=nothing,
        search::Tuple{Float64,Float64}=(0.0, 10.0),
        MaxSteps::Int=100)::Array{Float64,3}
        """
        Find the alpha shape corresponding to the 2D array of points
        X: [npoints,2], and the α value α.

        If α == nothing then a search over the range of values search is done
        to find the best value. The optimisation objective is the alpha shape
        area (minimise) subject to all points in X being included in the
        alpha shape triangulation.
        """
        T = GetTriangulation(X)
        if α == nothing
            println("Finding the optimum α value...\n")
            α = FindAlpha(X;search=search,MaxSteps=MaxSteps)
        end
        α2 = α^2.0
        A = [SimplexCircumRadiusSquared(T[i,:,:]) < α2 for i in 1:size(T,1)]
        return T[A.==1,:,:]
    end
end # module AlphaShapes
