module AlphaShapes
    import Triangle.basic_triangulation_vertices
    import BlackBoxOptim.bboptimize
    import BlackBoxOptim.best_candidate
    import LinearAlgebra.norm

    export AlphaShape, basic_triangulation_vertices, AlphaShapeArea

    function CircumCircle(trig_points)
        """
        Fairly efficiently compute the circumcircle for a given
        triangle, vertices are ordered as [vertex number, dimension]
        in 2D space
        """
        Bp = trig_points[2,:]-trig_points[1,:]
        Cp = trig_points[3,:]-trig_points[1,:]
        Dp = 2.0*(Bp[1]*Cp[2]-Bp[2]*Cp[1])

        B2 = Bp[1]^2. + Bp[2]^2.
        C2 = Cp[1]^2. + Cp[2]^2.
        Ux = (Cp[2]*(B2)-Bp[2]*(C2)) / Dp
        Uy = (Bp[1]*(C2)-Cp[1]*(B2)) / Dp
        return [Ux,Uy] .+ trig_points[1,:], norm([Ux,Uy])
    end

    function VertexInTriangle(x,T)
        if x == T[1,:] || x == T[2,:] || x == T[3,:]
            return true
        else
            return false
        end
    end
    function VertexInTriangulation(x,T)
        for i in 1:size(T,1)
            if VertexInTriangle(x,T[i])
                return true
            end
        end
        return false
    end
    function AllPointsInAlphaShape(X,A)
        for i in 1:size(X,1)
            if VertexInTriangulation(X[i,:],A) == false
                return false
            end
        end
        return true
    end
    function TriangleArea(T)
        a = norm(T[2,:]-T[1,:])
        b = norm(T[3,:]-T[2,:])
        c = norm(T[1,:]-T[3,:])
        s = (a+b+c)/2.
        return sqrt(s*(s-a)*(s-b)*(s-c))
    end
    function AlphaShapeArea(A)
        Area = 0.0
        for i in 1:size(A,1)
            Area += TriangleArea(A[i])
        end
        return Area
    end

    function FindAlpha(X;search=(0.0, 10.0),MaxSteps=100)
        objective = function(α)
            α = α[1]
            A = AlphaShape(X,α=α);
            t = AllPointsInAlphaShape(X,A)
            o =  AlphaShapeArea(A)
            # minimise the are but if not all points are
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

    function AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100)
        """
        Find the alpha shape corresponding to the 2D array of points
        X: [npoints,2], and the α value α.

        If α == nothing then a search over the range of values search is done
        to find the best value. The optimisation objective is the alpha shape
        area (minimise) subject to all points in X being included in the
        alpha shape triangulation.
        """
        T = basic_triangulation_vertices(X)
        if α == nothing
            println("Finding the optimum α value...\n")
            α = FindAlpha(X;search=search,MaxSteps=MaxSteps)
        end
        A = [CircumCircle(T[i])[2] < α for i in 1:size(T,1)]
        return T[A.==1]
    end
end # module AlphaShapes
