function newtonSolver(f,M,x0,TOL)
        x1=x0[1,1]
        x2=x0[2,1]
	x=x0
        while norm(f(x1,x2))>TOL
        x=x-inv(M(x1,x2))*f(x1,x2)
	x1=x[1,1]
	x2=x[2,1]
        end
        return x
        end