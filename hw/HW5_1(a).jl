function squareRoot(x0,A)
	f(x) = x^2 - A
	x = x0
	i = 0
	while i <= k
		x = 0.5*(x+A/x)
		i = i+1
	end
	sqrtA = x
	return sqrtA
end