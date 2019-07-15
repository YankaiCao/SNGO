using JuMP, SCIP
function test(iter)
for i = 1:iter
m= Model()
nvar = 1000
@variable(m, 0<=x[1:nvar]<=1)
@constraint(m, con[i in 1:(nvar-2)], x[i]+x[i+1]<=1)
@objective(m, :Min, sum{(x[i]+x[i+2])^2, i in 1:(nvar-2)})
m.solver=SCIPSolver("display/verblevel", 0)
solve(m)
end
end
test(10)