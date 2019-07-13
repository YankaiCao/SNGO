P = RandomStochasticModel(createModel, NS)   #RandomStochasticModel(createModel, nscen=10, nfirst=5, nparam=5)

# if solve with SNGO
m = copyStoModel(P)
branch_bound(m)




# if solve with SCIP
#=
m= extensiveSimplifiedModel(P)
m.solver = SCIPSolver("limits/gap", mingap, "limits/absgap", mingap, "limits/time", 43200.0)
solve(m)
=#


# if solve with Ipopt
#=
m= extensiveSimplifiedModel(P)
m = copyNLModel(m)
m.solver = IpoptSolver()
solve(m)
=#
