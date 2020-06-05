function preprocessSto!(P)
    scenarios = PlasmoOld.getchildren(P)	 
    # provide initial value if not defined
    for i in 1:length(P.colVal)
        if isnan(P.colVal[i])
           P.colVal[i] = 0
        end
    end

    # find index of first stage variables in each scenario. Linking constraints is stored in the master node. Each linking constraint involves 1 variable from master node and 1 varibale from a scenario node.
    ncols_first = P.numCols
    for (idx,scenario) in enumerate(scenarios)    
    	scenario.ext[:firstVarsId] = zeros(Int, ncols_first)
	scenario.ext[:firstVarsId][1:end]= -1
    end
    for c in 1:length(P.linconstr)
        coeffs = P.linconstr[c].terms.coeffs
        vars   = P.linconstr[c].terms.vars
	firstVarId = 0
        for (it,ind) in enumerate(coeffs)
            if (vars[it].m) == P
	        firstVarId = vars[it].col
		break
            end
        end
        for (it,ind) in enumerate(coeffs)
            if (vars[it].m) != P
               scenario = vars[it].m
	       scenario.ext[:firstVarsId][firstVarId] = vars[it].col 
            end
        end
    end

    # provide bounds if not defined
    ncols =  P.numCols
    for i in 1:ncols
        if P.colLower[i] == -Inf
            P.colLower[i] = default_lower_bound_value
        end
        if P.colUpper[i] == Inf
            P.colUpper[i] = default_upper_bound_value
        end
    end
    for (idx,scenario) in enumerate(scenarios)
    	for i in 1:scenario.numCols
            if scenario.colLower[i] == -Inf
                scenario.colLower[i] = default_lower_bound_value
            end
            if scenario.colUpper[i] == Inf
                scenario.colUpper[i] = default_upper_bound_value
            end
    	end
	if debug
	    println("first:  ", scenario.ext[:firstVarsId])
	end    
    end
    updateStoFirstBounds!(P)    


    nscen = length(scenarios)
    pr_children = []
    for (idx,scen) in enumerate(scenarios)
        #println("scen: ",scen)
        pr, scenarios[idx] = preprocess!(scen)
        push!(pr_children, pr)

        for i = 1:length(P.linconstr)
            Pcon = P.linconstr[i]
            Pvars = Pcon.terms.vars
            for (j, Pvar) in enumerate(Pvars)
                if (Pvar.m == scen)
                     Pvars[j] = Variable(scenarios[idx], Pvar.col)
                end
            end
        end
	scen = -1
    end
    return pr_children

end


function Stopreprocess!(P)
    pr_children = []
    scenarios = PlasmoOld.getchildren(P)
    for (idx,scenario) in enumerate(scenarios)
    	pr = preprocess!(scenario)
	push!(pr_children, pr)
    end
    return pr_children
end