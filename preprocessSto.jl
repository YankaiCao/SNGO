function preprocessSto!(P)
    scenarios = PlasmoOld.getchildren(P)	 
    # provide initial value if not defined
    for i in 1:length(P.colVal)
        if isnan(P.colVal[i])
           P.colVal[i] = 0
        end
    end

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
            P.colLower[i] = -1e8
        end
        if P.colUpper[i] == Inf
            P.colUpper[i] = 1e8
        end
    end
    for (idx,scenario) in enumerate(scenarios)
    	nsecond =  scenario.numCols
    	for i in 1:nsecond
            if scenario.colLower[i] == -Inf
                scenario.colLower[i] = -1e8
            end
            if scenario.colUpper[i] == Inf
                scenario.colUpper[i] = 1e8
            end
    	end
	println("first:  ", scenario.ext[:firstVarsId])
    end

    for (idx,scenario) in enumerate(scenarios)
    	firstVarsId = scenario.ext[:firstVarsId]
    	for i in 1:ncols	   
	    if firstVarsId[i] > 0
                if scenario.colLower[firstVarsId[i]] >  P.colLower[i] 
                    P.colLower[i] = scenario.colLower[firstVarsId[i]]
		end
                if scenario.colUpper[firstVarsId[i]] <  P.colUpper[i]
                    P.colUpper[i] = scenario.colUpper[firstVarsId[i]]
                end    
            end
        end    	
    end
    updateFirstBounds!(P, P.colLower, P.colUpper)
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