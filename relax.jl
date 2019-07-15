function relax(P, pr=nothing, U=1e10)
    #m = Model(solver=IpoptSolver(print_level = 0))
    #m = Model(solver = GLPKSolverLP())
    #m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
    m = Model(solver=GurobiSolver(Method=2, Crossover=0, Threads=1, LogToConsole=0,OutputFlag=0,DualReductions = 0))#, LogFile=string(runName,".txt")))
    #m = Model(solver=GurobiSolver(Method=1, Threads=1, LogToConsole=0,OutputFlag=0))
    
    # Variables
    m.numCols = P.numCols
    m.colNames = P.colNames[:]
    m.colNamesIJulia = P.colNamesIJulia[:]
    m.colLower = P.colLower[:]
    m.colUpper = P.colUpper[:]
    m.colCat = P.colCat[:]
    m.colVal = P.colVal[:]

    # Constraints
    m.linconstr  = map(c->copy(c, m), P.linconstr)
    #m.quadconstr = map(c->copy(c, m), P.quadconstr)
        
    # Objective
    m.obj = copy(P.obj, m)
    m.objSense = P.objSense

    # extension
    if !isempty(P.ext)
        m.ext = similar(P.ext)
        for (key, val) in P.ext
            m.ext[key] = try
                copy(P.ext[key])
            catch
                continue;  #error("Error copying extension dictionary. Is `copy` defined for all your user types?")
            end
        end
    end

    #if U != 1e10
    @constraint(m, m.obj.aff <= U)
    #println("obj", m.obj.aff)
    #end

    qbcId = Dict()
    bilinearConsId = []
 
    for i = 1:length(P.quadconstr)
            con = P.quadconstr[i]
            terms = con.terms
            qvars1 = copy(terms.qvars1, m)
            qvars2 = copy(terms.qvars2, m)
            qcoeffs = terms.qcoeffs
            aff = copy(terms.aff, m)

            newB = []
	    definedB = []
	    definedBilinearVars = []
	    bilinearVars_local_con = []
	    #qbvarsId_local_con = []

	    for j in 1:length(qvars1)
	    	if haskey(qbcId, (qvars1[j].col, qvars2[j].col)) 
		    push!(definedB, j)
		    push!(definedBilinearVars,  Variable(m, qbcId[qvars1[j].col, qvars2[j].col][1]))
		elseif haskey(qbcId, (qvars2[j].col, qvars1[j].col))		
                    push!(definedB, j)
                    push!(definedBilinearVars,  Variable(m, qbcId[qvars2[j].col, qvars1[j].col][1]))
		else
		       xi = P.colVal[qvars1[j].col]
            	       yi = P.colVal[qvars2[j].col]
		       push!(newB, j)
		       bilinear = @variable(m)
		       setvalue(bilinear, xi*yi)
                       varname = " bilinear_con$(i)_"*string(qvars1[j])*"_"*string(qvars2[j])*"_$(j)"
                       setname(bilinear, varname)
		       push!(bilinearVars_local_con, bilinear)
		       qbcId[qvars1[j].col, qvars2[j].col] = (bilinear.col, length(m.linconstr)+1)
		       
		       xl=getlowerbound(qvars1[j])
                       xu=getupperbound(qvars1[j])
                       yl=getlowerbound(qvars2[j])
                       yu=getupperbound(qvars2[j])

		       m.colLower[end] = min(xl*yl, xl*yu, xu*yl, xu*yu)
		       m.colUpper[end] = max(xl*yl, xl*yu, xu*yl, xu*yu)	

		       if qvars1[j] == qvars2[j]
		       	   m.colLower[end] = max(m.colLower[end], 0)
		       end		       

		       if qvars1[j] == qvars2[j]
		       	   if (xu !=Inf)  && (yu != Inf)
		       	      temp = yu+xu
                              @constraint(m, bilinear >= temp*qvars1[j] - xu*yu)
			   end
			   if (xl !=-Inf) && (yl != -Inf)   
			       temp = yl+xl
                               @constraint(m, bilinear >= temp*qvars1[j] - xl*yl)
			   end
			   if (xl !=-Inf) && (yu != Inf)    
			       temp = yu+xl
                               @constraint(m, bilinear <= temp*qvars1[j] - xl*yu)
			   end    
		       else
			   if (xu !=Inf) && (yu != Inf)	
		               @constraint(m, bilinear >= yu*qvars1[j] + xu*qvars2[j] - xu*yu)
			   end
                           if (xl !=-Inf) && (yl != -Inf)	    
		               @constraint(m, bilinear >= yl*qvars1[j] + xl*qvars2[j] - xl*yl)
			   end  
                           if (xl !=-Inf) && (yu != Inf)	  
		               @constraint(m, bilinear <= yu*qvars1[j] + xl*qvars2[j] - xl*yu)
			   end  
                           if (xu !=Inf) && (yl != -Inf)	  
                               @constraint(m, bilinear <= yl*qvars1[j] + xu*qvars2[j] - xu*yl)
			   end    
                       end
                end
	    end
	    
            if con.sense == :(<=)
                    @constraint(m, aff + sum(qcoeffs[definedB[j]]*definedBilinearVars[j] for j in 1:length(definedB)) + sum(qcoeffs[newB[j]]*bilinearVars_local_con[j] for j in 1:length(newB)) <= 0)
            elseif con.sense == :(>=)
                    @constraint(m, aff + sum(qcoeffs[definedB[j]]*definedBilinearVars[j] for j in 1:length(definedB)) + sum(qcoeffs[newB[j]]*bilinearVars_local_con[j] for j in 1:length(newB)) >= 0)
            else
                    @constraint(m, aff + sum(qcoeffs[definedB[j]]*definedBilinearVars[j] for j in 1:length(definedB)) + sum(qcoeffs[newB[j]]*bilinearVars_local_con[j] for j in 1:length(newB)) == 0)
            end
    end
    println("additional  variables:  ",length(m.colLower)-length(P.colLower))
    m.ext[:qbcId] = qbcId
    if pr != nothing
       EqVconstr = pr.EqVconstr
       JuMP.addVectorizedConstraint(m,  map(c->copy(c, m), EqVconstr))
    end   
    return m
end


