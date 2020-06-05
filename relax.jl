function relax(P, pr=nothing, U=1e10)
    #m = Model(solver=IpoptSolver(print_level = 0))
    #m = Model(solver = GLPKSolverLP())
    #m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
    m = Model(solver=GurobiSolver(Method=2, Crossover=0, Threads=1, LogToConsole=0,OutputFlag=0,DualReductions = 0))
    #m = Model(solver = GurobiSolver(Method=2, Threads=1, Crossover=0, LogToConsole=0, OutputFlag=0, DualReductions = 0, BarHomogeneous=1))
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
		       #println("before setvalue   "," bilinear  ", bilinear.col, "    ",length(m.colVal) )
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


    expVariable_list = pr.expVariable_list
    logVariable_list = pr.logVariable_list
    powerVariable_list = pr.powerVariable_list
    monomialVariable_list = pr.monomialVariable_list

    for (i, ev) in enumerate(expVariable_list)
    	    lvarId = ev.lvarId
    	    nlvarId = ev.nlvarId
	    op = ev.op
	    b = ev.b
	    ev.cid = [-1;-1;-1]
    	    xl = P.colLower[nlvarId]
	    xu = P.colUpper[nlvarId]

            if (xu - xl)<= small_bound_improve
                continue
	    end 
	       
            if op == :(<=) || op == :(==)
	        if (-1e8<= (exp(b*xu)-exp(b*xl))/(xu-xl) <= 1e8) &&  (-1e8 <= (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl)<=1e8) && bounded(xl) && bounded(xu)
	            @constraint(m, Variable(m, lvarId) <= (exp(b*xu)-exp(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl))
		    ev.cid[1] = length(m.linconstr)
		end    
            end 		
            if op == :(>=) || op == :(==)
	        if (-1e8<= b*exp(b*xl)<= 1e8) && (-1e8<= exp(b*xl)*(1-b*xl)<=1e8) && bounded(xl)
                    @constraint(m, Variable(m, lvarId) >= b*exp(b*xl) * Variable(m, nlvarId) + exp(b*xl)*(1-b*xl))
		    ev.cid[2] = length(m.linconstr)
		end    
		if (-1e8<= b*exp(b*xu)<= 1e8) && (-1e8<= exp(b*xu)*(1-b*xu)<=1e8)&& bounded(xu)
		    @constraint(m, Variable(m, lvarId) >= b*exp(b*xu) * Variable(m, nlvarId) + exp(b*xu)*(1-b*xu))
		    ev.cid[3] = length(m.linconstr)
		end
	    end 
    end


    for (i, ev) in enumerate(monomialVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
            d = ev.d
	    ev.cid = [-1;-1;-1]
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]

            if (xu - xl)<= small_bound_improve
                continue
            end

            if d <= 0
                error("warning:  a^x, a cannot be negative, please reformulate ")
            end
            b = b*log(d)

            if op == :(<=) || op == :(==)
                if (-1e8<= (exp(b*xu)-exp(b*xl))/(xu-xl) <= 1e8) &&  (-1e8 <= (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl)<=1e8) && bounded(xl) && bounded(xu)
                    @constraint(m, Variable(m, lvarId) <= (exp(b*xu)-exp(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl))
                    ev.cid[1] = length(m.linconstr)
                end
            end
            if op == :(>=) || op == :(==)
                if (-1e8<= b*exp(b*xl)<= 1e8) && (-1e8<= exp(b*xl)*(1-b*xl)<=1e8) && bounded(xl)
                    @constraint(m, Variable(m, lvarId) >= b*exp(b*xl) * Variable(m, nlvarId) + exp(b*xl)*(1-b*xl))
                    ev.cid[2] = length(m.linconstr)
                end
                if (-1e8<= b*exp(b*xu)<= 1e8) && (-1e8<= exp(b*xu)*(1-b*xu)<=1e8)&& bounded(xu)
                    @constraint(m, Variable(m, lvarId) >= b*exp(b*xu) * Variable(m, nlvarId) + exp(b*xu)*(1-b*xu))
                    ev.cid[3] = length(m.linconstr)
                end
            end
    end

    #lvar op log(b * nlvar)
    for (i, ev) in enumerate(logVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
	    b = ev.b
            ev.cid = [-1;-1;-1]
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]

            if (xu - xl)<= small_bound_improve
                continue
            end

	    #=
            if b >= 0 && xl <= 0
	        error("lvarId:  ", lvarId, "   ",xl, "    ", xu )
                xl = 1e-20
		P.colLower[nlvarId] = xl
            elseif (b>= 0 && xu <= 0)  || (b< 0 && xl >= 0)
                println("warning:  log(a), a cannot be negative ")
                return false
            elseif b< 0 && xu >= 0
	    	error("lvarId:  ", lvarId, "   ",xl, "    ", xu )   
                xu = -1e-20
		P.colUpper[nlvarId] = xu
            end
	    =#

            if op == :(<=) || op == :(==)
                if (-1e8<= 1/xl<= 1e8) && (-1e8<= log(b*xl) <=1e8) && bounded(xl)
                    @constraint(m, Variable(m, lvarId) <= 1/xl * Variable(m, nlvarId) + log(b*xl) - 1 )
                    ev.cid[1] = length(m.linconstr)
                end
                if (-1e8<= 1/xu<= 1e8) && (-1e8<= log(b*xu) <=1e8) && bounded(xu)
                    @constraint(m, Variable(m, lvarId) <= 1/xu * Variable(m, nlvarId) + log(b*xu) - 1 )
                    ev.cid[2] = length(m.linconstr)

		    num_cons = MathProgBase.numconstr(m)
		    #println("relax:    ", ev.cid, "totoal  ",num_cons, " linear  ", length(m.linconstr))
                end
            end
            if op == :(>=) || op == :(==) 
                if (-1e8<= (log(b*xu)-log(b*xl))/(xu-xl) <= 1e8) &&  (-1e8 <= (xu*log(b*xl)-xl*log(b*xu))/(xu-xl) <=1e8) && bounded(xl)&& bounded(xl)
                    @constraint(m, Variable(m, lvarId) >= (log(b*xu)-log(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*log(b*xl)-xl*log(b*xu))/(xu-xl))
                    ev.cid[3] = length(m.linconstr)
                end
            end
    end

    
    for (i, ev) in enumerate(powerVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
            d = ev.d
            ev.cid = [-1;-1;-1]
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]

            if (xu - xl)<= small_bound_improve
                continue
            end
	    
	    #println("power  variable:    ", ev)
	    #println("xl  xu", xl, "     ",xu)
            if positiveFrac(d) || negativeFrac(d)
                if b >= 0
		    xl = max(xl, 1e-20)
                else
                    xu = min(xu, -1e-20)
                end
            end


            if op == :(<=) || op == :(==)
                if positiveEven(d) || negativeFrac(d) || (negativeEven(d) && xu <= 0) || (negativeEven(d) && xl >= 0) || (Odd(d) && b>=0 && xl >= 0) || (Odd(d) && b<=0 && xu <= 0)
                    if -1e8 <= ((b*xu)^d-(b*xl)^d)/(xu-xl) <= 1e8 && -1e8<= (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl) <= 1e8  && bounded(xl)&& bounded(xu)
                        @constraint(m, Variable(m, lvarId) <= ((b*xu)^d-(b*xl)^d)/(xu-xl) * Variable(m, nlvarId) + (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl))
                        ev.cid[1] = length(m.linconstr)
                    end
                elseif (positiveFrac(d)) || (Odd(d) && b>=0 && xu <= 0) || (Odd(d) && b<=0 && xl >= 0) 
                    if (-1e8<= b*d*(b*xl)^(d-1)  <= 1e8) && (-1e8<= (- xl*b*d*(b*xl)^(d-1) + (b*xl)^d)  <=1e8) && bounded(xl)
                        @constraint(m, Variable(m, lvarId) <= b*d*(b*xl)^(d-1) * Variable(m, nlvarId) - xl*b*d*(b*xl)^(d-1) + (b*xl)^d)
                    end
                    if (-1e8<= b*d*(b*xu)^(d-1)  <= 1e8) && (-1e8<= (- xu*b*d*(b*xu)^(d-1) + (b*xu)^d)  <=1e8) && bounded(xu)
                        @constraint(m, Variable(m, lvarId) <= b*d*(b*xu)^(d-1) * Variable(m, nlvarId) - xu*b*d*(b*xu)^(d-1)     + (b*xu)^d)
                    end
                elseif positiveOdd(d) && xl < 0 && xu>0
                    println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    println("nothing to do")
                end
            end
	   

            if op == :(>=) || op == :(==)
                if positiveFrac(d) || (Odd(d) && b>=0 && xu <= 0) || (Odd(d) && b<=0 && xl >= 0) 
                    if -1e8 <= ((b*xu)^d-(b*xl)^d)/(xu-xl) <= 1e8 && -1e8<= (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl) <= 1e8 && bounded(xl) && bounded(xu)
                        @constraint(m, Variable(m, lvarId) >= ((b*xu)^d-(b*xl)^d)/(xu-xl) * Variable(m, nlvarId) + (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl))
			#println("m.linconstr[end]:     ", m.linconstr[end])
                        ev.cid[1] = length(m.linconstr)
                    end
                elseif positiveEven(d) || negativeFrac(d) || (negativeEven(d) && xu <= 0) || (negativeEven(d) && xl >= 0) || (Odd(d) && b>=0 && xl >= 0) || (Odd(d) && b<=0 && xu <= 0)
                    if (-1e8<= b*d*(b*xl)^(d-1)  <= 1e8) && (-1e8<= (- xl*b*d*(b*xl)^(d-1) + (b*xl)^d)  <=1e8)&& bounded(xl)
                        @constraint(m, Variable(m, lvarId) >= b*d*(b*xl)^(d-1) * Variable(m, nlvarId) - xl*b*d*(b*xl)^(d-1) + (b*xl)^d)
                    end
                    if (-1e8<= b*d*(b*xu)^(d-1)  <= 1e8) && (-1e8<= (- xu*b*d*(b*xu)^(d-1) + (b*xu)^d)  <=1e8)&& bounded(xu)
                        @constraint(m, Variable(m, lvarId) >= b*d*(b*xu)^(d-1) * Variable(m, nlvarId) - xu*b*d*(b*xu)^(d-1)     + (b*xu)^d)
                    end
                elseif positiveOdd(d)
                       println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    println("nothing to do")
                end
            end
    end
    #println(m)
    

    println("additional  variables:  ",length(m.colLower)-length(P.colLower))
    m.ext[:qbcId] = qbcId
    if pr != nothing
       EqVconstr = pr.EqVconstr
       JuMP.addVectorizedConstraint(m,  map(c->copy(c, m), EqVconstr))
    end   
    return m
end


