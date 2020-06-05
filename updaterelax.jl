function updaterelax(R, P, pr, U=1e10, initialValue = nothing)
    m = copyModel(R)
    #m.solver = R.solver
    #m.solver = GurobiSolver(Method=1, Threads=1, LogToConsole=0,OutputFlag=0)
    m.solver = GurobiSolver(Method=2, Threads=1, Crossover=0, LogToConsole=0,OutputFlag=0, DualReductions = 0)    
    #m.solver = GurobiSolver(Method=2, Threads=1, Crossover=0, LogToConsole=0, OutputFlag=0, DualReductions = 0, BarHomogeneous=1)
    n = P.numCols	 
    m.colLower[1:n] = copy(P.colLower)
    m.colUpper[1:n] = copy(P.colUpper)
    #m.colVal = R.colVal[:]
    #m.colVal[1:n] = P.colVal[:]
    #m.colVal[n+1:end] = NaN
    if initialValue != nothing
       m.colVal[1:length(initialValue)] = copy(initialValue)
    end
    qbcId = m.ext[:qbcId]   #pr.qbcId 
    con = m.linconstr[length(P.linconstr)+1]    
    con.ub = U - m.obj.aff.constant

    for (key, value) in qbcId
    	    xid = key[1]
	    yid = key[2]
	    bid = value[1]
	    cid = value[2]
            xl = P.colLower[xid]
	    xu = P.colUpper[xid]
	    yl = P.colLower[yid]
	    yu = P.colUpper[yid]
	    
	    m.colLower[bid] = min(xl*yl, xl*yu, xu*yl, xu*yu)
	    m.colUpper[bid] = max(xl*yl, xl*yu, xu*yl, xu*yu)	
	    if xid == yid
	        m.colLower[bid] = max(m.colLower[bid], 0)
	    end
	    if initialValue !=nothing && length(initialValue) == n
	        m.colVal[bid] = m.colVal[xid]*m.colVal[yid]
	    end

	    if (xu !=Inf)  && (yu != Inf)
	        updateCon(m, m.linconstr[cid], - xu*yu, Inf, xid, -yu, yid, -xu)
	    end
	    if (xl !=-Inf) && (yl != -Inf)			
	        updateCon(m, m.linconstr[cid+1], - xl*yl, Inf, xid, -yl, yid, -xl)
	    end
	    if (xl !=-Inf) && (yu != Inf)			
	        updateCon(m, m.linconstr[cid+2], -Inf, -xl*yu, xid, -yu, yid, -xl)
	    end			
	    if xid != yid
	        if (xu !=Inf) && (yl != -Inf)
 	            updateCon(m, m.linconstr[cid+3], -Inf, -xu*yl, xid, -yl, yid, -xu)
		end    
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
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]

	    if (xu - xl)<= small_bound_improve 
	       	continue  
	    end	
            if op == :(<=) || op == :(==)
	        if (-1e8<= (exp(b*xu)-exp(b*xl))/(xu-xl) <= 1e8) && (-1e8 <= (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl)<=1e8)
	            if ev.cid[1] != -1 
	                updateCon(m, m.linconstr[ev.cid[1]], -Inf, (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl), nlvarId, - (exp(b*xu)-exp(b*xl))/(xu-xl) )
		    else		    
                        @constraint(m, Variable(m, lvarId) <= (exp(b*xu)-exp(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl))
		    end
		end
            end
            if op == :(>=) || op == :(==)
                if (-1e8<= b*exp(b*xl)<= 1e8) && (-1e8<= exp(b*xl)*(1-b*xl)<=1e8)
		    if ev.cid[2] != -1
		        updateCon(m, m.linconstr[ev.cid[2]], exp(b*xl)*(1-b*xl), Inf, nlvarId, -b*exp(b*xl))
		    else
                        @constraint(m, Variable(m, lvarId) >= b*exp(b*xl) * Variable(m, nlvarId) + exp(b*xl)*(1-b*xl))
		    end	
                end 
  	     
		if (-1e8<= b*exp(b*xu)<= 1e8) && (-1e8<= exp(b*xu)*(1-b*xu)<=1e8)
		    if ev.cid[3] != -1
		        updateCon(m, m.linconstr[ev.cid[3]], exp(b*xu)*(1-b*xu), Inf, nlvarId, -b*exp(b*xu))  
		    else
                        @constraint(m, Variable(m, lvarId) >= b*exp(b*xu) * Variable(m, nlvarId) + exp(b*xu)*(1-b*xu))
		    end	
		end
            end
    end


    #lvar op log(b * nlvar)
    for (i, ev) in enumerate(logVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
            xl = P.colLower[nlvarId]
            xu = P.colUpper[nlvarId]

	    #=
            if b >= 0 && xl <= 0
                xl = 1e-20
            elseif (b>= 0 && xu <= 0)  || (b< 0 && xl >= 0)
	    	   println("warning:  log(a), a cannot be negative ")
                return false
            elseif b< 0 && xu >= 0
	    	   xu = -1e-20
            end
	    =#

            if (xu - xl)<= small_bound_improve
                continue
            end
            if op == :(<=) || op == :(==)
                if (-1e8<= 1/xl<= 1e8) && (-1e8<= log(b*xl) <=1e8)
		    if ev.cid[1] != -1
		        updateCon(m, m.linconstr[ev.cid[1]], -Inf, log(b*xl) - 1,  nlvarId, -1/xl)
		    else
                        @constraint(m, Variable(m, lvarId) <= 1/xl * Variable(m, nlvarId) + log(b*xl) - 1 )
   		    end	
                end
                if (-1e8<= 1/xu<= 1e8) && (-1e8<= log(b*xu) <=1e8)
		    if ev.cid[2] != -1
			num_cons = MathProgBase.numconstr(m) 			
		        updateCon(m, m.linconstr[ev.cid[2]], -Inf, log(b*xu) - 1, nlvarId, -1/xu)
		    else
                        @constraint(m, Variable(m, lvarId) <= 1/xu * Variable(m, nlvarId) + log(b*xu) - 1 )
		    end	
                end
            end
            if op == :(>=) || op == :(==)
                if (-1e8<= (log(b*xu)-log(b*xl))/(xu-xl) <= 1e8) &&  (-1e8 <= (xu*log(b*xl)-xl*log(b*xu))/(xu-xl) <=1e8)
		    if ev.cid[3] != -1
		       updateCon(m, m.linconstr[ev.cid[3]], (xu*log(b*xl)-xl*log(b*xu))/(xu-xl), Inf, nlvarId, - (log(b*xu)-log(b*xl))/(xu-xl))
		    else
                        @constraint(m, Variable(m, lvarId) >= (log(b*xu)-log(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*log(b*xl)-xl*log(b*xu))/(xu-xl))
		    end	
                end
            end
    end

    for (i, ev) in enumerate(monomialVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
	    d = ev.d
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
                if (-1e8<= (exp(b*xu)-exp(b*xl))/(xu-xl) <= 1e8) && (-1e8 <= (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl)<=1e8)
                    if ev.cid[1] != -1
                        updateCon(m, m.linconstr[ev.cid[1]], -Inf, (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl), nlvarId, - (exp(b*xu)-exp(b*xl))/(xu-xl) )
                    else
                        @constraint(m, Variable(m, lvarId) <= (exp(b*xu)-exp(b*xl))/(xu-xl) * Variable(m, nlvarId) + (xu*exp(b*xl)-xl*exp(b*xu))/(xu-xl))
                    end
                end
            end
            if op == :(>=) || op == :(==)
                if (-1e8<= b*exp(b*xl)<= 1e8) && (-1e8<= exp(b*xl)*(1-b*xl)<=1e8)
                    if ev.cid[2] != -1
                        updateCon(m, m.linconstr[ev.cid[2]], exp(b*xl)*(1-b*xl), Inf, nlvarId, -b*exp(b*xl))
                    else
                        @constraint(m, Variable(m, lvarId) >= b*exp(b*xl) * Variable(m, nlvarId) + exp(b*xl)*(1-b*xl))
                    end
                end

                if (-1e8<= b*exp(b*xu)<= 1e8) && (-1e8<= exp(b*xu)*(1-b*xu)<=1e8)
                    if ev.cid[3] != -1
                        updateCon(m, m.linconstr[ev.cid[3]], exp(b*xu)*(1-b*xu), Inf, nlvarId, -b*exp(b*xu))
                    else
                        @constraint(m, Variable(m, lvarId) >= b*exp(b*xu) * Variable(m, nlvarId) + exp(b*xu)*(1-b*xu))
                    end
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

            if positiveFrac(d) || negativeFrac(d)
                if b >= 0
                    @assert xl >= 0
                else
                    @assert xu <= 0
                end
            end

            if op == :(<=) || op == :(==)
	        if positiveEven(d) || negativeFrac(d) || (negativeEven(d) && xu <= 0) || (negativeEven(d) && xl >= 0) || (Odd(d) && b>=0 && xl >= 0) || (Odd(d) && b<=0 && xu <= 0)
                    if -1e8 <= ((b*xu)^d-(b*xl)^d)/(xu-xl) <= 1e8 && -1e8<= (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl) <= 1e8
                        @constraint(m, Variable(m, lvarId) <= ((b*xu)^d-(b*xl)^d)/(xu-xl) * Variable(m, nlvarId) + (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl))
                    end
                elseif (positiveFrac(d)) || (Odd(d) && b>=0 && xu <= 0) || (Odd(d) && b<=0 && xl >= 0)
                    if (-1e8<= b*d*(b*xl)^(d-1)  <= 1e8) && (-1e8<= (- xl*b*d*(b*xl)^(d-1) + (b*xl)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) <= b*d*(b*xl)^(d-1) * Variable(m, nlvarId) - xl*b*d*(b*xl)^(d-1) + (b*xl)^d)
                    end
                    if (-1e8<= b*d*(b*xu)^(d-1)  <= 1e8) && (-1e8<= (- xu*b*d*(b*xu)^(d-1) + (b*xu)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) <= b*d*(b*xu)^(d-1) * Variable(m, nlvarId) - xu*b*d*(b*xu)^(d-1)     + (b*xu)^d)
                    end
                elseif positiveOdd(d) && xl < 0 && xu>0
                    #println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    #println("nothing to do")
                end
            end
            if op == :(>=) || op == :(==)
                if positiveFrac(d) || (Odd(d) && b>=0 && xu <= 0) || (Odd(d) && b<=0 && xl >= 0)
                    if -1e8 <= ((b*xu)^d-(b*xl)^d)/(xu-xl) <= 1e8 && -1e8<= (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl) <= 1e8
                        @constraint(m, Variable(m, lvarId) >= ((b*xu)^d-(b*xl)^d)/(xu-xl) * Variable(m, nlvarId) + (xu*(b*xl)^d-xl*(b*xu)^d)/(xu-xl))
                    end
                elseif positiveEven(d) || negativeFrac(d) || (negativeEven(d) && xu <= 0) || (negativeEven(d) && xl >= 0) || (Odd(d) && b>=0 && xl >= 0) || (Odd(d) && b<=0 && xu <= 0)
                    if (-1e8<= b*d*(b*xl)^(d-1)  <= 1e8) && (-1e8<= (- xl*b*d*(b*xl)^(d-1) + (b*xl)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) >= b*d*(b*xl)^(d-1) * Variable(m, nlvarId) - xl*b*d*(b*xl)^(d-1) + (b*xl)^d)
                    end
                    if (-1e8<= b*d*(b*xu)^(d-1)  <= 1e8) && (-1e8<= (- xu*b*d*(b*xu)^(d-1) + (b*xu)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) >= b*d*(b*xu)^(d-1) * Variable(m, nlvarId) - xu*b*d*(b*xu)^(d-1)     + (b*xu)^d)
                    end
                elseif positiveOdd(d)
                       #println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    #println("nothing to do")
                end
            end
    end
    
    return m
end


function updateCon(m, con, lb, ub, xid, xcoeff, yid = 0, ycoeff=0)
    con.lb = lb
    con.ub = ub
    if xid == yid
       xcoeff = xcoeff + ycoeff
    end       
    definedindex = find( x->(x.col == xid), con.terms.vars)
    if length(definedindex) > 0
        con.terms.coeffs[definedindex[1]] = xcoeff
    else
        push!(con.terms.coeffs, xcoeff)
        push!(con.terms.vars, Variable(m, xid))
    end

    if xid != yid && yid != 0
        definedindex = find( x->(x.col == yid), con.terms.vars)
	if length(definedindex) > 0
            con.terms.coeffs[definedindex[1]] = ycoeff
        else
	    push!(con.terms.coeffs, ycoeff)
            push!(con.terms.vars, Variable(m, yid))
        end
    end
    return con
end

function addaBB!(m, pr)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId
    multiVariable_aBB = pr.multiVariable_aBB
    for mv in multiVariable_aBB
                terms = mv.terms
		alpha = mv.alpha
                qvars1 = terms.qvars1
                qvars2 = terms.qvars2
                qcoeffs = terms.qcoeffs
		#println("add aBB")
                newcon = LinearConstraint(AffExpr(), 0, Inf)
                aff = newcon.terms

                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
		    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]

                    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)
                    #bid = qbcId[index[1]][3]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - coeff*yv)
                    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                for k in 1:length(mv.qVars)
                    xid = mv.qVars[k].col
                    xv = m.colVal[xid]
		    bid = qbcId[xid,xid][1]
		    bv = m.colVal[bid]
		    coeff = alpha[k]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - 2*coeff*xv)
                    constant += coeff*xv*xv
                    qconstant += coeff*bv
                end
                newcon.lb = - constant
                #if qconstant <= (constant - sigma_violation)
                   JuMP.addconstraint(m, newcon)
                #end
    end
    newncon = length(m.linconstr) - oldncon
    return newncon
end



function addOuterApproximationGrid!(m, pr, ngrid = 10)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId
    multiVariable_convex = pr.multiVariable_convex
    for (key, value) in qbcId
        xid = key[1]
        yid = key[2]
        bid = value[1]
        if xid == yid
            xl = m.colLower[xid]
            xu = m.colUpper[xid]
            for i = 1:ngrid
               xv = xl + (xu-xl)*i/(ngrid+1)
               @constraint(m, Variable(m, bid) - 2*xv*Variable(m, xid) + xv*xv >= 0)
            end
        end
    end
end


function addOuterApproximation!(m, pr)
    oldncon = length(m.linconstr)
    qbcId  = m.ext[:qbcId]   #pr.qbcId 
    #multiVariable_list = pr.multiVariable_list  
    multiVariable_convex = pr.multiVariable_convex

    for (key, value) in qbcId
        xid = key[1]
        yid = key[2]
        bid = value[1]
	if xid == yid
	    xv = m.colVal[xid]
	    bv = m.colVal[bid]    	    
	    if bv <= (xv^2 - sigma_violation)	    
	       @constraint(m, Variable(m, bid) - 2*xv*Variable(m, xid) + xv*xv >= 0)
	       #println("add constraint")
	    end
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
            xl = m.colLower[nlvarId]
            xu = m.colUpper[nlvarId]
	    xv = m.colVal[nlvarId]
            if op == :(>=) || op == :(==) 
	        if (xv-xl)>= machine_error && (xu-xv)>= machine_error && (-1e8<= b*exp(b*xv)<=1e8) && (-1e8<= exp(b*xv)*(1-b*xv) <= 1e8)		    
                    @constraint(m, Variable(m, lvarId) >= b*exp(b*xv) * Variable(m, nlvarId) + exp(b*xv)*(1-b*xv))
		end
            end
    end

    for (i, ev) in enumerate(logVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
            xl = m.colLower[nlvarId]
            xu = m.colUpper[nlvarId]
            xv = m.colVal[nlvarId]
	    if op == :(<=) || op == :(==)
                if (xv-xl)>= machine_error && (xu-xv)>= machine_error && (-1e8<= 1/xv<= 1e8) && (-1e8<= log(b*xv) <=1e8)
		    @constraint(m, Variable(m, lvarId) <= 1/xv * Variable(m, nlvarId) + log(b*xv) - 1 )
                end
            end
    end


    for (i, ev) in enumerate(monomialVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
	    d = ev.d
            xl = m.colLower[nlvarId]
            xu = m.colUpper[nlvarId]
            xv = m.colVal[nlvarId]

            if d <= 0
                error("warning:  a^x, a cannot be negative, please reformulate ")
            end
            b = b*log(d)

            if op == :(>=) || op == :(==)
                if (xv-xl)>= machine_error && (xu-xv)>= machine_error && (-1e8<= b*exp(b*xv)<=1e8) && (-1e8<= exp(b*xv)*(1-b*xv) <= 1e8)
                    @constraint(m, Variable(m, lvarId) >= b*exp(b*xv) * Variable(m, nlvarId) + exp(b*xv)*(1-b*xv))
                end
            end
    end


    for (i, ev) in enumerate(powerVariable_list)
            lvarId = ev.lvarId
            nlvarId = ev.nlvarId
            op = ev.op
            b = ev.b
	    d = ev.d
            xl = m.colLower[nlvarId]
            xu = m.colUpper[nlvarId]
            xv = m.colVal[nlvarId]

            if op == :(<=) || op == :(==)
                if positiveFrac(d) || (Odd(d) && b>=0 && xu <= 0) || (Odd(d) && b<=0 && xl >= 0)
                    if (-1e8<= b*d*(b*xv)^(d-1)  <= 1e8) && (-1e8<= (- xv*b*d*(b*xv)^(d-1) + (b*xv)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) <= b*d*(b*xv)^(d-1) * Variable(m, nlvarId) - xv*b*d*(b*xv)^(d-1) + (b*xv)^d)
                    end
                elseif positiveOdd(d) && xl < 0 && xu>0
                    #println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    #println("nothing to do")
                end
            end
            if op == :(>=) || op == :(==)
                if positiveEven(d) || negativeFrac(d) || (negativeEven(d) && xu <= 0) || (negativeEven(d) && xl >= 0) || (Odd(d) && b>=0 && xl >= 0) || (Odd(d) && b<=0 &&xu <= 0)
                    if (-1e8<= b*d*(b*xv)^(d-1)  <= 1e8) && (-1e8<= (- xv*b*d*(b*xv)^(d-1) + (b*xv)^d)  <=1e8)
                        @constraint(m, Variable(m, lvarId) >= b*d*(b*xv)^(d-1) * Variable(m, nlvarId) - xv*b*d*(b*xv)^(d-1) + (b*xv)^d)
                    end
                elseif positiveOdd(d)
                       #println("to do")
                elseif  (negativeEven(d) || negativeOdd(d)) && xl < 0 && xu>0
                    #println("nothing to do")
                end
            end
    end


    for mv in multiVariable_convex
	    if mv.pd == 1 && length(mv.qVarsId) > 1
	        terms = mv.terms
            	qvars1 = terms.qvars1
            	qvars2 = terms.qvars2
            	qcoeffs = terms.qcoeffs
	        #println("add convex")
	        newcon = LinearConstraint(AffExpr(), 0, Inf)
                aff = newcon.terms
							
		constant = 0
		qconstant = 0
	        for k in 1:length(qvars1)
		    xid = qvars1[k].col
                    yid = qvars2[k].col	  
		    xv = m.colVal[xid]
		    yv = m.colVal[yid]      
		    coeff = qcoeffs[k]
		    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)		    
		    #bid = qbcId[index[1]][3]
		    bid = mv.bilinearVars[k].col		    
		    bv = m.colVal[bid]

		    push!(aff.vars, Variable(m, bid))
		    push!(aff.coeffs, coeff)	       		    
		    push!(aff.vars, Variable(m,xid))
		    push!(aff.coeffs, - coeff*yv)
		    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)		  
		    constant += coeff*xv*yv
		    qconstant += coeff*bv
		end
		newcon.lb = - constant 

		if qconstant <= (constant - sigma_violation)
		   JuMP.addconstraint(m, newcon)
		end   
	    elseif mv.pd == -1 && length(mv.qVarsId) > 1
                terms = mv.terms
                qvars1 = terms.qvars1
                qvars2 = terms.qvars2
                qcoeffs = terms.qcoeffs
                #println("add convex")
                newcon = LinearConstraint(AffExpr(), -Inf, Inf)
                aff = newcon.terms

                constant = 0
                qconstant = 0
                for k in 1:length(qvars1)
                    xid = qvars1[k].col
                    yid = qvars2[k].col
                    xv = m.colVal[xid]
                    yv = m.colVal[yid]
                    coeff = qcoeffs[k]
                    #index = find( x->(x[1:2] == [xid, yid] || x[1:2] == [yid, xid]), qbcId)
                    #bid = qbcId[index[1]][3]
                    bid = mv.bilinearVars[k].col
                    bv = m.colVal[bid]

                    push!(aff.vars, Variable(m, bid))
                    push!(aff.coeffs, coeff)
                    push!(aff.vars, Variable(m,xid))
                    push!(aff.coeffs, - coeff*yv)
                    push!(aff.vars, Variable(m,yid))
                    push!(aff.coeffs, - coeff*xv)
                    constant += coeff*xv*yv
                    qconstant += coeff*bv
                end
                newcon.ub = - constant

                if qconstant >= (constant + sigma_violation)
                   JuMP.addconstraint(m, newcon)
                end
	    end
    end
    newncon = length(m.linconstr) - oldncon
    return newncon
end


